"""RGS Network Dynamics — synthetic baseline prototype v0.1.

Compares three graph-generation rules using identical synthetic agents,
skill profiles, initial states, node count and edge count.

This is a toy model, not an empirically calibrated model of human learning,
organizational performance, or Sun Microsystems' Practices experience.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np

Edge = Tuple[int, int]
N_AGENTS = 30
N_SKILLS = 8
N_EDGES = 45
DEFAULT_SEED = 20261119
DEFAULT_STEPS = 160
DEFAULT_ETA = 0.22
DEFAULT_RHO = 0.008
DEFAULT_REMOVAL_FRACTION = 0.20
DEFAULT_ROBUSTNESS_REPETITIONS = 80


def make_skill_profiles(rng: np.random.Generator, n_agents: int = N_AGENTS,
                        n_skills: int = N_SKILLS) -> np.ndarray:
    """Give each synthetic agent two or three initially available skills."""
    profiles = np.zeros((n_agents, n_skills), dtype=float)
    for agent in range(n_agents):
        count = int(rng.integers(2, 4))
        profiles[agent, rng.choice(n_skills, size=count, replace=False)] = 1.0
    return profiles


def _tree_backbone(n_agents: int) -> Set[Edge]:
    return {(min(i, (i - 1) // 2), max(i, (i - 1) // 2))
            for i in range(1, n_agents)}


def _select_top_edges(candidates: Iterable[Edge], scores: Mapping[Edge, float],
                      target: int) -> Set[Edge]:
    ordered = sorted(candidates, key=lambda edge: (-scores[edge], edge[0], edge[1]))
    return set(ordered[:target])


def build_networks(profiles: np.ndarray, seed: int = DEFAULT_SEED,
                   n_edges: int = N_EDGES) -> Dict[str, Set[Edge]]:
    """Build equal-size hierarchical, matrix and Practices-like undirected graphs."""
    if profiles.ndim != 2 or profiles.shape[0] < 2:
        raise ValueError("profiles must be a 2D array with at least two agents")
    n_agents = int(profiles.shape[0])
    max_edges = n_agents * (n_agents - 1) // 2
    if not (n_agents - 1 <= n_edges <= max_edges):
        raise ValueError("n_edges must be between n_agents-1 and the simple-graph maximum")

    rng = np.random.default_rng(seed)
    pairs = [(i, j) for i in range(n_agents) for j in range(i + 1, n_agents)]

    # Hierarchical: binary-tree backbone plus links between similar tree depths.
    hierarchy = _tree_backbone(n_agents)
    depth = lambda i: int(np.floor(np.log2(i + 1)))
    candidates = [edge for edge in pairs if edge not in hierarchy]
    scores = {edge: 1.0 / (1.0 + abs(depth(edge[0]) - depth(edge[1])))
              + 0.02 * rng.random() for edge in candidates}
    hierarchy |= _select_top_edges(candidates, scores, n_edges - len(hierarchy))

    # Matrix: prioritize complementary skill profiles.
    matrix_scores: Dict[Edge, float] = {}
    for i, j in pairs:
        overlap = float(np.dot(profiles[i], profiles[j]))
        union = float(np.count_nonzero(profiles[i] + profiles[j]))
        similarity = overlap / union if union else 0.0
        matrix_scores[(i, j)] = 1.0 - similarity + 0.03 * rng.random()
    matrix = _select_top_edges(pairs, matrix_scores, n_edges)

    # Practices-like: prioritize shared competencies / practice membership.
    practices_scores = {
        (i, j): float(np.dot(profiles[i], profiles[j])) + 0.04 * rng.random()
        for i, j in pairs
    }
    practices = _select_top_edges(pairs, practices_scores, n_edges)
    return {"hierarchical": hierarchy, "matrix": matrix, "practices_like": practices}


def adjacency_matrix(n_agents: int, edges: Iterable[Edge]) -> np.ndarray:
    matrix = np.zeros((n_agents, n_agents), dtype=float)
    for i, j in edges:
        if i == j or min(i, j) < 0 or max(i, j) >= n_agents:
            raise ValueError(f"invalid edge {(i, j)}")
        matrix[i, j] = matrix[j, i] = 1.0
    return matrix


def neighbor_mixing_matrix(adjacency: np.ndarray) -> np.ndarray:
    """Return a row-stochastic matrix; isolated nodes retain their own state."""
    degree = adjacency.sum(axis=1)
    mixing = np.zeros_like(adjacency, dtype=float)
    for i, value in enumerate(degree):
        if value > 0:
            mixing[i] = adjacency[i] / value
        else:
            mixing[i, i] = 1.0
    return mixing


def simulate_knowledge(mixing: np.ndarray, initial: np.ndarray, profiles: np.ndarray,
                       steps: int = DEFAULT_STEPS, eta: float = DEFAULT_ETA,
                       rho: float = DEFAULT_RHO) -> np.ndarray:
    """Iterate the declared bounded knowledge-state update."""
    if not (0.0 <= eta <= 1.0 and 0.0 <= rho <= 1.0):
        raise ValueError("eta and rho must lie in [0, 1]")
    if steps < 0:
        raise ValueError("steps must be non-negative")
    state = np.asarray(initial, dtype=float).copy()
    target = np.asarray(profiles, dtype=float)
    if state.shape != target.shape or mixing.shape != (state.shape[0], state.shape[0]):
        raise ValueError("incompatible state, profile and mixing dimensions")
    for _ in range(steps):
        state = state + eta * (mixing @ state - state) + rho * (target - state)
        state = np.clip(state, 0.0, 1.0)
    return state


def connected_components(adjacency: np.ndarray) -> List[List[int]]:
    seen: Set[int] = set()
    components: List[List[int]] = []
    for start in range(adjacency.shape[0]):
        if start in seen:
            continue
        stack, component = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in np.flatnonzero(adjacency[node]):
                neighbor = int(neighbor)
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def robustness_score(adjacency: np.ndarray, rng: np.random.Generator,
                     removal_fraction: float = DEFAULT_REMOVAL_FRACTION,
                     repetitions: int = DEFAULT_ROBUSTNESS_REPETITIONS) -> float:
    """Mean largest-component fraction after uniform random node removal."""
    if not (0.0 <= removal_fraction < 1.0) or repetitions < 1:
        raise ValueError("invalid removal fraction or repetition count")
    n_agents = adjacency.shape[0]
    n_removed = int(round(removal_fraction * n_agents))
    scores: List[float] = []
    for _ in range(repetitions):
        removed = set(rng.choice(n_agents, size=n_removed, replace=False).tolist())
        retained = [i for i in range(n_agents) if i not in removed]
        if not retained:
            scores.append(0.0)
            continue
        reduced = adjacency[np.ix_(retained, retained)]
        largest = max(len(component) for component in connected_components(reduced))
        scores.append(largest / len(retained))
    return float(np.mean(scores))


def run_experiment(seed: int = DEFAULT_SEED, n_agents: int = N_AGENTS,
                   n_skills: int = N_SKILLS, n_edges: int = N_EDGES,
                   steps: int = DEFAULT_STEPS) -> List[dict]:
    rng = np.random.default_rng(seed)
    profiles = make_skill_profiles(rng, n_agents, n_skills)
    initial = profiles * rng.uniform(0.15, 0.55, size=profiles.shape)
    networks = build_networks(profiles, seed=seed + 1, n_edges=n_edges)
    rows: List[dict] = []
    for index, (name, edges) in enumerate(networks.items()):
        adjacency = adjacency_matrix(n_agents, edges)
        mixing = neighbor_mixing_matrix(adjacency)
        final_state = simulate_knowledge(mixing, initial, profiles, steps=steps)
        robustness_rng = np.random.default_rng(seed + 1000 + index)
        rows.append({
            "architecture": name,
            "nodes": n_agents,
            "skills": n_skills,
            "edges": len(edges),
            "mean_degree": float(adjacency.sum() / n_agents),
            "connected_components": len(connected_components(adjacency)),
            "skill_coverage_threshold_0_60": float(np.mean(final_state > 0.60)),
            "mean_knowledge_state": float(final_state.mean()),
            "within_skill_dispersion": float(np.mean(np.std(final_state, axis=0))),
            "robustness_largest_component_after_random_removal": robustness_score(
                adjacency, robustness_rng
            ),
        })
    return rows


def write_results(rows: Sequence[dict], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "summary.json").write_text(
        json.dumps(list(rows), indent=2) + "\n", encoding="utf-8"
    )
    if rows:
        with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--agents", type=int, default=N_AGENTS)
    parser.add_argument("--skills", type=int, default=N_SKILLS)
    parser.add_argument("--edges", type=int, default=N_EDGES)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    rows = run_experiment(args.seed, args.agents, args.skills, args.edges, args.steps)
    write_results(rows, args.output)
    if rows:
        print(",".join(rows[0].keys()))
        for row in rows:
            print(",".join(str(row[key]) for key in row))


if __name__ == "__main__":
    main()
