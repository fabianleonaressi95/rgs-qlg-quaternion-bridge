"""Finite A5-root Cayley tori with their faithful S6 action.

Coordinates are in the simple-root basis alpha_k=e_k-e_(k+1), k=0,...,4.
All arithmetic is integral; finite quotients reduce coordinates modulo m.
"""
from __future__ import annotations

from itertools import permutations, product


def from_root_coordinates(q: tuple[int, ...]) -> tuple[int, ...]:
    """Convert simple-root coordinates in Z^5 to ambient A5 coordinates."""
    if len(q) != 5:
        raise ValueError("expected five simple-root coordinates")
    y = [q[0]]
    y.extend(q[i] - q[i - 1] for i in range(1, 5))
    y.append(-q[4])
    return tuple(y)


def to_root_coordinates(y: tuple[int, ...]) -> tuple[int, ...]:
    """Convert an ambient integer vector of coordinate sum zero to A5 basis."""
    if len(y) != 6 or sum(y) != 0:
        raise ValueError("expected an integer six-vector with coordinate sum zero")
    return tuple(sum(y[:k + 1]) for k in range(5))


def root_step(i: int, j: int) -> tuple[int, ...]:
    """Coordinates of root e_i-e_j in the simple-root basis."""
    if not (0 <= i < 6 and 0 <= j < 6 and i != j):
        raise ValueError("i,j must be distinct elements of range(6)")
    y = [0] * 6
    y[i] += 1
    y[j] -= 1
    return to_root_coordinates(tuple(y))


def act_permutation(q: tuple[int, ...], sigma: tuple[int, ...]) -> tuple[int, ...]:
    """Apply coordinate permutation e_i -> e_{sigma(i)} to an A5 lattice vector."""
    if len(sigma) != 6 or sorted(sigma) != list(range(6)):
        raise ValueError("sigma must be a permutation of range(6)")
    y = from_root_coordinates(q)
    permuted = [0] * 6
    for i, value in enumerate(y):
        permuted[sigma[i]] += value
    return to_root_coordinates(tuple(permuted))


def root_steps() -> tuple[tuple[int, ...], ...]:
    """The 30 roots ±(e_i-e_j), represented in the simple-root basis."""
    return tuple(root_step(i, j) for i in range(6) for j in range(6) if i != j)


def add_mod(a: tuple[int, ...], b: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    return tuple((x + y) % modulus for x, y in zip(a, b))


def cayley_neighbors(vertex: tuple[int, ...], modulus: int) -> tuple[tuple[int, ...], ...]:
    if modulus < 3:
        raise ValueError("modulus must be at least 3")
    return tuple(add_mod(vertex, r, modulus) for r in root_steps())


def finite_a5_audit(modulus: int = 3) -> dict:
    """Exhaustively audit root steps and all S6 actions on the root set."""
    if modulus < 3:
        raise ValueError("modulus must be at least 3")
    roots = root_steps()
    assert len(roots) == 30
    assert len(set(tuple(x % modulus for x in r) for r in roots)) == 30
    assert all(any(x % modulus for x in r) for r in roots)

    root_action_checks = 0
    faithful_action_checks = 0
    identity = tuple(range(6))
    root_classes = {tuple(x % modulus for x in r) for r in roots}
    for sigma in permutations(range(6)):
        image = {tuple(x % modulus for x in act_permutation(r, sigma)) for r in roots}
        assert image == root_classes
        root_action_checks += len(roots)
        if sigma != identity:
            assert any(
                tuple(x % modulus for x in act_permutation(r, sigma))
                != tuple(x % modulus for x in r)
                for r in roots
            )
            faithful_action_checks += 1

    vertices = list(product(range(modulus), repeat=5))
    degree_checks = 0
    for v in vertices:
        neighbors = cayley_neighbors(v, modulus)
        assert len(neighbors) == 30
        assert len(set(neighbors)) == 30
        degree_checks += 1

    # Connectedness: the five simple roots are among the root steps and form a Z-basis.
    return {
        "model": "finite A5-root Cayley torus",
        "modulus": modulus,
        "vertex_count": modulus ** 5,
        "root_directions": len(roots),
        "degree": 30,
        "undirected_edge_count": 15 * modulus ** 5,
        "group_action": "faithful S6 by graph automorphisms; includes A6",
        "S6_order": 720,
        "root_action_checks": root_action_checks,
        "nonidentity_permutations_detected": faithful_action_checks,
        "vertex_degree_checks": degree_checks,
        "connectedness": "proved: simple roots generate A5/mA5",
        "status": "exact finite combinatorial audit",
    }


if __name__ == "__main__":
    import json
    import sys
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print(json.dumps(finite_a5_audit(m), indent=2, sort_keys=True))
