"""Non-degenerate integral face lift for the existing Kathara -> Q6 vertex chart.

Each source edge is sent to a shortest Q6 path. Each source face boundary is
filled by exact cubical square swaps. The fixed axis order is a declared
combinatorial convention selected by an exact finite rank audit; it is not
claimed to be RGS-canonical. Compatible Ellis/SU(2) transport is a separate
open gate.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from src.kathara_cellular_complex import audit_kathara_cellular_complex
from src.kathara_q6_chain_map import kathara_vertex_to_q6_label
from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data


# Fixed target-axis order. This is a declared path convention, not the
# canonical numerical axis order and not a metric-derived RGS selection.
AXIS_ORDER = (1, 2, 0, 3, 4, 5)
_AXIS_POSITION = {axis: i for i, axis in enumerate(AXIS_ORDER)}


@dataclass(frozen=True)
class NondegenerateKatharaQ6Audit:
    vertex_labels: tuple[int, ...]
    R0: np.ndarray
    R1: np.ndarray
    R2: np.ndarray
    chain_identity_0: bool
    chain_identity_1: bool
    R2_rank_mod_prime: int
    face_boundary_l1: tuple[int, ...]
    face_fill_l1: tuple[int, ...]
    mapped_rgs_coordinates: np.ndarray
    minimum_source_edge_length: float
    status: str


def _axis_mask(axis: int) -> int:
    return 1 << (5 - axis)


def _path_axes(start: int, end: int) -> list[int]:
    """Canonical shortest path axes in the declared AXIS_ORDER."""
    diff = start ^ end
    return [axis for axis in AXIS_ORDER if diff & _axis_mask(axis)]


def _oriented_edge_path_axes(u: int, v: int, labels: tuple[int, ...]) -> list[int]:
    """Path word for directed source edge u->v, matching R1's orientation.

    R1 stores columns for canonical source edges (min(u,v), max(u,v)).
    Traversing such an edge backwards reverses the entire path word.
    """
    lo, hi = sorted((u, v))
    axes = _path_axes(labels[lo - 1], labels[hi - 1])
    return axes if (u, v) == (lo, hi) else list(reversed(axes))


def _add_path(chain: np.ndarray, start: int, axes: list[int],
              edge_index: dict[tuple[int, int], int], sign: int = 1) -> int:
    current = start
    for axis in axes:
        nxt = current ^ _axis_mask(axis)
        edge = (min(current, nxt), max(current, nxt))
        chain[edge_index[edge]] += sign if current < nxt else -sign
        current = nxt
    return current


def _path_difference(start: int, first: int, second: int,
                     edge_index: dict[tuple[int, int], int],
                     edge_count: int) -> np.ndarray:
    out = np.zeros(edge_count, dtype=np.int64)
    _add_path(out, start, [first, second], edge_index, +1)
    _add_path(out, start, [second, first], edge_index, -1)
    return out


def _rank_mod_prime(matrix: np.ndarray, prime: int = 1_000_003) -> int:
    a = np.asarray(matrix, dtype=np.int64).copy() % prime
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r, col] % prime), None)
        if pivot is None:
            continue
        a[[rank, pivot]] = a[[pivot, rank]]
        inv = pow(int(a[rank, col]), -1, prime)
        a[rank] = (a[rank] * inv) % prime
        for r in range(rows):
            if r != rank and a[r, col]:
                a[r] = (a[r] - int(a[r, col]) * a[rank]) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def _fill_closed_axis_word(start: int, word: list[int], target,
                           edge_index: dict[tuple[int, int], int],
                           face_index: dict[tuple[int, int, int, int], int]
                           ) -> tuple[np.ndarray, int]:
    word = list(word)
    end = _add_path(np.zeros(len(target.edges), dtype=np.int64), start,
                    word, edge_index)
    if end != start:
        raise ValueError("face boundary word is not closed")
    coefficients = np.zeros(len(target.faces), dtype=np.int64)
    swaps = 0
    while True:
        changed = False
        current = start
        for pos in range(len(word) - 1):
            a, b = word[pos], word[pos + 1]
            if _AXIS_POSITION[a] > _AXIS_POSITION[b]:
                difference = _path_difference(
                    current, a, b, edge_index, len(target.edges)
                )
                low, high = sorted((a, b))
                base = current & ~_axis_mask(low) & ~_axis_mask(high)
                square = (
                    base, base ^ _axis_mask(low),
                    base ^ _axis_mask(low) ^ _axis_mask(high),
                    base ^ _axis_mask(high),
                )
                square_id = face_index[square]
                boundary = target.boundary_2[:, square_id].astype(np.int64)
                nonzero = np.flatnonzero(boundary)
                pivot = int(nonzero[0])
                if difference[pivot] % boundary[pivot]:
                    raise AssertionError("swap is not an integral square boundary")
                coefficient = int(difference[pivot] // boundary[pivot])
                if not np.array_equal(difference, coefficient * boundary):
                    raise AssertionError("axis swap failed exact square-boundary check")
                coefficients[square_id] += coefficient
                word[pos], word[pos + 1] = b, a
                swaps += 1
                changed = True
                break
            current ^= _axis_mask(a)
        if not changed:
            break
    for axis in range(6):
        if word.count(axis) % 2:
            raise AssertionError("closed Q6 word has odd parity in an axis")
    if any(word[2*k] != word[2*k+1] for k in range(len(word)//2)):
        raise AssertionError("sorted closed word did not reduce to backtracks")
    return coefficients, swaps


def audit_nondegenerate_kathara_q6_chain_map() -> NondegenerateKatharaQ6Audit:
    source = audit_kathara_cellular_complex()
    target = build_q6_complex()
    labels = tuple(kathara_vertex_to_q6_label(v) for v in source.vertices)
    if len(set(labels)) != len(labels):
        raise AssertionError("Kathara vertex chart is not injective")
    qedge_index = {edge: i for i, edge in enumerate(target.edges)}
    qface_index = {face: i for i, face in enumerate(target.faces)}
    R0 = np.zeros((len(target.vertices), len(source.vertices)), dtype=np.int64)
    for j, label in enumerate(labels):
        R0[label, j] = 1
    R1 = np.zeros((len(target.edges), len(source.edges)), dtype=np.int64)
    for j, (u, v) in enumerate(source.edges):
        start, end = labels[u-1], labels[v-1]
        _add_path(R1[:, j], start, _path_axes(start, end), qedge_index)
        if start != end and np.count_nonzero(R1[:, j]) == 0:
            raise AssertionError("nontrivial source edge mapped to zero")
    R2 = np.zeros((len(target.faces), len(source.faces)), dtype=np.int64)
    boundary_l1, fill_l1 = [], []
    for face_id, face in enumerate(source.faces):
        word = []
        for u, v in zip(face, face[1:] + face[:1]):
            word.extend(_oriented_edge_path_axes(u, v, labels))
        fill, _ = _fill_closed_axis_word(
            labels[face[0]-1], word, target, qedge_index, qface_index
        )
        R2[:, face_id] = fill
        mapped_boundary = R1 @ np.asarray(source.b2, dtype=np.int64)[:, face_id]
        filled_boundary = target.boundary_2.astype(np.int64) @ fill
        if not np.array_equal(mapped_boundary, filled_boundary):
            raise AssertionError(f"face {face_id} failed exact filling identity")
        boundary_l1.append(int(np.sum(np.abs(mapped_boundary))))
        fill_l1.append(int(np.sum(np.abs(fill))))
    lhs0 = target.boundary_1.astype(np.int64) @ R1
    rhs0 = R0 @ np.asarray(source.b1, dtype=np.int64)
    lhs1 = target.boundary_2.astype(np.int64) @ R2
    rhs1 = R1 @ np.asarray(source.b2, dtype=np.int64)
    identity0, identity1 = np.array_equal(lhs0, rhs0), np.array_equal(lhs1, rhs1)
    if not identity0 or not identity1:
        raise AssertionError("non-degenerate chain-map identities failed")
    rank_mod = _rank_mod_prime(R2)
    if rank_mod != len(source.faces):
        raise AssertionError("R2 is not full column rank over the certified prime field")
    seed = load_rgs_seed_data()
    coordinates = seed.coordinates[np.asarray(labels, dtype=np.int64)].copy()
    lengths = [float(np.linalg.norm(coordinates[v-1]-coordinates[u-1]))
               for u, v in source.edges]
    if min(lengths) <= 0:
        raise AssertionError("mapped source edge has zero RGS chord length")
    return NondegenerateKatharaQ6Audit(
        labels, R0, R1, R2, identity0, identity1, rank_mod,
        tuple(boundary_l1), tuple(fill_l1), coordinates, min(lengths),
        "INTEGRAL NON-DEGENERATE CHAIN MAP; fixed rank-selected axis convention; "
        "RGS-equivariant canonicality and compatible SU(2) transport remain open."
    )
