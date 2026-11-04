"""Canonical incidence data for the six-dimensional hypercube Q6.

Vertex labels use six-bit big-endian lexicographic indexing:
    vertex = int(bitstring, 2)
Axis 0 is the most significant bit and axis 5 the least significant bit.

Edges are oriented from the smaller integer label to the larger label.
Each square face is oriented by the ordered axis pair (a, b), a < b,
with boundary v -> v^a -> v^(a,b) -> v^b -> v.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Q6Complex:
    vertices: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    faces: tuple[tuple[int, int, int, int], ...]
    boundary_1: np.ndarray
    boundary_2: np.ndarray


def _axis_mask(axis: int) -> int:
    if not 0 <= axis < 6:
        raise ValueError("Q6 axis must be in range 0..5.")
    return 1 << (5 - axis)


def build_q6_complex() -> Q6Complex:
    """Build Q6 vertices, canonically oriented edges and square faces."""
    vertices = tuple(range(64))

    edge_set: set[tuple[int, int]] = set()
    for v in vertices:
        for axis in range(6):
            w = v ^ _axis_mask(axis)
            if v < w:
                edge_set.add((v, w))
    edges = tuple(sorted(edge_set))
    edge_index = {edge: k for k, edge in enumerate(edges)}

    b1 = np.zeros((64, len(edges)), dtype=np.int8)
    for k, (u, v) in enumerate(edges):
        b1[u, k] = -1
        b1[v, k] = 1

    faces_list: list[tuple[int, int, int, int]] = []
    face_boundaries: list[tuple[tuple[int, int], ...]] = []

    for a in range(6):
        ma = _axis_mask(a)
        for b in range(a + 1, 6):
            mb = _axis_mask(b)
            for v in vertices:
                # The base vertex has zero bits along both face axes.
                if (v & ma) or (v & mb):
                    continue
                va = v ^ ma
                vab = va ^ mb
                vb = v ^ mb
                faces_list.append((v, va, vab, vb))
                face_boundaries.append(((v, va), (va, vab), (vab, vb), (vb, v)))

    faces = tuple(faces_list)
    b2 = np.zeros((len(edges), len(faces)), dtype=np.int8)

    for face_id, boundary in enumerate(face_boundaries):
        for source, target in boundary:
            canonical = (min(source, target), max(source, target))
            sign = 1 if (source, target) == canonical else -1
            b2[edge_index[canonical], face_id] += sign

    return Q6Complex(
        vertices=vertices,
        edges=edges,
        faces=faces,
        boundary_1=b1,
        boundary_2=b2,
    )


def audit_q6_complex() -> dict[str, int | float | bool]:
    """Return structural counts and the exact chain-complex residual."""
    complex_ = build_q6_complex()
    residual = complex_.boundary_1.astype(np.int64) @ complex_.boundary_2.astype(np.int64)
    return {
        "vertices": len(complex_.vertices),
        "edges": len(complex_.edges),
        "faces": len(complex_.faces),
        "degree_min": int(np.min(np.sum(np.abs(complex_.boundary_1), axis=1))),
        "degree_max": int(np.max(np.sum(np.abs(complex_.boundary_1), axis=1))),
        "boundary_squared_max_abs": int(np.max(np.abs(residual))),
        "boundary_squared_exact_zero": bool(np.count_nonzero(residual) == 0),
    }

def export_q6_incidence_bundle(path: str) -> dict[str, object]:
    """Export the source-defined canonical Q6 incidence ordering as NPZ + hash.

    This schema is explicit and reproducible, but is not asserted byte-compatible
    with the separately recovered historical NPZ unless an array-by-array
    ordering comparison is performed.
    """
    import hashlib
    import json
    from pathlib import Path

    complex_ = build_q6_complex()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "schema": "RGS_Q6_CANONICAL_INCIDENCE_V1",
        "vertex_order": "integer labels 0..63; six-bit big-endian",
        "edge_order": "lexicographically sorted (u,v), u<v",
        "face_order": "axis pairs a<b, then base vertex ascending; cyclic boundary v,va,vab,vb",
        "boundary_convention": "B1:C1->C0, B2:C2->C1; B1@B2=0",
        "historical_npz_byte_compatibility": "not asserted",
    }
    np.savez_compressed(
        output,
        vertices=np.asarray(complex_.vertices, dtype=np.int16),
        edges=np.asarray(complex_.edges, dtype=np.int16),
        faces=np.asarray(complex_.faces, dtype=np.int16),
        B1=np.asarray(complex_.boundary_1, dtype=np.int8),
        B2=np.asarray(complex_.boundary_2, dtype=np.int8),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    return {
        **metadata,
        "path": str(output),
        "sha256": digest,
        "shapes": {
            "vertices": [64],
            "edges": [192, 2],
            "faces": [240, 4],
            "B1": [64, 192],
            "B2": [192, 240],
        },
    }
