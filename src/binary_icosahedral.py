"""Explicit binary icosahedral group 2I as 120 unit quaternions.

Convention: scalar-first q=(w,x,y,z); Hamilton product; active rotations
on column vectors. The vertex set is the standard 600-cell vertex set:
8 coordinate units, 16 half-coordinate vertices, and 96 even-permutation
vertices built from (0, 1/2, phi/2, 1/(2 phi)).
"""
from __future__ import annotations

from itertools import permutations, product
import numpy as np

PHI = (1.0 + np.sqrt(5.0)) / 2.0
IDENTITY = np.array([1.0, 0.0, 0.0, 0.0])


def _parity(p: tuple[int, ...]) -> int:
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2


def binary_icosahedral_elements() -> np.ndarray:
    """Return a deterministic (120,4) array containing every element of 2I."""
    points: set[tuple[float, float, float, float]] = set()

    # 8 coordinate units: (±1,0,0,0) and coordinate permutations.
    for axis in range(4):
        for sign in (-1.0, 1.0):
            q = np.zeros(4)
            q[axis] = sign
            points.add(tuple(q))

    # 16 vertices with all coordinates ±1/2.
    for signs in product((-1.0, 1.0), repeat=4):
        points.add(tuple(0.5 * np.asarray(signs)))

    # 96 vertices: even permutations, independent signs on nonzero entries.
    base = (0.0, 0.5, PHI / 2.0, 1.0 / (2.0 * PHI))
    for p in permutations(range(4)):
        if _parity(p):
            continue
        permuted = tuple(base[i] for i in p)
        nonzero = [i for i, value in enumerate(permuted) if value != 0.0]
        for signs in product((-1.0, 1.0), repeat=3):
            q = list(permuted)
            for i, sign in zip(nonzero, signs):
                q[i] *= sign
            points.add(tuple(float(x) for x in q))

    arr = np.asarray(sorted(points), dtype=float)
    if arr.shape != (120, 4):
        raise RuntimeError(f"Expected 120 distinct vertices, got {arr.shape}")
    return arr


def qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product, scalar-first."""
    aw, ax, ay, az = np.asarray(a, dtype=float)
    bw, bx, by, bz = np.asarray(b, dtype=float)
    return np.array([
        aw*bw - ax*bx - ay*by - az*bz,
        aw*bx + ax*bw + ay*bz - az*by,
        aw*by - ax*bz + ay*bw + az*bx,
        aw*bz + ax*by - ay*bx + az*bw,
    ])


def qconj(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    return q * np.array([1.0, -1.0, -1.0, -1.0])


def quaternion_to_so3(q: np.ndarray) -> np.ndarray:
    """Active SO(3) rotation matrix induced by a unit scalar-first quaternion."""
    w, x, y, z = np.asarray(q, dtype=float)
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
        [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)],
    ])


def audit_binary_icosahedral(tol: float = 2e-10) -> dict:
    """Exhaustively audit cardinality, S3 norm, closure, inverses and SO(3) quotient."""
    elements = binary_icosahedral_elements()
    keys = {tuple(np.round(q, 10)) for q in elements}
    lookup = {key: i for i, key in enumerate(elements)}
    identity_key = tuple(np.round(IDENTITY, 10))
    if identity_key not in lookup:
        raise AssertionError("Identity missing")

    norm_error = float(np.max(np.abs(np.linalg.norm(elements, axis=1) - 1.0)))
    closure_error = 0.0
    missing_products = 0
    for a in elements:
        for b in elements:
            c = qmul(a, b)
            key = tuple(np.round(c, 10))
            if key not in keys:
                missing_products += 1
            else:
                closure_error = max(closure_error, float(np.linalg.norm(c - elements[lookup[key]])))

    inverse_missing = 0
    for q in elements:
        if tuple(np.round(qconj(q), 10)) not in keys:
            inverse_missing += 1

    # q and -q induce the same SO(3) rotation; quotient should have 60 classes.
    rotation_keys = {
        tuple(np.round(quaternion_to_so3(q), 8).ravel())
        for q in elements
    }
    orth_error = 0.0
    det_error = 0.0
    for q in elements:
        r = quaternion_to_so3(q)
        orth_error = max(orth_error, float(np.linalg.norm(r.T @ r - np.eye(3), ord="fro")))
        det_error = max(det_error, abs(float(np.linalg.det(r)) - 1.0))

    result = {
        "order_2I": int(len(elements)),
        "unit_norm_max_error": norm_error,
        "closure_missing_products": int(missing_products),
        "closure_max_matching_error": closure_error,
        "inverse_missing": int(inverse_missing),
        "distinct_SO3_rotations": int(len(rotation_keys)),
        "SO3_orthogonality_max_fro": orth_error,
        "SO3_determinant_max_error": det_error,
        "quaternion_convention": "scalar-first (w,x,y,z)",
        "rotation_convention": "active, column vectors",
    }
    if (len(elements) != 120 or norm_error > tol or missing_products or
            inverse_missing or len(rotation_keys) != 60 or
            orth_error > tol or det_error > tol):
        raise AssertionError(f"2I audit failed: {result}")
    return result
