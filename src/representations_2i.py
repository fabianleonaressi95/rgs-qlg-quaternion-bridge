"""Explicit 2I representations and numerical intertwiner spaces.

The three representations implemented here are:
  * trivial: the 1-dimensional trivial representation;
  * fundamental: the defining 2-dimensional complex SU(2) representation;
  * vector: the 3-dimensional real SO(3) representation.

These are representations of 2I, not a complete catalogue of all irreducibles.
The vector representation factors through 2I/{+/-1} ~= A5. The fundamental
representation does not: rho(-1)=-I_2. No reflection is represented here.
"""
from __future__ import annotations

from typing import Literal

import numpy as np

from src.binary_icosahedral import (
    IDENTITY,
    binary_icosahedral_elements,
    qmul,
    quaternion_to_so3,
)

RepresentationName = Literal["trivial", "fundamental", "vector"]


def quaternion_to_fundamental_su2(q: np.ndarray) -> np.ndarray:
    """Defining 2x2 complex representation, consistent with src.su2_lifts."""
    w, x, y, z = np.asarray(q, dtype=float)
    return np.array(
        [[w - 1j * z, -y - 1j * x],
         [y - 1j * x, w + 1j * z]],
        dtype=np.complex128,
    )


def representation_matrix(
    q: np.ndarray,
    name: RepresentationName,
) -> np.ndarray:
    """Evaluate one of the explicitly supported 2I representations."""
    q = np.asarray(q, dtype=float)
    if q.shape != (4,):
        raise ValueError("A group element must be a scalar-first 4-vector.")
    if name == "trivial":
        return np.ones((1, 1), dtype=np.complex128)
    if name == "fundamental":
        return quaternion_to_fundamental_su2(q)
    if name == "vector":
        return quaternion_to_so3(q).astype(np.complex128)
    raise ValueError(f"Unsupported representation {name!r}.")


def _dimension(name: RepresentationName) -> int:
    return {"trivial": 1, "fundamental": 2, "vector": 3}[name]


def representation_audit(tol: float = 2e-10) -> dict[str, object]:
    """Audit homomorphism, unitarity/orthogonality and central element action."""
    elements = binary_icosahedral_elements()
    rng = np.random.default_rng(20261103)
    indices = rng.integers(0, len(elements), size=(512, 2))
    homomorphism = {name: 0.0 for name in ("trivial", "fundamental", "vector")}
    unitarity = {name: 0.0 for name in ("trivial", "fundamental", "vector")}
    for i, j in indices:
        q1, q2 = elements[i], elements[j]
        product = qmul(q1, q2)
        for name in homomorphism:
            lhs = representation_matrix(product, name)
            rhs = representation_matrix(q1, name) @ representation_matrix(q2, name)
            homomorphism[name] = max(
                homomorphism[name], float(np.linalg.norm(lhs-rhs, ord="fro"))
            )
    for q in elements:
        for name in unitarity:
            matrix = representation_matrix(q, name)
            identity = np.eye(matrix.shape[0], dtype=complex)
            unitarity[name] = max(
                unitarity[name],
                float(np.linalg.norm(matrix.conj().T @ matrix - identity, ord="fro")),
            )

    minus_identity = -IDENTITY
    central_action = {
        name: representation_matrix(minus_identity, name).tolist()
        for name in ("trivial", "fundamental", "vector")
    }
    passed = (
        max(homomorphism.values()) <= tol
        and max(unitarity.values()) <= tol
        and np.allclose(representation_matrix(IDENTITY, "fundamental"), np.eye(2))
    )
    return {
        "audit": "2I_EXPLICIT_REPRESENTATIONS",
        "group_order": len(elements),
        "representations": {
            "trivial": {"dimension": 1, "factors_through_A5": True},
            "fundamental": {"dimension": 2, "factors_through_A5": False},
            "vector": {"dimension": 3, "factors_through_A5": True},
        },
        "homomorphism_max_residual": homomorphism,
        "unitarity_max_residual": unitarity,
        "central_minus_one_action": central_action,
        "passed": bool(passed),
    }


def intertwiner_basis(
    output_rep: RepresentationName,
    input_rep: RepresentationName,
    elements: np.ndarray | None = None,
    rtol: float = 1e-10,
) -> tuple[np.ndarray, int, float]:
    """Return an orthonormal basis for K satisfying Rout(g)K=K Rin(g).

    Basis shape is (nullity, dim_out, dim_in), over the complex numbers.
    The reported residual is the maximum Frobenius residual over all supplied
    group elements for every returned basis vector.
    """
    group = binary_icosahedral_elements() if elements is None else np.asarray(elements)
    no, ni = _dimension(output_rep), _dimension(input_rep)
    columns = []
    for q in group:
        rout = representation_matrix(q, output_rep)
        rin = representation_matrix(q, input_rep)
        block = np.empty((no * ni, no * ni), dtype=np.complex128)
        for k in range(no * ni):
            e = np.zeros((no, ni), dtype=np.complex128)
            e.flat[k] = 1.0
            block[:, k] = (rout @ e - e @ rin).ravel()
        columns.append(block)
    constraint = np.vstack(columns)
    _, singular, vh = np.linalg.svd(constraint, full_matrices=True)
    threshold = rtol * (singular[0] if singular.size else 1.0)
    rank = int(np.count_nonzero(singular > threshold))
    null_vectors = vh.conj().T[:, rank:]
    basis = null_vectors.T.reshape((-1, no, ni))

    residual = 0.0
    for matrix in basis:
        for q in group:
            residual = max(
                residual,
                float(np.linalg.norm(
                    representation_matrix(q, output_rep) @ matrix
                    - matrix @ representation_matrix(q, input_rep),
                    ord="fro",
                )),
            )
    return basis, int(basis.shape[0]), residual


def intertwiner_dimension_audit(rtol: float = 1e-10) -> dict[str, object]:
    """Compute selected Hom-space dimensions and Schur-lemma checks."""
    names = ("trivial", "fundamental", "vector")
    dimensions = {}
    residuals = {}
    for out_name in names:
        for in_name in names:
            _, dim, residual = intertwiner_basis(out_name, in_name, rtol=rtol)
            dimensions[f"{in_name}->{out_name}"] = dim
            residuals[f"{in_name}->{out_name}"] = residual
    expected = {
        "trivial->trivial": 1,
        "fundamental->fundamental": 1,
        "vector->vector": 1,
        "trivial->fundamental": 0,
        "fundamental->trivial": 0,
        "trivial->vector": 0,
        "vector->trivial": 0,
        "fundamental->vector": 0,
        "vector->fundamental": 0,
    }
    return {
        "audit": "2I_INTERTWINER_DIMENSIONS",
        "dimensions": dimensions,
        "expected_by_irreducibility": expected,
        "dimension_match": dimensions == expected,
        "max_intertwiner_residual": max(residuals.values(), default=0.0),
        "residuals": residuals,
        "note": "Selected irreducible representations only; not the full 2I character table.",
    }
