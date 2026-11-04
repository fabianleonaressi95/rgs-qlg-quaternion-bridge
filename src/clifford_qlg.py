"""Lorentzian Clifford algebra and 2I-equivariant spinor kernels.

Convention: eta = diag(-1,+1,+1,+1), gamma^mu gamma^nu +
gamma^nu gamma^mu = 2 eta^{mu nu} I_4. The matrices are the complex
4x4 Dirac representation obtained by multiplying the standard (+---)
Dirac matrices by i.

The finite 2I action is the restriction of the spatial Spin(3) action:
S(q) = diag(U(q), U(q)), where U is the repository's scalar-first SU(2)
lift. This is a spatial rotation subgroup action, not a representation
of Lorentz boosts and not a lift of Möbius reflections.
"""
from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np

from src.binary_icosahedral import binary_icosahedral_elements, quaternion_to_so3
from src.representations_2i import quaternion_to_fundamental_su2


ETA = np.diag([-1.0, 1.0, 1.0, 1.0]).astype(float)
I2 = np.eye(2, dtype=np.complex128)
I4 = np.eye(4, dtype=np.complex128)

PAULI = (
    np.array([[0, 1], [1, 0]], dtype=np.complex128),
    np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    np.array([[1, 0], [0, -1]], dtype=np.complex128),
)


def gamma_matrices() -> tuple[np.ndarray, ...]:
    """Return gamma^0,...,gamma^3 for signature (-,+,+,+)."""
    gamma0 = 1j * np.block([[I2, np.zeros((2, 2))],
                            [np.zeros((2, 2)), -I2]])
    spatial = tuple(
        1j * np.block([[np.zeros((2, 2)), sigma],
                       [-sigma, np.zeros((2, 2))]])
        for sigma in PAULI
    )
    return (gamma0, *spatial)


def clifford_blade_basis(
    gammas: Iterable[np.ndarray] | None = None,
) -> tuple[tuple[tuple[int, ...], ...], np.ndarray]:
    """Return the 16 ordered gamma blades and their flattened column basis."""
    gs = tuple(gamma_matrices() if gammas is None else gammas)
    if len(gs) != 4 or any(np.asarray(g).shape != (4, 4) for g in gs):
        raise ValueError("Exactly four 4x4 gamma matrices are required.")
    labels = []
    blades = []
    for grade in range(5):
        for indices in combinations(range(4), grade):
            blade = I4.copy()
            for index in indices:
                blade = blade @ gs[index]
            labels.append(indices)
            blades.append(blade)
    matrix = np.column_stack([blade.reshape(-1) for blade in blades])
    if np.linalg.matrix_rank(matrix, tol=1e-11) != 16:
        raise ArithmeticError("The 16 Clifford blades are not linearly independent.")
    return tuple(labels), matrix


def clifford_grade_projections(
    value: np.ndarray,
    gammas: Iterable[np.ndarray] | None = None,
) -> dict[int, np.ndarray]:
    """Project a 4x4 complex matrix onto grades 0 through 4."""
    value = np.asarray(value, dtype=np.complex128)
    if value.shape != (4, 4):
        raise ValueError("A Clifford-algebra element must be a 4x4 matrix.")
    gs = tuple(gamma_matrices() if gammas is None else gammas)
    labels, basis = clifford_blade_basis(gs)
    coefficients, *_ = np.linalg.lstsq(basis, value.reshape(-1), rcond=None)
    projections = {grade: np.zeros((4, 4), dtype=np.complex128)
                   for grade in range(5)}
    for coefficient, indices in zip(coefficients, labels):
        blade = I4.copy()
        for index in indices:
            blade = blade @ gs[index]
        projections[len(indices)] += coefficient * blade
    return projections


def clifford_relation_audit(tol: float = 2e-12) -> dict[str, object]:
    """Check all 16 Clifford anticommutator relations."""
    gs = gamma_matrices()
    residuals = {}
    for mu in range(4):
        for nu in range(4):
            lhs = gs[mu] @ gs[nu] + gs[nu] @ gs[mu]
            rhs = 2.0 * ETA[mu, nu] * I4
            residuals[f"{mu},{nu}"] = float(np.linalg.norm(lhs-rhs, ord="fro"))
    labels, basis = clifford_blade_basis(gs)
    return {
        "audit": "CLIFFORD_CL_1_3",
        "signature": [-1, 1, 1, 1],
        "gamma_shape": [4, 4],
        "blade_count": len(labels),
        "blade_basis_rank": int(np.linalg.matrix_rank(basis, tol=1e-11)),
        "max_anticommutator_residual": max(residuals.values()),
        "residuals": residuals,
        "passed": max(residuals.values()) <= tol and len(labels) == 16,
    }


def spinor_lift_2i(q: np.ndarray) -> np.ndarray:
    """Embed the defining SU(2) representation as spatial Spin(3) on Dirac spinors."""
    u = quaternion_to_fundamental_su2(np.asarray(q, dtype=float))
    return np.block([[u, np.zeros((2, 2), dtype=complex)],
                     [np.zeros((2, 2), dtype=complex), u]])


def spinor_group_action_audit(tol: float = 3e-10) -> dict[str, object]:
    """Check homomorphism, unitarity, central action and gamma-vector covariance."""
    elements = binary_icosahedral_elements()
    gammas = gamma_matrices()
    hom = 0.0
    unit = 0.0
    covariance = 0.0
    for q in elements:
        s = spinor_lift_2i(q)
        unit = max(unit, float(np.linalg.norm(s.conj().T @ s-I4, ord="fro")))
        inverse = s.conj().T
        rotation = quaternion_to_so3(q)
        for i in range(3):
            transformed = s @ gammas[i+1] @ inverse
            expected = sum(rotation[j, i] * gammas[j+1] for j in range(3))
            covariance = max(covariance, float(np.linalg.norm(transformed-expected, ord="fro")))
    # Exhaustive 2I closure is audited in binary_icosahedral; sample deterministic
    # products here to check that the 4x4 spin lift preserves multiplication.
    from src.binary_icosahedral import qmul
    rng = np.random.default_rng(20261103)
    for i, j in rng.integers(0, len(elements), size=(512, 2)):
        lhs = spinor_lift_2i(qmul(elements[i], elements[j]))
        rhs = spinor_lift_2i(elements[i]) @ spinor_lift_2i(elements[j])
        hom = max(hom, float(np.linalg.norm(lhs-rhs, ord="fro")))
    central = spinor_lift_2i(-np.array([1.0, 0.0, 0.0, 0.0]))
    return {
        "audit": "2I_SPINOR_ACTION_ON_CL_1_3",
        "group_order": len(elements),
        "spinor_dimension": 4,
        "homomorphism_max_residual": hom,
        "unitarity_max_residual": unit,
        "gamma_spatial_covariance_max_residual": covariance,
        "central_minus_one_residual_to_minus_I4": float(np.linalg.norm(central+I4, ord="fro")),
        "passed": max(hom, unit, covariance, float(np.linalg.norm(central+I4))) <= tol,
        "scope": "spatial Spin(3) subgroup only; no boosts or reflections",
    }


def intertwiner_basis(
    output_matrices: Iterable[np.ndarray],
    input_matrices: Iterable[np.ndarray],
    rtol: float = 1e-10,
) -> tuple[np.ndarray, int, float]:
    """Compute K satisfying R_out(g) K = K R_in(g) for paired group elements."""
    outputs = tuple(np.asarray(x, dtype=complex) for x in output_matrices)
    inputs = tuple(np.asarray(x, dtype=complex) for x in input_matrices)
    if len(outputs) == 0 or len(outputs) != len(inputs):
        raise ValueError("Output and input representations need equal nonzero lengths.")
    no, ni = outputs[0].shape[0], inputs[0].shape[0]
    if any(x.shape != (no, no) for x in outputs) or any(
        x.shape != (ni, ni) for x in inputs
    ):
        raise ValueError("Each representation must have a fixed square matrix shape.")
    blocks = []
    for rout, rin in zip(outputs, inputs):
        block = np.empty((no*ni, no*ni), dtype=complex)
        for k in range(no*ni):
            e = np.zeros((no, ni), dtype=complex)
            e.flat[k] = 1.0
            block[:, k] = (rout @ e - e @ rin).reshape(-1)
        blocks.append(block)
    constraint = np.vstack(blocks)
    _, singular, vh = np.linalg.svd(constraint, full_matrices=True)
    threshold = rtol * (singular[0] if singular.size else 1.0)
    rank = int(np.count_nonzero(singular > threshold))
    basis = vh.conj().T[:, rank:].T.reshape((-1, no, ni))
    residual = 0.0
    for k in basis:
        for rout, rin in zip(outputs, inputs):
            residual = max(residual, float(np.linalg.norm(rout@k-k@rin, ord="fro")))
    return basis, int(len(basis)), residual


def spinor_steerable_kernel_audit(rtol: float = 1e-10) -> dict[str, object]:
    """Find all 2I-equivariant linear maps from Dirac spinors to Dirac spinors."""
    elements = binary_icosahedral_elements()
    reps = tuple(spinor_lift_2i(q) for q in elements)
    basis, dimension, residual = intertwiner_basis(reps, reps, rtol=rtol)
    return {
        "audit": "2I_SPINOR_STEERABLE_KERNELS",
        "group_order": len(elements),
        "input_dimension": 4,
        "output_dimension": 4,
        "intertwiner_dimension": dimension,
        "basis_shape": list(basis.shape),
        "max_intertwining_residual": residual,
        "expected_dimension": 4,
        "interpretation": "Dirac spinor restricted to Spin(3) is two equivalent SU(2) doublets; its commutant is Mat_2(C).",
        "passed": dimension == 4 and residual <= 2e-9,
    }

def clifford_grade_action_audit(tol: float = 3e-9) -> dict[str, object]:
    """Check that the 2I spatial spin action preserves Clifford grades 0..4.

    This is the adjoint action S B S^{-1} on the matrix Clifford algebra.
    It tests the spatial Spin(3) subgroup only; it does not implement Pin
    reflections or the orientation-reversing Möbius frame transitions.
    """
    labels, _ = clifford_blade_basis()
    elements = binary_icosahedral_elements()
    maximum_leakage = 0.0
    maximum_reconstruction = 0.0
    per_grade = {grade: 0.0 for grade in range(5)}
    for q in elements:
        spin = spinor_lift_2i(q)
        inverse = spin.conj().T
        for indices in labels:
            grade = len(indices)
            blade = I4.copy()
            for index in indices:
                blade = blade @ gamma_matrices()[index]
            transformed = spin @ blade @ inverse
            parts = clifford_grade_projections(transformed)
            leakage = max(
                (float(np.linalg.norm(value, ord="fro"))
                 for other_grade, value in parts.items() if other_grade != grade),
                default=0.0,
            )
            reconstruction = float(np.linalg.norm(
                sum(parts.values()) - transformed, ord="fro"
            ))
            per_grade[grade] = max(per_grade[grade], leakage)
            maximum_leakage = max(maximum_leakage, leakage)
            maximum_reconstruction = max(maximum_reconstruction, reconstruction)
    return {
        "audit": "2I_CLIFFORD_GRADE_ACTION",
        "group_order": len(elements),
        "grades_checked": [0, 1, 2, 3, 4],
        "max_grade_leakage": maximum_leakage,
        "max_grade_reconstruction_residual": maximum_reconstruction,
        "per_grade_max_leakage": {str(k): v for k, v in per_grade.items()},
        "passed": maximum_leakage <= tol and maximum_reconstruction <= tol,
        "scope": "spatial Spin(3) adjoint action; no reflection/Pin lift",
    }
