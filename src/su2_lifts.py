"""SU(2) representation and lifts of SO(3) rotations.

Quaternion convention:
    q = (w, x, y, z), scalar first.

Matrix representation:
    U(q) = [[w - i*z, -y - i*x],
            [y - i*x,  w + i*z]]

For unit quaternions:
    U(q) is in SU(2)
    U(q1*q2) = U(q1) @ U(q2)
    U(-q) = -U(q)

The SO(3) projection identifies q and -q.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

try:
    from .quaternions import Quaternion, DEFAULT_TOLERANCE
    from .rotations import quaternion_to_matrix
except ImportError:
    from quaternions import Quaternion, DEFAULT_TOLERANCE
    from rotations import quaternion_to_matrix


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def quaternion_to_su2(q: Quaternion) -> ComplexArray:
    """Map a unit quaternion to its 2x2 SU(2) matrix.

    The input is normalized internally.
    """
    q = q.normalized()
    w, x, y, z = q.as_tuple()

    return np.array(
        [
            [w - 1j * z, -y - 1j * x],
            [y - 1j * x, w + 1j * z],
        ],
        dtype=np.complex128,
    )


def su2_diagnostics(
    matrix: NDArray[np.complex128],
) -> dict[str, float]:
    """Return unitarity and determinant diagnostics."""
    matrix = np.asarray(matrix, dtype=np.complex128)

    if matrix.shape != (2, 2):
        raise ValueError("An SU(2) matrix must have shape (2, 2).")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix entries must be finite.")

    identity = np.eye(2, dtype=np.complex128)
    determinant = np.linalg.det(matrix)

    return {
        "unitarity_error": float(
            np.linalg.norm(matrix.conj().T @ matrix - identity, ord="fro")
        ),
        "determinant_real": float(determinant.real),
        "determinant_imag": float(determinant.imag),
        "determinant_error": float(abs(determinant - 1.0)),
    }


def is_su2(
    matrix: NDArray[np.complex128],
    tolerance: float = 1e-10,
) -> bool:
    """Check unitarity and determinant one."""
    if tolerance < 0:
        raise ValueError("Tolerance must be non-negative.")

    try:
        diagnostics = su2_diagnostics(matrix)
    except (TypeError, ValueError):
        return False

    return (
        diagnostics["unitarity_error"] <= tolerance
        and diagnostics["determinant_error"] <= tolerance
    )


def su2_inverse(
    matrix: NDArray[np.complex128],
) -> ComplexArray:
    """Return the group inverse U^{-1} = U dagger."""
    matrix = np.asarray(matrix, dtype=np.complex128)

    if matrix.shape != (2, 2):
        raise ValueError("Expected a 2x2 matrix.")

    return matrix.conj().T


def su2_to_quaternion(
    matrix: NDArray[np.complex128],
    tolerance: float = 1e-10,
) -> Quaternion:
    """Recover the scalar-first unit quaternion from an SU(2) matrix.

    The matrix must follow this module's representation convention.
    """
    matrix = np.asarray(matrix, dtype=np.complex128)

    if not is_su2(matrix, tolerance):
        raise ValueError("Input matrix is not SU(2) within tolerance.")

    # U = [[w - iz, -y - ix],
    #      [y - ix,  w + iz]]
    w = 0.5 * (matrix[0, 0].real + matrix[1, 1].real)
    z = 0.5 * (matrix[1, 1].imag - matrix[0, 0].imag)

    x = -0.5 * (matrix[0, 1].imag + matrix[1, 0].imag)
    y = 0.5 * (matrix[1, 0].real - matrix[0, 1].real)

    q = Quaternion(w, x, y, z)

    # Remove numerical round-off while preserving the lift sign.
    return q.normalized()


def opposite_lift(q: Quaternion) -> Quaternion:
    """Return the other SU(2) lift: -q.

    Both q and -q project to the same SO(3) rotation.
    """
    return -q


def su2_lifts_of_rotation(
    q: Quaternion,
) -> tuple[ComplexArray, ComplexArray]:
    """Return both SU(2) lifts of the SO(3) rotation represented by q."""
    q_unit = q.normalized()
    return (
        quaternion_to_su2(q_unit),
        quaternion_to_su2(-q_unit),
    )


def su2_to_so3(
    matrix: NDArray[np.complex128],
    tolerance: float = 1e-10,
) -> FloatArray:
    """Project an SU(2) matrix to its associated SO(3) rotation.

    The projection is performed through the corresponding quaternion.
    """
    q = su2_to_quaternion(matrix, tolerance)
    return quaternion_to_matrix(q)


def su2_multiplication_residual(
    q1: Quaternion,
    q2: Quaternion,
) -> float:
    """Test the representation homomorphism numerically."""
    left = quaternion_to_su2(q1 * q2)
    right = quaternion_to_su2(q1) @ quaternion_to_su2(q2)

    return float(np.linalg.norm(left - right, ord="fro"))


def lift_sign_residual(
    q: Quaternion,
) -> float:
    """Verify that opposite lifts are negatives as SU(2) matrices."""
    u_plus = quaternion_to_su2(q)
    u_minus = quaternion_to_su2(-q)

    return float(np.linalg.norm(u_minus + u_plus, ord="fro"))


def projection_sign_residual(
    q: Quaternion,
) -> float:
    """Verify that q and -q define the same SO(3) matrix."""
    r_plus = quaternion_to_matrix(q)
    r_minus = quaternion_to_matrix(-q)

    return float(np.linalg.norm(r_plus - r_minus, ord="fro"))


def su2_so3_roundtrip_residual(
    q: Quaternion,
) -> float:
    """Measure the SU(2) -> SO(3) projection against R(q)."""
    expected = quaternion_to_matrix(q)
    recovered = su2_to_so3(quaternion_to_su2(q))

    return float(np.linalg.norm(expected - recovered, ord="fro"))


def quaternion_inverse_representation_residual(
    q: Quaternion,
) -> float:
    """Verify U(q^{-1}) = U(q) dagger."""
    q_inverse = q.inverse()

    left = quaternion_to_su2(q_inverse)
    right = quaternion_to_su2(q).conj().T

    return float(np.linalg.norm(left - right, ord="fro"))


def validate_su2_lift(
    q: Quaternion,
    tolerance: float = 1e-10,
) -> dict[str, float | bool]:
    """Run a compact diagnostic suite for a quaternionic SU(2) lift."""
    u = quaternion_to_su2(q)
    diagnostics = su2_diagnostics(u)

    multiplication_test = su2_multiplication_residual(
        q,
        Quaternion(0.9238795325112867, 0.0, 0.3826834323650898, 0.0),
    )

    roundtrip = su2_so3_roundtrip_residual(q)
    inverse = quaternion_inverse_representation_residual(q)
    sign = projection_sign_residual(q)

    return {
        **diagnostics,
        "is_su2": is_su2(u, tolerance),
        "multiplication_residual": multiplication_test,
        "so3_roundtrip_residual": roundtrip,
        "inverse_representation_residual": inverse,
        "projection_sign_residual": sign,
        "all_checks_pass": (
            is_su2(u, tolerance)
            and multiplication_test <= tolerance
            and roundtrip <= tolerance
            and inverse <= tolerance
            and sign <= tolerance
        ),
    }
