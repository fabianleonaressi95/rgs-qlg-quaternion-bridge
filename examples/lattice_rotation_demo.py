
"""Demonstration of SU(2) lattice holonomies and gauge invariance.

Run from the repository root:

    python examples/lattice_rotation_demo.py
"""

import sys
from pathlib import Path

import numpy as np

# Allow execution directly from the repository root.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.holonomies import SU2Lattice, dagger, is_su2


def su2_rotation(axis, angle):
    """Return an SU(2) matrix for a rotation by angle around axis."""
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)

    if norm == 0.0:
        raise ValueError("Rotation axis cannot be zero.")

    x, y, z = axis / norm
    c = np.cos(angle / 2.0)
    s = np.sin(angle / 2.0)

    return np.array(
        [
            [c - 1j * s * z, -s * y - 1j * s * x],
            [s * y - 1j * s * x, c + 1j * s * z],
        ],
        dtype=complex,
    )


def main():
    np.set_printoptions(precision=8, suppress=True)

    # 1. Build a triangular lattice.
    lattice = SU2Lattice(vertices=[0, 1, 2])

    lattice.add_link(
        0, 1,
        su2_rotation([1, 0, 0], 0.37),
    )
    lattice.add_link(
        1, 2,
        su2_rotation([0, 1, 0], 0.61),
    )
    lattice.add_link(
        2, 0,
        su2_rotation([0, 0, 1], 0.83),
    )

    path = [0, 1, 2, 0]

    # 2. Calculate the plaquette holonomy.
    holonomy = lattice.holonomy(path)
    wilson = lattice.wilson_trace(path)

    print("=" * 68)
    print("SU(2) LATTICE HOLONOMY DEMONSTRATION")
    print("=" * 68)

    print("\nPlaquette path:", path)
    print("\nHolonomy:")
    print(holonomy)

    print("\nHolonomy checks:")
    print("  SU(2):", is_su2(holonomy))
    print(
        "  Unitarity residual:",
        np.linalg.norm(
            dagger(holonomy) @ holonomy - np.eye(2)
        ),
    )
    print("  Determinant:", np.linalg.det(holonomy))

    print("\nNormalized Wilson trace:", wilson)

    # 3. Choose an independent SU(2) gauge transformation at each vertex.
    gauge = {
        0: su2_rotation([1, 1, 0], 0.23),
        1: su2_rotation([0, 1, 1], 0.47),
        2: su2_rotation([1, 0, 1], 0.71),
    }

    transformed_lattice = lattice.gauge_transform(gauge)

    transformed_holonomy = transformed_lattice.holonomy(path)
    transformed_wilson = transformed_lattice.wilson_trace(path)

    # 4. Verify the based holonomy conjugation law.
    expected_holonomy = (
        gauge[0] @ holonomy @ dagger(gauge[0])
    )

    conjugation_residual = np.linalg.norm(
        transformed_holonomy - expected_holonomy
    )

    wilson_residual = abs(wilson - transformed_wilson)

    print("\nGauge covariance:")
    print("  Holonomy conjugation residual:", conjugation_residual)
    print("  Original Wilson trace:", wilson)
    print("  Transformed Wilson trace:", transformed_wilson)
    print("  Wilson invariance residual:", wilson_residual)

    # 5. Numerical assertions.
    tolerance = 1e-10

    assert is_su2(holonomy)
    assert is_su2(transformed_holonomy)
    assert conjugation_residual < tolerance
    assert wilson_residual < tolerance

    print("\n" + "=" * 68)
    print("ALL DEMONSTRATION CHECKS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    main()
```
