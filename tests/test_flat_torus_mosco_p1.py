import numpy as np

from src.flat_torus_refinement import build_flat_torus_mesh
from src.flat_torus_embedding import J_h


def test_p1_mosco_recovery_sequence_for_smooth_periodic_mode():
    side = 1.0
    levels = (4, 8, 16)
    energy_errors = []
    l2_errors = []
    k = 2.0 * np.pi / side
    exact_energy = 0.5 * (k * k + 0.25 * k * k) * side**2

    for n in levels:
        mesh = build_flat_torus_mesh(n, dimension=2, side_length=side)
        x = mesh.vertices
        values = np.sin(k * x[:, 0]) + 0.5 * np.cos(k * x[:, 1])
        energy = float(values @ mesh.stiffness0 @ values)
        energy_errors.append(abs(energy - exact_energy) / exact_energy)

        count = 4 * n
        axis = (np.arange(count) + 0.5) * side / count
        xx, yy = np.meshgrid(axis, axis, indexing="ij")
        points = np.column_stack((xx.ravel(), yy.ravel()))
        reconstructed = J_h(mesh, values).evaluate(points)
        exact = np.sin(k * points[:, 0]) + 0.5 * np.cos(k * points[:, 1])
        l2_errors.append(float(np.sqrt(np.mean((reconstructed - exact) ** 2))))

    assert energy_errors[-1] < energy_errors[0]
    assert l2_errors[-1] < l2_errors[0]
    assert all(np.isfinite(energy_errors))
    assert all(np.isfinite(l2_errors))
