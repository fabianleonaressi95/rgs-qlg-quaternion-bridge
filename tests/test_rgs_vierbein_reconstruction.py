import numpy as np

from examples.reconstruct_rgs_vierbein_tetrads import (
    ETA,
    ellis_diagonal_tetrad,
    kerr_schild_tetrad,
    metric_from_coframe,
    reconstruct_spatial_vierbein,
    static_3plus1_tetrad,
)
from src.rgs_ellis_refinement import load_rgs_seed_data


def test_reconstructs_archived_64_site_seed_and_positive_spatial_frames():
    seed = load_rgs_seed_data()
    assert seed.coordinates.shape == (64, 3)
    e, c = reconstruct_spatial_vierbein(seed.coordinates)
    assert e.shape == (64, 3, 3)
    assert np.all(np.linalg.det(e) > 0)
    assert np.allclose(e @ np.swapaxes(e, 1, 2), c, atol=2e-12)


def test_static_lorentzian_lift_reproduces_site_metrics():
    seed = load_rgs_seed_data()
    e, _ = reconstruct_spatial_vierbein(seed.coordinates)
    tetrad = static_3plus1_tetrad(e)
    metric = metric_from_coframe(tetrad)
    expected = np.zeros_like(metric)
    expected[:, 0, 0] = -1.0
    expected[:, 1:, 1:] = np.einsum("...am,...an->...mn", e, e)
    assert tetrad.shape == (64, 4, 4)
    assert np.allclose(metric, expected, atol=2e-12)
    assert np.all(np.linalg.det(tetrad) > 0)


def test_ellis_tetrad_reproduces_wormhole_metric():
    u = np.linspace(-2.0, 2.0, 17)
    theta = np.linspace(0.2, np.pi - 0.2, 17)
    r0 = 0.8292
    e = ellis_diagonal_tetrad(u, theta, r0=r0)
    g = metric_from_coframe(e)
    expected = np.zeros_like(g)
    radius = np.sqrt(u*u + r0*r0)
    expected[:, 0, 0] = -1.0
    expected[:, 1, 1] = 1.0
    expected[:, 2, 2] = radius*radius
    expected[:, 3, 3] = radius*radius*np.sin(theta)**2
    assert np.allclose(g, expected, atol=2e-12)


def test_kerr_schild_tetrad_reproduces_exact_kerr_schild_metric():
    rng = np.random.default_rng(71)
    xyz = rng.normal(size=(32, 3))
    xyz *= np.linspace(2.0, 5.0, 32)[:, None] / np.linalg.norm(xyz, axis=1)[:, None]
    e, g_exact, aux = kerr_schild_tetrad(xyz, mass=1.0, spin=0.6)
    assert e.shape == (32, 4, 4)
    assert np.allclose(metric_from_coframe(e), g_exact, atol=2e-12)
    assert np.allclose(np.linalg.det(e), 1.0, atol=2e-12)
    assert np.all(aux[:, 0] > 0.0)
    assert np.allclose(ETA, np.diag([-1.0, 1.0, 1.0, 1.0]))
