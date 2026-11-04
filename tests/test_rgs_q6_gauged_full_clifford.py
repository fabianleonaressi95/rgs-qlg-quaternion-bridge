import numpy as np

from src.rgs_q6_gauged_full_clifford import (
    audit_q6_gauged_full_clifford,
    covariant_full_cubical_differential,
    gauged_full_dirac_kaehler_operator,
)
from src.rgs_q6_clifford_module import build_full_q6_cubical_complex
from src.rgs_q6_su2_connection import random_su2_matrices


def test_full_cubical_gauged_differential_has_all_degrees_and_is_hermitian():
    cube = build_full_q6_cubical_complex()
    links = random_su2_matrices(192, seed=201)
    d = covariant_full_cubical_differential(links, cube=cube)
    dk = gauged_full_dirac_kaehler_operator(links, cube=cube)
    assert d.shape == (1458, 1458)
    assert dk.shape == (1458, 1458)
    assert np.max(np.abs((dk - dk.getH()).data), initial=0.0) < 1e-12


def test_flatness_and_curvature_obstruction_to_covariant_d_squared():
    cube = build_full_q6_cubical_complex()
    identity = np.repeat(np.eye(2, dtype=complex)[None, :, :], 192, axis=0)
    flat = covariant_full_cubical_differential(identity, cube=cube)
    assert (flat @ flat).nnz == 0

    curved = covariant_full_cubical_differential(
        random_su2_matrices(192, seed=202), cube=cube
    )
    assert np.max(np.abs((curved @ curved).data), initial=0.0) > 1e-8


def test_gauged_full_clifford_transport_audit():
    report = audit_q6_gauged_full_clifford(seed=20261104)
    assert report["all_checks_passed"], report
    assert report["cell_counts_by_degree"] == [64, 192, 240, 160, 60, 12, 1]
    assert report["square_holonomies"] == 240
    assert report["gauge_covariance_max_abs_residual"] < 1e-10
    assert report["flat_connection_D_squared_max_abs"] < 1e-12
    assert report["curved_connection_D_squared_max_abs"] > 1e-8
