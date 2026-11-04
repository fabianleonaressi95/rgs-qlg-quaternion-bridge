import numpy as np

from src.rgs_q6_a4_stf_intertwiner import build_q6_a4_stf_intertwiner


def test_q6_a4_stf_linear_intertwiner():
    report = build_q6_a4_stf_intertwiner()
    assert report["a4_order"] == 12
    assert report["intertwiner_nullity"] == 3
    assert report["all_checks_passed"]
    assert np.linalg.matrix_rank(report["intertwiner"]) == 5
    assert report["equivariance_residual"] < 2e-7
    assert report["orthogonality_residual"] < 2e-7
