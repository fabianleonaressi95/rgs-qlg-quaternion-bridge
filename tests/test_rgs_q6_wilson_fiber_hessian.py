import numpy as np

from src.rgs_q6_wilson_fiber_hessian import q6_wilson_fiber_hessian_report


def test_q6_wilson_quadratic_fibre_marginalization():
    for beta, lam in ((0.25, 12.0), (1.0, 12.0), (3.5, 12.0)):
        report = q6_wilson_fiber_hessian_report(beta, lam)
        assert report["all_checks_passed"], report["checks"]
        assert report["integrated_fibre_dimension"] == 124
        assert report["cross_hessian_frobenius"] < 1e-8
        assert report["schur_residual"] < 1e-8
        assert np.allclose(
            np.asarray(report["schur_complement"]),
            beta * lam * np.eye(5),
            atol=1e-8,
        )
