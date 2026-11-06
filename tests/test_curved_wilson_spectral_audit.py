import numpy as np
import pytest
from src.curved_wilson_spectral_audit import curved_graph_wilson_spectral_audit


def _rot_z(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0.0],
                     [np.sin(theta), np.cos(theta), 0.0],
                     [0.0, 0.0, 1.0]])


def test_curved_graph_wilson_spectral_audit_reports_finite_checks():
    x = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., .2]])
    e = np.stack([_rot_z(t) @ np.diag([1.0 + .05*i, 1.0, 1.0])
                  for i, t in enumerate([0.0, .08, .16, .23])])
    edges = [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2)]
    report = curved_graph_wilson_spectral_audit(
        x, e, edges, wilson_r=0.7, physical_cutoff=1e-7
    )
    assert report["scope"] == "FINITE_CURVED_GRAPH_ONLY"
    assert report["matrix_dimension"] == 16
    assert report["wilson_kernel"]["minimum_eigenvalue"] > -1e-9
    assert report["euclidean_dirac"]["gamma5_hermiticity_residual"] < 1e-9
    assert report["all_finite_checks_passed"]
    assert "no curved continuum" in " ".join(report["limitations"])


def test_curved_graph_audit_rejects_invalid_cutoff():
    x = np.array([[0., 0., 0.], [1., 0., 0.]])
    e = np.broadcast_to(np.eye(3), (2, 3, 3)).copy()
    with pytest.raises(ValueError, match="physical_cutoff"):
        curved_graph_wilson_spectral_audit(x, e, [(0, 1)], physical_cutoff=-1)
