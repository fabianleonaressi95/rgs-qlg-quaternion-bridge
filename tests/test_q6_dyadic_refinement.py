import numpy as np

from src.q6_complex import build_q6_complex
from src.q6_dyadic_refinement import build_q6_dyadic_refinement, audit_q6_dyadic_refinement


def test_q6_dyadic_refinement_cell_counts_and_map_shapes():
    r = build_q6_dyadic_refinement()
    assert (len(r.vertices), len(r.edges), len(r.faces)) == (496, 1344, 960)
    assert r.S0.shape == (496, 64)
    assert r.S1.shape == (1344, 192)
    assert r.S2.shape == (960, 240)


def test_q6_dyadic_refinement_is_an_exact_chain_map():
    c = build_q6_complex()
    f = build_q6_dyadic_refinement()
    assert np.array_equal(f.boundary_1 @ f.boundary_2, np.zeros((496, 960), dtype=np.int64))
    assert np.array_equal(f.boundary_1 @ f.S1, f.S0 @ c.boundary_1)
    assert np.array_equal(f.boundary_2 @ f.S2, f.S1 @ c.boundary_2)


def test_each_parent_cell_has_the_expected_oriented_children():
    f = build_q6_dyadic_refinement()
    assert np.all(np.sum(np.abs(f.S0), axis=0) == 1)
    assert np.all(np.sum(np.abs(f.S1), axis=0) == 2)
    assert np.all(np.sum(np.abs(f.S2), axis=0) == 4)
    assert set(np.unique(f.S0)) <= {0, 1}
    assert set(np.unique(f.S1)) <= {-1, 0, 1}
    assert set(np.unique(f.S2)) <= {0, 1}


def test_audit_reports_zero_chain_map_residuals():
    report = audit_q6_dyadic_refinement()
    assert report["fine_counts"] == [496, 1344, 960]
    assert report["all_exact"]
    assert report["chain_map_S1_residual_max_abs"] == 0
    assert report["chain_map_S2_residual_max_abs"] == 0
