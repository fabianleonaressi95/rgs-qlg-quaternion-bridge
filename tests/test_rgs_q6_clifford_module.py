import numpy as np

from src.rgs_q6_clifford_module import (
    audit_q6_clifford_module, build_full_q6_cubical_complex,
    exterior_clifford_generators, form_degree_chirality,
)


def test_full_q6_cube_has_all_cell_dimensions_and_boundary_squared_zero():
    c = build_full_q6_cubical_complex()
    assert c["cell_counts"] == (64, 192, 240, 160, 60, 12, 1)
    # The global cochain differential is degree-raising and has the
    # expected rectangular block orientation (C^k -> C^(k+1)).
    assert c["d"].shape == (sum(c["cell_counts"]), sum(c["cell_counts"]))
    for k in range(1, 6):
        assert np.array_equal(
            c["coboundaries"][k+1] @ c["coboundaries"][k],
            np.zeros((c["cell_counts"][k+1], c["cell_counts"][k-1]), dtype=np.int8),
        )


def test_six_dimensional_exterior_clifford_relations_and_chirality():
    gammas = exterior_clifford_generators(6)
    chi = form_degree_chirality(6)
    eye = np.eye(64)
    for i, gi in enumerate(gammas):
        assert np.allclose(gi, gi.conj().T)
        assert np.allclose(chi @ gi + gi @ chi, 0)
        for j, gj in enumerate(gammas):
            assert np.allclose(gi @ gj + gj @ gi, 2*eye if i == j else 0)


def test_complete_clifford_audit():
    report = audit_q6_clifford_module()
    assert report["all_checks_passed"]
    assert report["boundary_ranks"] == [63, 129, 111, 49, 11, 1]
