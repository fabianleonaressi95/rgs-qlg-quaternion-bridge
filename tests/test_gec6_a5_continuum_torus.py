from src.gec6_a5_continuum_torus import (
    Q, quadratic_form, shortest_dual_vectors,
    shortest_vector_certificate, continuum_spectral_data,
    de_aquino_comparison,
)


def test_q_is_symmetric_and_shortest_value_is_five():
    assert all(Q[i][j] == Q[j][i] for i in range(5) for j in range(5))
    assert quadratic_form((1, 0, 0, 0, 0)) == 5
    assert shortest_vector_certificate()["minimum"] == 5


def test_first_eigenvalue_has_twelve_shortest_vectors():
    assert len(shortest_dual_vectors()) == 12
    data = continuum_spectral_data()
    assert data["first_eigenvalue"] == "10*pi^2/3"
    assert data["first_eigenvalue_multiplicity"] == 12


def test_de_aquino_bridge_is_not_overstated():
    data = de_aquino_comparison()
    assert data["zero_threshold_squared"] == "5/4"
    assert data["candidate_normalization"] == "xi_k^2 = (k^T Q k)/4"
    assert data["status"].startswith("CONJECTURAL NORMALIZATION ONLY")
