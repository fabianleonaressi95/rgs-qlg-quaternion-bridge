from src.gec6_a5_spectral import (
    exact_fourier_spectrum_m3,
    fourier_laplacian_eigenvalue,
    numerical_matrix_spectrum,
    normalized_comparison,
)


def test_exact_fourier_laplacian_spectrum_m3():
    assert exact_fourier_spectrum_m3() == {
        0: 1,
        15: 12,
        24: 30,
        27: 50,
        33: 120,
        36: 30,
    }


def test_fourier_multiplicities_sum_to_number_of_vertices():
    assert sum(exact_fourier_spectrum_m3().values()) == 3**5


def test_trivial_character_gives_zero_laplacian_eigenvalue():
    assert fourier_laplacian_eigenvalue((0, 0, 0, 0, 0)) == 0


def test_numerical_matrix_diagonalization_matches_fourier():
    assert numerical_matrix_spectrum() == exact_fourier_spectrum_m3()


def test_normalized_gaps_are_recorded_correctly():
    comparison = normalized_comparison()
    assert comparison["GEC6_m3_normalized_spectral_gap"] == 0.5
    assert comparison["Q6_normalized_spectral_gap"] == 1 / 3
    assert 0 < comparison["PRISM_C64xC6xC6xC6xC5_normalized_spectral_gap"] < 0.001


def test_general_modulus_gap_values_m3_to_m8():
    from math import isclose
    from src.gec6_a5_spectral import normalized_spectral_gap
    expected = {
        3: 0.5,
        4: 1.0 / 3.0,
        5: 0.23032766854168418,
        6: 1.0 / 6.0,
        7: 0.12550339938042213,
        8: 0.09763107293781749,
    }
    for modulus, value in expected.items():
        assert isclose(normalized_spectral_gap(modulus), value, rel_tol=1e-11, abs_tol=1e-12)


def test_gap_scaling_audit_approaches_quadratic_limit():
    from math import pi
    from src.gec6_a5_spectral import gap_scaling_audit
    rows = gap_scaling_audit(8)["rows"]
    assert len(rows) == 6
    assert rows[-1]["m2_times_gap"] > rows[0]["m2_times_gap"]
    assert rows[-1]["m2_times_gap"] < 2 * pi * pi / 3
