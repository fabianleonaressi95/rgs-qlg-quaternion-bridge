import numpy as np
import pytest

from src.rgs_q6_d5_a5_angular_measure import (
    a5_invariant_harmonic_multiplicity,
    a5_invariant_harmonic_spectrum,
    build_q6_d5_a5_angular_measure_audit,
    so3_spherical_harmonic_character,
)


def test_a5_character_average_has_no_low_degree_angular_harmonics():
    assert a5_invariant_harmonic_spectrum(6) == [1, 0, 0, 0, 0, 0, 1]
    assert a5_invariant_harmonic_multiplicity(0) == 1
    assert a5_invariant_harmonic_multiplicity(6) == 1


def test_so3_spherical_harmonic_character_identity_and_validation():
    assert so3_spherical_harmonic_character(0, 0.7) == pytest.approx(1.0)
    assert so3_spherical_harmonic_character(3, 0.0) == pytest.approx(7.0)
    with pytest.raises(ValueError):
        so3_spherical_harmonic_character(-1, 0.2)
    with pytest.raises(ValueError):
        a5_invariant_harmonic_spectrum(-1)


def test_q6_d5_phase_shell_pushforward_has_a5_symmetry_theorem():
    report = build_q6_d5_a5_angular_measure_audit(max_degree=8)
    assert report["A5_invariant_spherical_harmonic_multiplicities"][:7] == [
        1, 0, 0, 0, 0, 0, 1
    ]
    assert report["first_possible_angular_anisotropy_degree"] == 6
    assert report["degree_six_invariant_dimension"] == 1
    assert report["zero_locus_witness_matrix_rank"] == 3
    assert report["zero_locus_has_surface_measure_zero"]
    assert report["all_checks_passed"]
    assert "surface measure zero" in report["theorem"]
    assert any("zero locus" in item for item in report["conditional_or_open"])
