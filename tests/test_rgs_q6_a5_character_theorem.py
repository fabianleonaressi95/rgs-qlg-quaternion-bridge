from src.rgs_q6_a5_character_theorem import q6_top_coexact_a5_character_report


def test_q6_top_coexact_five_modes_are_irreducible_a5_five():
    report = q6_top_coexact_a5_character_report()
    assert report["group_order"] == 60
    assert report["lambda12_full_one_cochain_dimension"] == 6
    assert report["lambda12_exact_one_cochain_dimension"] == 1
    assert report["lambda12_coexact_dimension"] == 5
    assert report["all_checks_passed"], report
    assert report["coexact_lambda12_character"] == __import__("pytest").approx(
        [5.0, 1.0, -1.0, 0.0, 0.0], abs=2e-7
    )
    assert report["character_inner_product"] == __import__("pytest").approx(
        1.0, abs=2e-7
    )
    assert report["irreducible_multiplicities"]["5"] == __import__("pytest").approx(
        1.0, abs=2e-7
    )
