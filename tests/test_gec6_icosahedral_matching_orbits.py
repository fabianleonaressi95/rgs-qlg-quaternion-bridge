from src.gec6_icosahedral_matching_orbits import audit_icosahedral_matching_orbits


def test_a5_orbits_on_pairs_partitions_and_perfect_matchings():
    result = audit_icosahedral_matching_orbits()
    assert result["group_order"] == 60
    assert result["pair_orbit_sizes"] == [15]
    assert result["three_plus_three_partition_orbit_sizes"] == [10]
    assert result["perfect_matching_orbit_sizes"] == [5, 10]


def test_c6_tensor_selects_the_five_matching_orbit_not_one_matching():
    result = audit_icosahedral_matching_orbits()
    assert result["special_matching_count_C6_product_plus"] == 5
    assert result["generic_matching_count_C6_product_minus"] == 10
    assert result["C6_product_classes_equal_A5_matching_orbits"]
    assert result["diagonal_matching_in_special_orbit"]
    assert result["diagonal_matching_orbit_size"] == 5


def test_each_trigram_partition_has_two_special_cross_matchings():
    result = audit_icosahedral_matching_orbits()
    assert result["special_cross_matchings_per_three_plus_three_partition"] == [2]
    assert len(result["special_cross_matchings_for_lower_upper_partition"]) == 2


def test_full_signed_a5_does_not_descend_to_gec6_quotient():
    result = audit_icosahedral_matching_orbits()
    assert result["block_partition_stabilizer_order"] == 6
    assert result["GEC6_contraction_kernel_preserving_signed_actions"] == 6
    assert result["quotient_preserving_group_equals_block_partition_stabilizer"]
    assert result["diagonal_matching_orbit_under_quotient_preserving_subgroup"] == 2
