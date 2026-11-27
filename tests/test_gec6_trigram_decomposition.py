from src.gec6_trigram_decomposition import (
    LOWER, UPPER, audit, block_permutations, pair_sectors, wedge_image,
    mixed_pair_orbits_independent, mixed_pair_orbits_diagonal
)


def test_pair_sector_partition_has_dimensions_3_9_3():
    lower, mixed, upper = pair_sectors()
    assert (len(lower), len(mixed), len(upper)) == (3, 9, 3)
    assert len(set(lower) | set(mixed) | set(upper)) == 15


def test_block_wreath_group_has_order_72_and_preserves_sectors():
    group = block_permutations()
    assert len(group) == 72
    lower, mixed, upper = pair_sectors()
    internal = set(lower) | set(upper)
    cross = set(mixed)
    for sigma in group:
        assert all(wedge_image(pair, sigma)[1] in internal for pair in lower + upper)
        assert all(wedge_image(pair, sigma)[1] in cross for pair in mixed)


def test_mixed_sector_orbits_depend_on_symmetry_choice():
    # Independent relabelling of lower and upper coordinates is transitive.
    assert tuple(map(len, mixed_pair_orbits_independent())) == (9,)
    # A chosen lower/upper coordinate matching reduces symmetry to diagonal S3.
    assert tuple(map(len, mixed_pair_orbits_diagonal())) == (3, 6)


def test_antipodal_audit():
    result = audit()
    assert result["antipodal_vertex_checks"] == 64
    assert result["antipodal_coordinate_pair_checks"] == 64 * 15
