from src.gec6_trigram_decomposition import (
    LOWER, UPPER, audit, block_permutations, pair_sectors, wedge_image
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


def test_antipodal_audit():
    result = audit()
    assert result["antipodal_vertex_checks"] == 64
    assert result["antipodal_coordinate_pair_checks"] == 64 * 15


def test_bit_translation_orientation_sign_and_antipodal_triviality():
    from itertools import combinations
    from src.gec6_trigram_decomposition import translation_square_sign

    for translation in range(64):
        for i, j in combinations(range(6), 2):
            expected = -1 if (((translation >> i) & 1) + ((translation >> j) & 1)) % 2 else 1
            assert translation_square_sign(translation, i, j) == expected
    assert all(translation_square_sign(63, i, j) == 1 for i, j in combinations(range(6), 2))
    assert translation_square_sign(1, 0, 1) == -1
