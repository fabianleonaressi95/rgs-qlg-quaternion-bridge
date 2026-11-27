from itertools import permutations

from src.q6_prism_trigram_equivariant_map import (
    BASE, LENGTHS, audit, edge_path, edge_winding, q6_edges,
    selected_unimodular_square_images, square_winding, symmetry_audit,
)


def test_q6_edge_counts_and_based_path_assignment():
    edges = q6_edges()
    assert len(edges) == 192
    for v, _w, d in edges:
        path = edge_path(v, d)
        assert path[0] == BASE
        assert path[-1] == BASE
        for a, b in zip(path, path[1:]):
            changed = [i for i in range(5) if a[i] != b[i]]
            assert len(changed) == 1
            i = changed[0]
            assert (b[i] - a[i]) % LENGTHS[i] == 1


def test_five_square_images_are_unimodular():
    import sympy as sp
    columns = selected_unimodular_square_images()
    assert columns == [
        (0, -1, 0, 0, 0),
        (0, 0, -1, 0, 0),
        (0, 0, 0, -1, 0),
        (-1, -1, 0, 0, 0),
        (0, -1, 0, 0, -1),
    ]
    matrix = sp.Matrix.hstack(*(sp.Matrix(c) for c in columns))
    assert abs(int(matrix.det())) == 1


def test_square_winding_formula_for_all_upper_weights():
    for k in range(3):
        mask = (1 << k) - 1
        for d in range(3):
            image = square_winding(d, 3, mask)
            assert image[d + 1] == -1
            assert image[0] == (-1 if k == 1 else 0)
            assert image[4] == (-1 if k == 2 else 0)


def test_full_trigram_block_permutation_equivariance():
    result = symmetry_audit()
    assert result["source_group_order"] == 36
    assert result["equivariant"]
    assert result["edge_path_equivariance_checks"] == 36 * 192


def test_audit_claims_integral_surjection():
    result = audit()
    assert result["induced_H1_map_is_surjective"]
    assert abs(result["selected_minor_determinant"]) == 1
