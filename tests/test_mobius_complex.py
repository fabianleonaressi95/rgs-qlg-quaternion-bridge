"""Tests for the canonical six-square Möbius-band subcomplex."""
import numpy as np

from src.mobius_complex import audit_mobius_complex, build_mobius_complex


def test_mobius_cardinalities_and_cell_structure():
    audit = audit_mobius_complex()
    assert (audit["vertices"], audit["edges"], audit["faces"]) == (12, 18, 6)
    assert audit["boundary_edges"] == 12
    assert audit["interior_edges"] == 6
    assert audit["boundary_components"] == 1


def test_mobius_chain_complex_identity_is_exact():
    complex_ = build_mobius_complex()
    product = (
        complex_.boundary_1.astype(np.int64)
        @ complex_.boundary_2.astype(np.int64)
    )
    assert product.shape == (12, 6)
    assert np.count_nonzero(product) == 0


def test_mobius_euler_characteristic_and_betti_numbers():
    audit = audit_mobius_complex()
    assert audit["euler_characteristic"] == 0
    assert audit["rank_boundary_1"] == 11
    assert audit["rank_boundary_2"] == 6
    assert (audit["betti_0"], audit["betti_1"], audit["betti_2"]) == (1, 1, 0)


def test_orientation_propagation_detects_nonorientability():
    audit = audit_mobius_complex()
    assert audit["orientable"] is False
    assert audit["orientation_conflicts"] > 0


def test_faces_have_cyclic_order_and_four_distinct_vertices():
    complex_ = build_mobius_complex()
    assert all(len(set(face)) == 4 for face in complex_.faces)
    assert complex_.faces[0] == (12, 44, 60, 28)


def test_mobius_vertices_retain_big_endian_q6_labels():
    complex_ = build_mobius_complex()
    assert complex_.vertices == (0, 4, 8, 12, 16, 24, 28, 40, 44, 48, 56, 60)


def test_incidence_matrices_are_integer_valued():
    complex_ = build_mobius_complex()
    assert np.issubdtype(complex_.boundary_1.dtype, np.integer)
    assert np.issubdtype(complex_.boundary_2.dtype, np.integer)
