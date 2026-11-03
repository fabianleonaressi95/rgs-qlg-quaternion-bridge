"""Structural and incidence tests for the canonical Q6 complex."""
import numpy as np
import pytest

from src.q6_complex import audit_q6_complex, build_q6_complex


def test_q6_cardinalities_and_regular_degree():
    audit = audit_q6_complex()
    assert audit["vertices"] == 64
    assert audit["edges"] == 192
    assert audit["faces"] == 240
    assert audit["degree_min"] == 6
    assert audit["degree_max"] == 6


def test_q6_chain_complex_identity_is_exact():
    q6 = build_q6_complex()
    product = (
        q6.boundary_1.astype(np.int64)
        @ q6.boundary_2.astype(np.int64)
    )
    assert product.shape == (64, 240)
    assert np.count_nonzero(product) == 0


def test_q6_boundary_matrices_have_expected_shapes():
    q6 = build_q6_complex()
    assert q6.boundary_1.shape == (64, 192)
    assert q6.boundary_2.shape == (192, 240)


def test_every_square_face_has_four_distinct_vertices():
    q6 = build_q6_complex()
    assert all(len(set(face)) == 4 for face in q6.faces)


def test_big_endian_axis_convention_is_fixed():
    q6 = build_q6_complex()
    # Axis 0 flips the most significant bit; axis 5 flips the least.
    assert (0, 32) in q6.edges
    assert (0, 1) in q6.edges
    assert (0, 2) not in q6.edges


def test_incidence_matrices_are_integer_valued():
    q6 = build_q6_complex()
    assert np.issubdtype(q6.boundary_1.dtype, np.integer)
    assert np.issubdtype(q6.boundary_2.dtype, np.integer)
