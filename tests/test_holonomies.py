
import numpy as np
import pytest

from src.holonomies import (
    SU2Lattice,
    dagger,
    is_su2,
)


def su2_rotation(axis, angle):
    """Build an SU(2) matrix using the Pauli-matrix convention."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)

    x, y, z = axis
    c = np.cos(angle / 2.0)
    s = np.sin(angle / 2.0)

    return np.array(
        [
            [c - 1j * s * z, -s * y - 1j * s * x],
            [s * y - 1j * s * x, c + 1j * s * z],
        ],
        dtype=complex,
    )


def make_triangle():
    """Construct a three-vertex lattice with nontrivial SU(2) links."""
    lattice = SU2Lattice([0, 1, 2])

    u01 = su2_rotation([1, 0, 0], 0.37)
    u12 = su2_rotation([0, 1, 0], 0.61)
    u20 = su2_rotation([0, 0, 1], 0.83)

    lattice.add_link(0, 1, u01)
    lattice.add_link(1, 2, u12)
    lattice.add_link(2, 0, u20)

    return lattice, u01, u12, u20


def test_link_matrices_are_in_su2():
    lattice, u01, u12, u20 = make_triangle()

    assert is_su2(u01)
    assert is_su2(u12)
    assert is_su2(u20)

    for source, target in lattice.edges():
        assert is_su2(lattice.link(source, target))
        assert is_su2(lattice.link(target, source))


def test_reverse_link_is_hermitian_conjugate():
    lattice, _, _, _ = make_triangle()

    for source, target in lattice.edges():
        forward = lattice.link(source, target)
        reverse = lattice.link(target, source)

        assert np.allclose(reverse, dagger(forward), atol=1e-12)
        assert np.allclose(
            forward @ reverse,
            np.eye(2),
            atol=1e-12,
        )


def test_path_transport_uses_column_vector_order():
    lattice, u01, u12, _ = make_triangle()

    expected = u12 @ u01
    actual = lattice.path_transport([0, 1, 2])

    assert np.allclose(actual, expected, atol=1e-12)


def test_reverse_path_transport_is_inverse():
    lattice, _, _, _ = make_triangle()

    forward = lattice.path_transport([0, 1, 2])
    reverse = lattice.path_transport([2, 1, 0])

    assert np.allclose(reverse, dagger(forward), atol=1e-12)


def test_closed_path_holonomy_is_in_su2():
    lattice, _, _, _ = make_triangle()

    holonomy = lattice.holonomy([0, 1, 2, 0])

    assert is_su2(holonomy)


def test_wilson_trace_is_real_and_normalized():
    lattice, _, _, _ = make_triangle()

    trace = lattice.wilson_trace([0, 1, 2, 0])
    full_trace = lattice.wilson_trace(
        [0, 1, 2, 0],
        normalized=False,
    )

    assert np.isclose(full_trace, 2.0 * trace, atol=1e-12)
    assert np.isclose(trace.imag if isinstance(trace, complex) else 0.0, 0.0)
    assert -1.0 - 1e-12 <= trace <= 1.0 + 1e-12


def test_gauge_transformation_preserves_su2_links():
    lattice, _, _, _ = make_triangle()

    gauge = {
        0: su2_rotation([1, 1, 0], 0.23),
        1: su2_rotation([0, 1, 1], 0.47),
        2: su2_rotation([1, 0, 1], 0.71),
    }

    transformed = lattice.gauge_transform(gauge)

    for source, target in transformed.edges():
        assert is_su2(transformed.link(source, target))


def test_closed_holonomy_transforms_by_conjugation():
    lattice, _, _, _ = make_triangle()

    gauge = {
        0: su2_rotation([1, 1, 0], 0.23),
        1: su2_rotation([0, 1, 1], 0.47),
        2: su2_rotation([1, 0, 1], 0.71),
    }

    original_holonomy = lattice.holonomy([0, 1, 2, 0])
    transformed = lattice.gauge_transform(gauge)
    transformed_holonomy = transformed.holonomy([0, 1, 2, 0])

    expected = (
        gauge[0]
        @ original_holonomy
        @ dagger(gauge[0])
    )

    assert np.allclose(
        transformed_holonomy,
        expected,
        atol=1e-11,
    )


def test_wilson_trace_is_gauge_invariant():
    lattice, _, _, _ = make_triangle()

    gauge = {
        0: su2_rotation([1, 1, 0], 0.23),
        1: su2_rotation([0, 1, 1], 0.47),
        2: su2_rotation([1, 0, 1], 0.71),
    }

    path = [0, 1, 2, 0]

    original = lattice.wilson_trace(path)
    transformed = lattice.gauge_transform(gauge).wilson_trace(path)

    assert np.isclose(original, transformed, atol=1e-11)


def test_gauge_transform_requires_all_vertices():
    lattice, _, _, _ = make_triangle()

    incomplete_gauge = {
        0: np.eye(2, dtype=complex),
        1: np.eye(2, dtype=complex),
    }

    with pytest.raises(ValueError):
        lattice.gauge_transform(incomplete_gauge)


def test_invalid_link_matrix_is_rejected():
    lattice = SU2Lattice([0, 1])

    invalid = np.array(
        [
            [1.0, 0.0],
            [0.0, 2.0],
        ],
        dtype=complex,
    )

    with pytest.raises(ValueError):
        lattice.add_link(0, 1, invalid)


def test_non_closed_path_is_rejected_as_holonomy():
    lattice, _, _, _ = make_triangle()

    with pytest.raises(ValueError):
        lattice.holonomy([0, 1, 2])


def test_path_with_unknown_vertex_is_rejected():
    lattice, _, _, _ = make_triangle()

    with pytest.raises(KeyError):
        lattice.path_transport([0, 1, 99])
```
