import numpy as np

from scripts.discrete_cartan_core import (
    assert_lorentz,
    face_holonomy,
    face_torsion,
    gauge_transform,
    lorentz_boost_x,
    minkowski_metric,
)


def triangle_data():
    # One stored orientation for each edge of (0,1,2).
    links = {
        (0, 1): lorentz_boost_x(0.13),
        (1, 2): lorentz_boost_x(-0.21),
        (0, 2): lorentz_boost_x(0.07),
    }
    coframes = {
        (0, 1): np.array([0.4, 0.1, 0.2, -0.1]),
        (1, 2): np.array([0.2, -0.3, 0.1, 0.2]),
        (0, 2): np.array([0.1, 0.2, -0.2, 0.3]),
    }
    return links, coframes


def test_boosts_preserve_lorentz_metric():
    for rapidity in (-1.1, -0.2, 0.0, 0.4, 1.3):
        assert_lorentz(lorentz_boost_x(rapidity))


def test_face_holonomy_and_torsion_transform_covariantly():
    links, coframes = triangle_data()
    face = (0, 1, 2)
    H = face_holonomy(face, links)
    T = face_torsion(face, links, coframes)

    frames = {
        0: lorentz_boost_x(0.31),
        1: lorentz_boost_x(-0.17),
        2: lorentz_boost_x(0.23),
    }
    links2, coframes2 = gauge_transform(links, coframes, frames)
    H2 = face_holonomy(face, links2, coframes2)
    T2 = face_torsion(face, links2, coframes2)

    h0 = frames[0]
    assert np.allclose(H2, h0 @ H @ np.linalg.inv(h0), atol=1e-11)
    assert np.allclose(T2, h0 @ T, atol=1e-11)


def test_exact_coframe_closes_when_connection_is_identity():
    points = {
        0: np.array([0.0, 0.0, 0.0, 0.0]),
        1: np.array([1.0, 0.2, 0.0, 0.0]),
        2: np.array([0.3, 1.0, 0.1, 0.0]),
    }
    I = np.eye(4)
    links = {(0, 1): I, (1, 2): I, (0, 2): I}
    coframes = {
        (0, 1): points[1] - points[0],
        (1, 2): points[2] - points[1],
        (0, 2): points[2] - points[0],
    }
    assert np.allclose(face_torsion((0, 1, 2), links, coframes), 0.0, atol=1e-14)
    assert np.allclose(face_holonomy((0, 1, 2), links), I, atol=1e-14)


def test_flat_connection_can_have_nonzero_coframe_closure():
    I = np.eye(4)
    links = {(0, 1): I, (1, 2): I, (0, 2): I}
    coframes = {
        (0, 1): np.array([1.0, 0.0, 0.0, 0.0]),
        (1, 2): np.array([0.0, 1.0, 0.0, 0.0]),
        (0, 2): np.array([0.0, 0.0, 0.0, 0.0]),
    }
    H = face_holonomy((0, 1, 2), links, coframes)
    T = face_torsion((0, 1, 2), links, coframes)
    assert np.allclose(H, I)
    assert np.linalg.norm(T) > 0.0


def test_curved_holonomy_can_be_paired_with_zero_face_closure():
    # Choose E_20 so that the transported sum vanishes for a nontrivial H.
    L01 = lorentz_boost_x(0.2)
    L12 = lorentz_boost_x(-0.1)
    L02 = lorentz_boost_x(0.05)
    links = {(0, 1): L01, (1, 2): L12, (0, 2): L02}
    E01 = np.array([0.2, 0.1, 0.0, 0.0])
    E12 = np.array([0.1, -0.2, 0.0, 0.0])
    L20 = np.linalg.inv(L02)
    E20 = -np.linalg.inv(L01 @ L12) @ (E01 + L01 @ E12)
    # Store E_02 by reversing E_20 = -L_02^-1 E_02.
    E02 = -L02 @ E20
    coframes = {(0, 1): E01, (1, 2): E12, (0, 2): E02}

    H = face_holonomy((0, 1, 2), links, coframes)
    T = face_torsion((0, 1, 2), links, coframes)
    assert not np.allclose(H, np.eye(4))
    assert np.allclose(T, 0.0, atol=1e-12)
