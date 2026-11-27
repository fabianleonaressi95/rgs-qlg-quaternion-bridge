import numpy as np

from scripts.discrete_cartan_core import (
    assert_lorentz,
    face_holonomy,
    face_torsion,
    gauge_transform,
    lorentz_boost_x,
    minkowski_metric,
    tetrahedron_boundary_bianchi,
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
    H2 = face_holonomy(face, links2)
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
    H = face_holonomy((0, 1, 2), links)
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

    H = face_holonomy((0, 1, 2), links)
    T = face_torsion((0, 1, 2), links, coframes)
    assert not np.allclose(H, np.eye(4))
    assert np.allclose(T, 0.0, atol=1e-12)



def lorentz_boost_y(rapidity):
    c = float(np.cosh(rapidity))
    s = float(np.sinh(rapidity))
    return np.array(
        [[c, 0.0, s, 0.0],
         [0.0, 1.0, 0.0, 0.0],
         [s, 0.0, c, 0.0],
         [0.0, 0.0, 0.0, 1.0]],
        dtype=float,
    )


def test_tetrahedron_affine_boundary_bianchi_identity():
    # Store all six positively oriented edges of the tetrahedron.
    # Noncommuting boosts ensure the test is not restricted to an Abelian subgroup.
    links = {
        (0, 1): lorentz_boost_x(0.17),
        (0, 2): lorentz_boost_y(-0.23),
        (0, 3): lorentz_boost_x(-0.11) @ lorentz_boost_y(0.19),
        (1, 2): lorentz_boost_y(0.31) @ lorentz_boost_x(-0.07),
        (1, 3): lorentz_boost_x(0.13) @ lorentz_boost_y(-0.29),
        (2, 3): lorentz_boost_y(0.09) @ lorentz_boost_x(0.21),
    }
    rng = np.random.default_rng(20261127)
    coframes = {
        edge: rng.normal(size=4)
        for edge in links
    }

    result = tetrahedron_boundary_bianchi(links, coframes, atol=2e-10)
    assert result["curvature_norm"] < 2e-10
    assert result["torsion_norm"] < 2e-10
    assert result["expanded_torsion_norm"] < 2e-10
    assert result["passed"] == 1.0



def test_all_five_tetrahedral_facets_of_a_4_simplex():
    # A 4-simplex has five tetrahedral boundary cells. Relabel each facet
    # to the local [0,1,2,3] convention and audit its affine boundary identity.
    vertices = range(5)
    all_edges = [(i, j) for i in vertices for j in vertices if i < j]
    rng = np.random.default_rng(20261128)
    links = {}
    coframes = {}
    for edge in all_edges:
        links[edge] = lorentz_boost_x(float(rng.uniform(-0.3, 0.3))) @ lorentz_boost_y(
            float(rng.uniform(-0.3, 0.3))
        )
        coframes[edge] = rng.normal(size=4)

    for omitted in vertices:
        facet = [v for v in vertices if v != omitted]
        relabel = {old: new for new, old in enumerate(facet)}
        facet_links = {
            (relabel[i], relabel[j]): value
            for (i, j), value in links.items()
            if i in relabel and j in relabel
        }
        facet_coframes = {
            (relabel[i], relabel[j]): value
            for (i, j), value in coframes.items()
            if i in relabel and j in relabel
        }
        result = tetrahedron_boundary_bianchi(
            facet_links, facet_coframes, atol=5e-10
        )
        assert result["curvature_norm"] < 5e-10
        assert result["torsion_norm"] < 5e-10
        assert result["passed"] == 1.0



def test_tetrahedron_affine_bianchi_is_local_lorentz_gauge_invariant():
    links = {
        (0, 1): lorentz_boost_x(0.17),
        (0, 2): lorentz_boost_y(-0.23),
        (0, 3): lorentz_boost_x(-0.11) @ lorentz_boost_y(0.19),
        (1, 2): lorentz_boost_y(0.31) @ lorentz_boost_x(-0.07),
        (1, 3): lorentz_boost_x(0.13) @ lorentz_boost_y(-0.29),
        (2, 3): lorentz_boost_y(0.09) @ lorentz_boost_x(0.21),
    }
    rng = np.random.default_rng(20261129)
    coframes = {edge: rng.normal(size=4) for edge in links}
    frames = {
        0: lorentz_boost_x(0.12),
        1: lorentz_boost_y(-0.18),
        2: lorentz_boost_x(-0.09) @ lorentz_boost_y(0.07),
        3: lorentz_boost_y(0.15),
    }

    original = tetrahedron_boundary_bianchi(links, coframes)
    links_g, coframes_g = gauge_transform(links, coframes, frames)
    transformed = tetrahedron_boundary_bianchi(links_g, coframes_g)

    assert original["curvature_norm"] < 2e-10
    assert original["torsion_norm"] < 2e-10
    assert transformed["curvature_norm"] < 2e-10
    assert transformed["torsion_norm"] < 2e-10
