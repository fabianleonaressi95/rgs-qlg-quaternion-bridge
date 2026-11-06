import numpy as np
import pytest

from src.ellis_spin_connection import (
    ellis_edge_spin_link,
    ellis_edge_spin_links,
    ellis_orthonormal_frame,
)
from src.rgs_ellis_refinement import ellis_embedding, build_ellis_rgs_refinement


def _point(u, n, r0):
    return ellis_embedding(np.array([u]), np.asarray(n, dtype=float)[None, :], r0)[0]


def test_ellis_frame_is_orthonormal_and_tangent():
    r0 = 0.83
    n = np.array([1., 2., -1.])
    n /= np.linalg.norm(n)
    u = 0.37
    frame = ellis_orthonormal_frame(u, n, r0)
    R = np.sqrt(u*u+r0*r0)
    normal = np.r_[(r0/R)*n, -u/R]
    assert np.allclose(frame.T@frame, np.eye(3), atol=2e-12)
    assert np.linalg.norm(frame.T@normal) < 2e-12


def test_radial_parallel_transport_is_identity_in_adapted_frame():
    r0 = 0.8292
    n = np.array([1., 2., 3.])
    n /= np.linalg.norm(n)
    p = _point(-0.4, n, r0)
    q = _point(0.7, n, r0)
    link, audit = ellis_edge_spin_link(p, q, r0, steps=8)
    assert np.allclose(link, np.eye(2), atol=2e-11)
    assert audit["spin_unitarity_residual"] < 1e-12


def test_edge_transport_is_reversible():
    r0 = 0.8292
    n0 = np.array([1., 0., 0.])
    n1 = np.array([0.8, 0.6, 0.])
    p, q = _point(-0.2, n0, r0), _point(0.35, n1, r0)
    forward, _ = ellis_edge_spin_link(p, q, r0, steps=64)
    backward, _ = ellis_edge_spin_link(q, p, r0, steps=64)
    assert np.allclose(backward, forward.conj().T, atol=2e-9)


def test_spin_link_is_unitary_and_lift_matches_so3_rotation():
    from src.binary_icosahedral import quaternion_to_so3
    from src.discrete_dirac import _rotation_to_quaternion
    r0 = 0.8292
    p = _point(0.1, [1., 0., 0.], r0)
    n1 = np.array([0.9, 0.3, 0.2])
    n1 /= np.linalg.norm(n1)
    q = _point(0.2, n1, r0)
    link, audit = ellis_edge_spin_link(p, q, r0, steps=64)
    # The scalar-first lift must reproduce the endpoint frame rotation.
    assert np.allclose(
        quaternion_to_so3(np.asarray(audit["spin_quaternion"])),
        np.asarray(audit["so3_rotation"]), atol=2e-10
    )
    assert audit["so3_orthogonality_residual"] < 1e-10
    assert audit["spin_unitarity_residual"] < 1e-10


def test_ellis_mesh_links_are_built_on_refinement_edges():
    mesh = build_ellis_rgs_refinement(level=0, r0=0.8292, radial_base=2,
                                      include_chain=True)
    edges = mesh.chain["edges"]
    links, audit = ellis_edge_spin_links(mesh.vertices_r4, edges, mesh.r0, steps=8)
    assert len(links) == len(edges)
    assert audit["connection"] == "Levi-Civita of exact Ellis-RGS hypersurface"
    assert audit["max_transport_residual"] < 2e-9
    for link in links.values():
        assert np.allclose(link.conj().T@link, np.eye(2), atol=2e-10)


def test_antipodal_edge_is_rejected_as_ambiguous():
    r0 = 1.0
    p = _point(0.0, [1., 0., 0.], r0)
    q = _point(0.0, [-1., 0., 0.], r0)
    with pytest.raises(ValueError, match="antipodal"):
        ellis_edge_spin_link(p, q, r0)



def test_ellis_spin_holonomy_detects_nonzero_curvature_at_throat():
    r0 = 0.8292
    eps = 0.12
    directions = [
        np.array([1., 0., 0.]),
        np.array([np.cos(eps), np.sin(eps), 0.]),
        np.array([np.cos(eps), 0., np.sin(eps)]),
    ]
    points = [_point(0.0, n, r0) for n in directions]
    u01, _ = ellis_edge_spin_link(points[0], points[1], r0, steps=96)
    u12, _ = ellis_edge_spin_link(points[1], points[2], r0, steps=96)
    u20, _ = ellis_edge_spin_link(points[2], points[0], r0, steps=96)
    holonomy = u20 @ u12 @ u01
    distance_from_identity = float(np.linalg.norm(holonomy-np.eye(2), ord="fro"))
    assert distance_from_identity > 1e-5
    assert np.linalg.norm(holonomy.conj().T@holonomy-np.eye(2)) < 1e-10


def test_rk4_edge_transport_converges_under_step_refinement():
    r0 = 0.8292
    p = _point(-0.3, [1., 0., 0.], r0)
    q = _point(0.45, [0.7, 0.5, np.sqrt(0.26)], r0)
    u16, _ = ellis_edge_spin_link(p, q, r0, steps=16)
    u32, _ = ellis_edge_spin_link(p, q, r0, steps=32)
    u64, _ = ellis_edge_spin_link(p, q, r0, steps=64)
    e16 = np.linalg.norm(u16-u64, ord="fro")
    e32 = np.linalg.norm(u32-u64, ord="fro")
    assert e32 < e16
