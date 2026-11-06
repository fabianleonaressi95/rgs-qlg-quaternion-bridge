import numpy as np
from src.ellis_spin_connection import ellis_edge_spin_links
from src.ellis_rgs_spin_dirac import assemble_ellis_spin_wilson_dirac
from src.rgs_ellis_refinement import build_ellis_rgs_refinement

def small_mesh():
    return build_ellis_rgs_refinement(level=0,r0=0.8292,radial_base=2,include_chain=True)

def tetra_edges(mesh):
    tet=tuple(map(int,mesh.tetrahedra[0]))
    return [(min(tet[a],tet[b]),max(tet[a],tet[b])) for a in range(4) for b in range(a+1,4)]

def test_intrinsic_ellis_operator_gamma5_hermiticity_and_wilson_positivity():
    mesh=small_mesh(); edges=tetra_edges(mesh)
    links,_=ellis_edge_spin_links(mesh.vertices_r4,edges,mesh.r0,steps=12)
    op,a=assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links=links,wilson_r=.7)
    assert op.shape==(4*len(mesh.vertices_r4),)*2
    assert a["dirac_hermiticity_residual"]<2e-10
    assert a["wilson_hermiticity_residual"]<2e-10
    assert a["gamma5_hermiticity_residual"]<2e-9
    assert a["wilson_minimum_eigenvalue"]>-2e-9
    assert a["edge_metric"]=="intrinsic Ellis path length"

def test_intrinsic_ellis_operator_is_local_spin_frame_gauge_covariant():
    mesh=small_mesh(); edges=tetra_edges(mesh)
    links,_=ellis_edge_spin_links(mesh.vertices_r4,edges,mesh.r0,steps=12)
    rng=np.random.default_rng(20261106); gauges=[]
    for _ in mesh.vertices_r4:
        q=rng.normal(size=4); q/=np.linalg.norm(q)
        w,x,y,z=q
        gauges.append([[w-1j*z,-y-1j*x],[y-1j*x,w+1j*z]])
    gauges=np.asarray(gauges,dtype=complex)
    base,_=assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links=links,wilson_r=.4)
    changed,a=assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links=links,wilson_r=.4,local_spin_frames=gauges)
    G=a["gauge_covariance_transform"]
    assert np.allclose(changed,G@base@G.conj().T,atol=2e-9)

def test_intrinsic_ellis_zero_wilson_and_positive_wilson_are_distinct():
    mesh=small_mesh(); edges=tetra_edges(mesh)
    links,_=ellis_edge_spin_links(mesh.vertices_r4,edges,mesh.r0,steps=8)
    d0,a0=assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links=links,wilson_r=0)
    d1,a1=assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links=links,wilson_r=.8)
    assert a0["wilson_minimum_eigenvalue"]==0.0
    assert a1["wilson_minimum_eigenvalue"]>-1e-9
    assert np.linalg.norm(d1-d0)>0

def test_intrinsic_ellis_assembler_rejects_incomplete_links():
    mesh=small_mesh(); edges=tetra_edges(mesh)
    with np.testing.assert_raises_regex(ValueError,"keys mismatch"):
        assemble_ellis_spin_wilson_dirac(mesh,edges=edges,spin_links={edges[0]:np.eye(2)})
