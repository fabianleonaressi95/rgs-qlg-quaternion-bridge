import numpy as np

from src.overlap_dirac import assemble_overlap_from_wilson


def test_finite_overlap_operator_satisfies_ginsparg_wilson_relation():
    gamma5=np.diag([1.0,-1.0]).astype(complex)
    naive=np.array([[0.0,1.0],[1.0,0.0]],dtype=complex)
    wilson=0.3*np.eye(2,dtype=complex)
    overlap,audit=assemble_overlap_from_wilson(
        naive,wilson,gamma5,lattice_spacing=1.0,rho=1.0
    )
    assert overlap.shape==(2,2)
    assert audit["wilson_kernel_hermiticity_residual"]<1e-12
    assert audit["wilson_kernel_minimum_absolute_eigenvalue"]>0.0
    assert audit["ginsparg_wilson_residual"]<1e-12
    assert audit["gamma5_hermiticity_residual"]<1e-12


def test_overlap_rejects_kernel_without_gamma5_hermiticity():
    gamma5=np.diag([1.0,-1.0]).astype(complex)
    naive=np.eye(2,dtype=complex)
    wilson=0.3*np.eye(2,dtype=complex)
    try:
        assemble_overlap_from_wilson(naive,wilson,gamma5)
    except ValueError as exc:
        assert "requires" in str(exc)
    else:
        raise AssertionError("invalid naive operator must be rejected")


def test_overlap_index_trace_resolves_a_single_negative_chirality_zero_mode():
    gamma5=np.diag([1.0,-1.0]).astype(complex)
    naive=np.zeros((2,2),dtype=complex)
    wilson=np.diag([1.5,0.5]).astype(complex)
    overlap,audit=assemble_overlap_from_wilson(
        naive,wilson,gamma5,lattice_spacing=1.0,rho=1.0
    )
    from src.overlap_dirac import overlap_index
    result=overlap_index(overlap,gamma5,1.0,1.0)
    assert result["nearest_integer_candidate"]==-1
    assert result["integer_within_tolerance"]
    assert result["zero_mode_count"]==1
    assert np.isclose(result["zero_mode_chirality_sum"],-1.0,atol=1e-12)


def test_overlap_spectral_flow_detects_sampled_sign_count_change():
    from src.overlap_dirac import overlap_kernel_spectral_flow
    gamma5=np.diag([1.0,-1.0]).astype(complex)
    naive=np.zeros((2,2),dtype=complex)
    wilson=np.diag([1.5,0.5]).astype(complex)
    result=overlap_kernel_spectral_flow(naive,wilson,gamma5,1.0,[0.25,0.75,1.25,1.75])
    assert len(result["samples"])==4
    assert len(result["sampled_sign_count_changes"])>=1


def test_overlap_locality_profile_uses_graph_shortest_path_shells():
    from src.overlap_dirac import overlap_locality_profile
    d=np.zeros((12,12),dtype=complex)
    for i in range(3):
        for j in range(3):
            d[4*i:4*i+4,4*j:4*j+4]=np.eye(4)/(1+abs(i-j))
    result=overlap_locality_profile(d,[(0,1),(1,2)],3,4)
    assert set(result["shells"])=={"1","2"}
    assert result["status"].startswith("FINITE_GRAPH_DIAGNOSTIC")
