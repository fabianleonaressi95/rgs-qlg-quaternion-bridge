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
