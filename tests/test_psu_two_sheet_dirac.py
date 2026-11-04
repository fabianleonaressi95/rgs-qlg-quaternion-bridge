import numpy as np
import pytest
from src.psu_two_sheet_dirac import assemble_two_sheet_dirac, symmetric_sheet_spectrum, two_sheet_spectrum

def test_identical_sheets_split_exactly():
    d=np.array([[0.,1.],[1.,0.]],complex)
    h,a=assemble_two_sheet_dirac(d,coupling=.37)
    assert np.allclose(two_sheet_spectrum(h)["eigenvalues"],symmetric_sheet_spectrum(d,.37),atol=1e-12)
    assert a["hermiticity_residual"]<1e-12 and a["total_dimension"]==4

def test_independent_sheets_and_unitary_transporter():
    dp=np.diag([-2.,.5,3.]).astype(complex); dm=np.diag([-1.,1.5,2.]).astype(complex)
    w=np.diag(np.exp(.23j*np.arange(3)))
    h,a=assemble_two_sheet_dirac(dp,dm,w,coupling=.4,mass_plus=.1,mass_minus=.2)
    assert np.allclose(h,h.conj().T,atol=1e-12)
    assert a["transporter_unitarity_residual"]<1e-12
    assert np.all(np.isfinite(two_sheet_spectrum(h)["eigenvalues"]))

def test_nonunitary_transporter_rejected():
    with pytest.raises(ValueError,match="unitary"):
        assemble_two_sheet_dirac(np.eye(2),transporter=2*np.eye(2),coupling=.1)

def test_scalar_sheet_coupling_is_mass_like_for_chirality():
    # A valid massless one-sheet Dirac block anticommutes with gamma5.
    # Use sigma_x against sigma_z; identical diagonal matrices commute.
    d=np.array([[0.,1.],[1.,0.]],complex)
    g5=np.diag([1.,-1.]).astype(complex)
    h,a=assemble_two_sheet_dirac(d,coupling=.2,gamma5=g5)
    assert a["sheet_chiral_anticommutator_residuals"][0]==0.
    assert a["inter_sheet_chiral_commutator_residual"]==0.
    assert a["full_chiral_anticommutator_residual"]>0.
    assert np.allclose(h,h.conj().T)

def test_nonhermitian_input_rejected():
    d=np.array([[0.,1.],[0.,0.]],complex)
    with pytest.raises(ValueError,match="Hermitian"):
        assemble_two_sheet_dirac(d)
