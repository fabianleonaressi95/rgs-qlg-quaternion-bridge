"""Two-sheet PSU extension of the repository's finite Dirac operators.

Couples two already-assembled Hermitian one-sheet operators through a unitary
map on the complete site-spinor space. It does not assemble a replacement
spatial Dirac discretisation. Spatial SU(2) spin links, electromagnetic U(1),
colour SU(3), and sheet transport remain distinct.
"""
from __future__ import annotations
import numpy as np

def _square_matrix(value, name):
    a=np.asarray(value,dtype=complex)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be a finite square matrix")
    return a

def assemble_two_sheet_dirac(dirac_plus, dirac_minus=None, transporter=None,
                              coupling=0.0, mass_plus=0.0, mass_minus=0.0,
                              gamma5=None, tolerance=1e-9):
    """Build [[D+ + m+I, gW], [gW†, D- + m-I]].

    W must be unitary on the full one-sheet site-spinor space. The routine
    checks, rather than assumes, Hermiticity and optional chirality relations.
    """
    dp=_square_matrix(dirac_plus,"dirac_plus")
    dm=dp.copy() if dirac_minus is None else _square_matrix(dirac_minus,"dirac_minus")
    if dm.shape!=dp.shape: raise ValueError("sheet Dirac blocks must have equal shape")
    n=len(dp); g=float(coupling); mp=float(mass_plus); mm=float(mass_minus); tol=float(tolerance)
    if not all(np.isfinite(v) for v in (g,mp,mm,tol)) or tol<=0:
        raise ValueError("couplings must be finite and tolerance positive")
    hp=float(np.linalg.norm(dp-dp.conj().T,"fro"))
    hm=float(np.linalg.norm(dm-dm.conj().T,"fro"))
    if max(hp,hm)>tol: raise ValueError("input Dirac blocks must be Hermitian")
    w=np.eye(n,dtype=complex) if transporter is None else _square_matrix(transporter,"transporter")
    if w.shape!=(n,n): raise ValueError("transporter must act on the complete one-sheet space")
    wu=float(np.linalg.norm(w.conj().T@w-np.eye(n),"fro"))
    if wu>tol: raise ValueError("inter-sheet transporter must be unitary")
    eye=np.eye(n,dtype=complex)
    h=np.block([[dp+mp*eye,g*w],[g*w.conj().T,dm+mm*eye]])
    report={"operator_kind":"FINITE_TWO_SHEET_HERMITIAN_DIRAC_BLOCK",
            "one_sheet_dimension":int(n),"total_dimension":int(2*n),
            "coupling":g,"mass_plus":mp,"mass_minus":mm,
            "input_plus_hermiticity_residual":hp,"input_minus_hermiticity_residual":hm,
            "transporter_unitarity_residual":wu,
            "hermiticity_residual":float(np.linalg.norm(h-h.conj().T,"fro")),
            "physical_status":"FINITE_OPERATOR_CANDIDATE"}
    if gamma5 is not None:
        ch=_square_matrix(gamma5,"gamma5")
        if ch.shape!=(n,n): raise ValueError("gamma5 has incompatible shape")
        z=np.zeros_like(ch); ch2=np.block([[ch,z],[z,ch]])
        report.update({
            "gamma5_involution_residual":float(np.linalg.norm(ch@ch-eye,"fro")),
            "gamma5_hermiticity_residual":float(np.linalg.norm(ch-ch.conj().T,"fro")),
            "sheet_chiral_anticommutator_residuals":[
                float(np.linalg.norm(dp@ch+ch@dp,"fro")),
                float(np.linalg.norm(dm@ch+ch@dm,"fro"))],
            "inter_sheet_chiral_commutator_residual":float(np.linalg.norm(ch@w-w@ch,"fro")),
            "full_chiral_anticommutator_residual":float(np.linalg.norm(h@ch2+ch2@h,"fro"))})
    return h,report

def two_sheet_spectrum(operator,tolerance=1e-9):
    h=_square_matrix(operator,"operator")
    herm=float(np.linalg.norm(h-h.conj().T,"fro"))
    if herm>tolerance: raise ValueError(f"operator is not Hermitian: {herm}")
    v=np.linalg.eigvalsh(h)
    return {"eigenvalues":v,"dimension":int(len(v)),"minimum":float(v[0]),
            "maximum":float(v[-1]),"near_zero_count":int(np.count_nonzero(np.abs(v)<=tolerance)),
            "hermiticity_residual":herm}

def symmetric_sheet_spectrum(single_sheet_operator,coupling):
    """Exact sorted union spec(D+gI) U spec(D-gI) for identical sheets, W=I."""
    d=_square_matrix(single_sheet_operator,"single_sheet_operator")
    if np.linalg.norm(d-d.conj().T,"fro")>1e-9: raise ValueError("D must be Hermitian")
    g=float(coupling)
    if not np.isfinite(g): raise ValueError("coupling must be finite")
    v=np.linalg.eigvalsh(d)
    return np.sort(np.concatenate((v+g,v-g)))
