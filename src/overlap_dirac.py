"""Finite-dimensional overlap construction from the repository's Wilson pair.

The repository stores a Hermitian Hamiltonian-like naive operator H0 with
{Gamma5,H0}=0. The Euclidean gamma5-Hermitian Wilson kernel used here is
D_E = i H0 + W, and H_W = Gamma5 (a D_E - rho I).
The overlap operator is D_ov = (rho/a) [I + Gamma5 sign(H_W)].

This is an explicit algebraic conversion, not a claim that the Q6 operator
already has a continuum Euclidean Dirac interpretation.
"""
from __future__ import annotations
import numpy as np


def assemble_overlap_from_wilson(
    naive_hermitian, wilson_laplacian, gamma5, lattice_spacing=1.0,
    rho=1.0, zero_tolerance=1e-11
):
    h0=np.asarray(naive_hermitian,dtype=complex)
    w=np.asarray(wilson_laplacian,dtype=complex)
    g=np.asarray(gamma5,dtype=complex)
    a=float(lattice_spacing); rho=float(rho)
    if h0.ndim!=2 or h0.shape[0]!=h0.shape[1] or w.shape!=h0.shape or g.shape!=h0.shape:
        raise ValueError("naive, Wilson and Gamma5 matrices must be equal-size square matrices")
    if not np.isfinite(a) or a<=0 or not np.isfinite(rho) or rho<=0:
        raise ValueError("lattice_spacing and rho must be positive and finite")
    if np.linalg.norm(h0-h0.conj().T)>1e-9 or np.linalg.norm(w-w.conj().T)>1e-9:
        raise ValueError("naive and Wilson operators must be Hermitian")
    if np.linalg.norm(g-g.conj().T)>1e-9 or np.linalg.norm(g@g-np.eye(len(g)))>1e-9:
        raise ValueError("Gamma5 must be a Hermitian involution")
    naive_chiral=float(np.linalg.norm(h0@g+g@h0,ord="fro"))
    wilson_comm=float(np.linalg.norm(w@g-g@w,ord="fro"))
    if naive_chiral>1e-8 or wilson_comm>1e-8:
        raise ValueError("overlap kernel requires {Gamma5,H0}=0 and [Gamma5,W]=0")
    d_e=1j*h0+w
    kernel=g@(a*d_e-rho*np.eye(len(g)))
    kernel_herm=float(np.linalg.norm(kernel-kernel.conj().T,ord="fro"))
    if kernel_herm>1e-8:
        raise ValueError(f"gamma5 Wilson kernel is not Hermitian: {kernel_herm}")
    values,vectors=np.linalg.eigh((kernel+kernel.conj().T)/2)
    min_gap=float(np.min(np.abs(values)))
    if min_gap<=zero_tolerance:
        raise ValueError("Wilson kernel has a zero/near-zero eigenvalue; sign is undefined")
    sign=(vectors*np.sign(values))@vectors.conj().T
    identity=np.eye(len(g),dtype=complex)
    overlap=(rho/a)*(identity+g@sign)
    gw=float(np.linalg.norm(g@overlap+overlap@g-(a/rho)*overlap@g@overlap,ord="fro"))
    gamma5_herm=float(np.linalg.norm(overlap.conj().T-g@overlap@g,ord="fro"))
    return overlap,{
        "dimension":len(g),
        "lattice_spacing":a,
        "rho":rho,
        "wilson_kernel_hermiticity_residual":kernel_herm,
        "wilson_kernel_minimum_absolute_eigenvalue":min_gap,
        "ginsparg_wilson_residual":gw,
        "gamma5_hermiticity_residual":gamma5_herm,
        "naive_chiral_residual":naive_chiral,
        "wilson_gamma5_commutator_residual":wilson_comm,
        "status":"FINITE_MATRIX_OVERLAP_GW_OPERATOR",
        "continuum_locality_and_index_theorem":"OPEN"
    }
