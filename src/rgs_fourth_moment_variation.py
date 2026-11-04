"""Variational coframe coupling for the RGS fourth-moment anisotropy.

Algebraic 3+1 constitutive sector. Its coframe Euler derivative is exact;
it deliberately does not identify this term with teleparallel torsion.
"""
from __future__ import annotations
import numpy as np


def anisotropic_potential(residual, q, coupling=1.0):
    r, q = np.asarray(residual, float), np.asarray(q, float)
    if r.shape != (3,3,3,3) or q.shape != (3,3):
        raise ValueError("expected residual (3,3,3,3) and q (3,3)")
    if not np.allclose(q, q.T, atol=1e-12, rtol=0.0) or abs(np.trace(q)) > 1e-10:
        raise ValueError("q must be symmetric and trace-free")
    if not np.isfinite(coupling) or not np.all(np.isfinite(r)) or not np.all(np.isfinite(q)):
        raise ValueError("inputs must be finite")
    return 0.5 * float(coupling) * float(np.einsum("abcd,ab,cd", r, q, q))


def coframe_action_density(coframe, lapse, residual, q, coupling=1.0):
    e = np.asarray(coframe, float)
    if e.shape != (3,3) or not np.all(np.isfinite(e)):
        raise ValueError("coframe must be finite 3x3")
    det_e = float(np.linalg.det(e))
    if det_e <= 0:
        raise ValueError("coframe must be orientation-preserving and nonsingular")
    if not np.isfinite(lapse) or lapse <= 0:
        raise ValueError("lapse must be finite and positive")
    return float(lapse) * det_e * anisotropic_potential(residual, q, coupling)


def coframe_euler_derivative(coframe, lapse, residual, q, coupling=1.0):
    """dL/d e[a,i] = N det(e) V (e^{-1})[i,a]; returned array is indexed [a,i], hence inv(e).T."""
    e = np.asarray(coframe, float)
    return coframe_action_density(e, lapse, residual, q, coupling) * np.linalg.inv(e).T


def lapse_euler_derivative(coframe, residual, q, coupling=1.0):
    e = np.asarray(coframe, float)
    if e.shape != (3,3) or np.linalg.det(e) <= 0:
        raise ValueError("coframe must be orientation-preserving and nonsingular")
    return float(np.linalg.det(e)) * anisotropic_potential(residual, q, coupling)


def q_euler_derivative(residual, q, coupling=1.0):
    return float(coupling) * np.einsum("abcd,cd->ab", np.asarray(residual,float), np.asarray(q,float))


def rotate_rank4(tensor, rotation):
    t, u = np.asarray(tensor,float), np.asarray(rotation,float)
    if t.shape != (3,3,3,3) or u.shape != (3,3):
        raise ValueError("expected tensor (3,3,3,3), rotation (3,3)")
    if not np.allclose(u.T@u, np.eye(3), atol=1e-10, rtol=0.0):
        raise ValueError("rotation must be orthogonal")
    return np.einsum("ia,jb,kc,ld,abcd->ijkl",u,u,u,u,t)


def rotate_rank2(tensor, rotation):
    t, u = np.asarray(tensor,float), np.asarray(rotation,float)
    if t.shape != (3,3) or u.shape != (3,3):
        raise ValueError("expected tensor and rotation (3,3)")
    return u @ t @ u.T


def audit_variational_bridge(residual, q, coframe=None, lapse=1.2, coupling=.7, step=1e-6):
    e = np.asarray(coframe if coframe is not None else
                   [[1.1,.1,0.],[0.,.9,.05],[.02,0.,1.2]], float)
    grad = coframe_euler_derivative(e,lapse,residual,q,coupling)
    fd = np.zeros((3,3))
    for a in range(3):
        for i in range(3):
            d=np.zeros((3,3)); d[a,i]=step
            fd[a,i]=(coframe_action_density(e+d,lapse,residual,q,coupling)
                     -coframe_action_density(e-d,lapse,residual,q,coupling))/(2*step)
    fd_n=(coframe_action_density(e,lapse+step,residual,q,coupling)
          -coframe_action_density(e,lapse-step,residual,q,coupling))/(2*step)
    from src.rgs_fourth_moment_bridge import symmetric_traceless_basis
    gq=q_euler_derivative(residual,q,coupling)
    q_errors=[]
    for b in symmetric_traceless_basis():
        f=(anisotropic_potential(residual,q+step*b,coupling)
           -anisotropic_potential(residual,q-step*b,coupling))/(2*step)
        q_errors.append(abs(f-float(np.einsum("ab,ab",gq,b))))
    return {
        "action_density": coframe_action_density(e,lapse,residual,q,coupling),
        "coframe_gradient_fd_max_abs": float(np.max(np.abs(grad-fd))),
        "lapse_gradient_fd_abs": abs(lapse_euler_derivative(e,residual,q,coupling)-fd_n),
        "q_gradient_sym0_fd_max_abs": float(max(q_errors)),
        "coframe_gradient": grad.tolist(), "q_gradient": gq.tolist(),
        "status": {
            "algebraic_coframe_variation":"NUMERICALLY_VERIFIED",
            "q_constitutive_response":"NUMERICALLY_VERIFIED",
            "local_frame_covariance":"COVARIANT_IF_R_AND_Q_TRANSFORM_AS_TENSORS",
            "teleparallel_torsion_generation":"NOT_GENERATED_BY_THIS_ALGEBRAIC_TERM",
            "full_coframe_field_equation":"OPEN",
            "flat_spin_connection_variation":"OPEN",
            "continuum_limit":"OPEN"
        }
    }
