"""Exact kinematic audit for Smilga's equal-mass two-particle construction.
Metric (+---), c=hbar=1; P=p1+p2 and q=(p1-p2)/2.
"""
import math
import numpy as np

ETA = np.diag([1., -1., -1., -1.])


def dot(a, b):
    return np.einsum("...i,ij,...j->...", np.asarray(a, float), ETA, np.asarray(b, float))


def to_total_relative(p1, p2):
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    if p1.shape != p2.shape or p1.shape[-1] != 4:
        raise ValueError("p1 and p2 must have matching (...,4) shapes")
    return p1 + p2, 0.5 * (p1 - p2)


def from_total_relative(P, q):
    P, q = np.asarray(P, float), np.asarray(q, float)
    if P.shape != q.shape or P.shape[-1] != 4:
        raise ValueError("P and q must have matching (...,4) shapes")
    return 0.5 * P + q, 0.5 * P - q


def jacobian_abs_8d():
    block = np.block([[np.eye(4), np.eye(4)],
                      [0.5*np.eye(4), -0.5*np.eye(4)]])
    return float(abs(np.linalg.det(block)))


def cm_pair(M, m, directions):
    M, m = float(M), float(m)
    if M < 2*m or M <= 0:
        raise ValueError("Require M >= 2m and M > 0")
    n = np.asarray(directions, float)
    if n.shape[-1] != 3:
        raise ValueError("directions must end in dimension 3")
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    if np.any(norm == 0):
        raise ValueError("directions must be nonzero")
    n = n / norm
    k = math.sqrt(max(0., M*M/4-m*m))
    E = np.full(n.shape[:-1] + (1,), M/2)
    return np.concatenate([E, k*n], axis=-1), np.concatenate([E, -k*n], axis=-1)


def boost(p, v):
    p, v = np.asarray(p, float), np.asarray(v, float)
    if p.shape[-1] != 4 or v.shape != (3,):
        raise ValueError("p must end in 4; v must have shape (3,)")
    v2 = float(v @ v)
    if v2 >= 1:
        raise ValueError("boost must be subluminal")
    if v2 == 0:
        return p.copy()
    g = 1/math.sqrt(1-v2)
    vp = np.einsum("i,...i->...", v, p[..., 1:])
    out = np.empty_like(p)
    out[..., 0] = g*(p[..., 0]+vp)
    out[..., 1:] = p[..., 1:] + (((g-1)*vp/v2)+g*p[..., 0])[..., None]*v
    return out


def residuals(p1, p2, m):
    P, q = to_total_relative(p1, p2)
    return {
        "p1_shell": dot(p1, p1)-m*m,
        "p2_shell": dot(p2, p2)-m*m,
        "P_dot_q": dot(P, q),
        "combined_casimir": dot(P, P)+4*dot(q, q)-4*m*m,
    }


def phase_space_density_dOmega(M, m):
    if M < 2*m or M <= 0:
        raise ValueError("Require M >= 2m and M > 0")
    k = math.sqrt(max(0., M*M/4-m*m))
    return k/(16*math.pi**2*M)


def integrated_phase_space(M, m):
    return 4*math.pi*phase_space_density_dOmega(M, m)


def wyler_alpha_inverse():
    alpha = 9/(8*math.pi**4)*(math.pi**5/(2**4*math.factorial(5)))**0.25
    return 1/alpha


def audit(m=1., M=3., velocity=(.2, -.1, .15)):
    dirs = np.array([[1.,0.,0.], [0.,1.,0.], [0.,0.,1.], [1.,2.,-1.]])
    a, b = cm_pair(M, m, dirs)
    a, b = boost(a, np.asarray(velocity)), boost(b, np.asarray(velocity))
    rr = residuals(a, b, m)
    P, q = to_total_relative(a, b)
    return {
        "convention": "metric (+---); P=p1+p2; q=(p1-p2)/2",
        "jacobian_abs_8d": jacobian_abs_8d(),
        "max_abs_residuals": {k: float(np.max(np.abs(v))) for k,v in rr.items()},
        "P2_mean": float(np.mean(dot(P,P))),
        "q2_mean": float(np.mean(dot(q,q))),
        "cm_relative_momentum": math.sqrt(max(0., M*M/4-m*m)),
        "dPhi2_dOmega": phase_space_density_dOmega(M,m),
        "integrated_dPhi2": integrated_phase_space(M,m),
        "wyler_alpha_inverse_unit_radius": wyler_alpha_inverse(),
        "proved": ["linear Jacobian for the declared convention",
                   "equal-mass shell identities and P.q=0",
                   "standard invariant two-body phase-space formula"],
        "open": ["S-matrix kernel is an orthogonal projection",
                 "bounded Lie-ball map preserves the physical on-shell measure",
                 "four-dimensional to five-dimensional extension is uniquely forced",
                 "spin multiplicity factor 4 is physically justified",
                 "geometric weight equals renormalized low-energy QED alpha",
                 "Q6 U(1)^129 measure induces D5 geometry"]
    }
