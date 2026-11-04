"""Analytic coframe jet on the Ellis-RGS continuum refinement chart.

Coordinates are (t,u,theta,phi), with
 ds^2=-dt^2+du^2+(u^2+r0^2)(dtheta^2+sin(theta)^2 dphi^2).
The chart is a continuum reference geometry associated with the RGS angular
refinement; it is not a reconstruction of derivatives of the 64 independent
site covariance frames.
"""
from __future__ import annotations

import numpy as np

from src.rgs_spin_connection import levi_civita_spin_connection


def ellis_coframe_jet(u: float, theta: float, r0: float = 0.8292):
    """Return e[a,mu] and de[a,nu,mu]=partial_mu e[a,nu] analytically."""
    u = float(u)
    theta = float(theta)
    r0 = float(r0)
    if not np.all(np.isfinite([u, theta, r0])) or r0 <= 0.0:
        raise ValueError("u, theta and r0 must be finite; r0 must be positive")
    if not 0.0 < theta < np.pi:
        raise ValueError("theta must lie strictly between 0 and pi (regular chart)")
    radius = float(np.hypot(u, r0))
    e = np.diag([1.0, 1.0, radius, radius * np.sin(theta)])
    de = np.zeros((4, 4, 4), dtype=float)
    # Coordinate order: (t,u,theta,phi); derivative index is last.
    de[2, 2, 1] = u / radius
    de[3, 3, 1] = (u / radius) * np.sin(theta)
    de[3, 3, 2] = radius * np.cos(theta)
    return e, de


def ellis_rgs_spin_connection(u: float, theta: float, r0: float = 0.8292):
    """Evaluate the torsion-free LC spin connection on the Ellis chart."""
    e, de = ellis_coframe_jet(u, theta, r0)
    report = levi_civita_spin_connection(e, de)
    report["coframe"] = e
    report["coframe_derivatives"] = de
    report["geometry"] = "ELLIS_RGS_ANALYTIC_CONTINUUM_CHART"
    report["coordinate_order"] = ["t", "u", "theta", "phi"]
    report["scope"] = (
        "analytic Ellis continuum chart sampled over the RGS angular refinement; "
        "not yet the derivative of the 64 site-dependent RGS covariance coframes"
    )
    return report
