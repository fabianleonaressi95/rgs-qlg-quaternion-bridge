"""Vierbein construction and metric reconstruction for RGS-PSU.

This module implements:
- Vierbein E_i from covariance matrices C_i
- Metric reconstruction g_μν = η_AB e^A_μ e^B_ν
- Killing vector analysis
- Symmetry verification
"""

from .construction import build_vierbein, reconstruct_metric
from .killing import KillingOperator, find_killing_vectors, verify_killing_conservation
from .audit import audit_vierbein, audit_metric, audit_killing_symmetries

__all__ = [
    "build_vierbein",
    "reconstruct_metric",
    "KillingOperator",
    "find_killing_vectors",
    "verify_killing_conservation",
    "audit_vierbein",
    "audit_metric",
    "audit_killing_symmetries",
]
