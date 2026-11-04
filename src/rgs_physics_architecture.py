"""Central compatibility registry for the RGS–PSU–PRISM–QLG programme.

This module does not claim that the listed sectors are already unified.
It checks shared algebraic conventions and records remaining coupling gates,
so independently developed mechanics cannot be silently identified.
"""
from __future__ import annotations
from typing import Any
import numpy as np

METRIC_SIGNATURE = (-1, 1, 1, 1)
ARCHITECTURE_SCHEMA = "RGS_PSU_PHYSICS_ARCHITECTURE_AUDIT_V1"

SECTOR_REGISTRY: tuple[dict[str, str], ...] = (
    {"sector": "RGS_discrete_geometry", "status": "NUMERICALLY_VALIDATED_COMPONENTS",
     "role": "coordinate seed, tangent spaces, coframes and geometric constraints",
     "missing_bridge": "versioned refining family and covariant continuum limit"},
    {"sector": "Q6_and_PRISM", "status": "FINITE_COMBINATORIAL_CONSTRUCTIONS",
     "role": "incidence complexes, product graphs and discrete spectra",
     "missing_bridge": "common geometric realization and refinement maps"},
    {"sector": "quaternions_2I_Clifford", "status": "FINITE_REPRESENTATION_AUDITS",
     "role": "rotation lifts, selected spin representations and Clifford algebra",
     "missing_bridge": "full Spin(3,1) connection, boosts and global bundle data"},
    {"sector": "Dirac_Einstein_Cartan", "status": "ACTION_DERIVED_LOCAL_OPERATORS",
     "role": "Dirac kinetic Euler operator and algebraic axial contact interaction",
     "missing_bridge": "construct D_a psi from a solved tetrad and spin connection"},
    {"sector": "TEGR", "status": "SEPARATE_GEOMETRIC_FORMULATION",
     "role": "teleparallel connection and torsion-based gravity",
     "missing_bridge": "derive inertial spin connection and boundary identity in one convention"},
    {"sector": "U1_Maxwell_Qball", "status": "CANDIDATE_ACTION_AND_RADIAL_EQUATIONS",
     "role": "two-complex-scalar gauged Q-ball with diagonal U(1)",
     "missing_bridge": "couple its stress/current consistently to the dynamical tetrad and other matter"},
    {"sector": "SU2_lattice_gauge", "status": "FINITE_GAUGE_COVARIANCE_AUDITS",
     "role": "oriented links, path transport and Wilson holonomies",
     "missing_bridge": "derive links from the same continuum connection as the spinor derivative"},
    {"sector": "SU3_QCD", "status": "REPRESENTATION_AND_FLAVOUR_COMPONENTS_ONLY",
     "role": "selected SU(3) representation algebra",
     "missing_bridge": "SU(3) link dynamics, fermion representations, chiral limit and continuum matching"},
    {"sector": "QLG_Dirac_Kahler", "status": "FINITE_OPERATOR_CONSTRUCTIONS",
     "role": "spatial graph Dirac, Wilson/overlap and cochain Dirac–Kähler operators",
     "missing_bridge": "curved refinement, 4D dynamics, compactness and Mosco/spectral convergence"},
    {"sector": "MHD_fluids_superfluids_elasticity", "status": "EFFECTIVE_LIMITS_NOT_YET_UNIFIED",
     "role": "collective matter and continuum mechanics",
     "missing_bridge": "derive constitutive laws and stress exchange from a common matter action"},
    {"sector": "QAT_thermodynamics", "status": "PARTIAL_FINITE_STATE_AUDITS",
     "role": "two-sheet state models, entropy and spectral diagnostics",
     "missing_bridge": "derive physical time, partition function and gravitational backreaction"},
)

COUPLING_GATES: tuple[dict[str, str], ...] = (
    {"gate": "geometry_to_spin_connection",
     "requirement": "Construct omega_mu from the chosen tetrad and connection field equations; do not infer it from a polar-frame SU(2) link alone.",
     "status": "OPEN"},
    {"gate": "spinor_covariant_derivative",
     "requirement": "Implement D_mu psi = partial_mu psi + Omega_mu psi + internal gauge action with explicit index and representation conventions.",
     "status": "OPEN"},
    {"gate": "EC_to_TEGR_separation",
     "requirement": "Keep algebraic Einstein–Cartan torsion distinct from TEGR torsion and any propagating-torsion model.",
     "status": "CONVENTION_GUARD"},
    {"gate": "U1_to_nonabelian_gauge",
     "requirement": "Specify independent U(1), SU(2)_L and SU(3)_C connections, generators, representations and couplings.",
     "status": "OPEN"},
    {"gate": "matter_to_hilbert_tensor",
     "requirement": "Vary the complete matter action before imposing field equations or a fluid/ensemble closure.",
     "status": "OPEN"},
    {"gate": "Noether_and_Bianchi",
     "requirement": "Verify gauge identities, diffeomorphism/local-Lorentz identities and total covariant stress conservation.",
     "status": "OPEN"},
    {"gate": "discrete_to_continuum",
     "requirement": "Use a refining family, common embeddings, uniform geometric bounds and separate convergence proofs for scalar, Dirac and gauge operators.",
     "status": "OPEN"},
    {"gate": "quantum_consistency",
     "requirement": "Check fermion doubling/chirality, gauge and gravitational anomalies, positivity and quantum measure.",
     "status": "OPEN"},
)


def _max_matrix_residual(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(left) - np.asarray(right), ord="fro"))


def audit_core_conventions(*, tolerance: float = 2e-12) -> dict[str, Any]:
    """Cross-check the shared metric, gamma matrices and Q6 chain complex."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")

    from src.clifford_qlg import ETA, gamma_matrices
    from src.discrete_dirac import chirality_matrix
    from src.q6_complex import build_q6_complex
    from src.rgs_ec_algebraic_torsion import (
        gamma5_matrix, minkowski_metric, dirac_gamma_matrices_minus_plus,
    )

    eta_clifford = np.asarray(ETA, dtype=float)
    eta_ec = minkowski_metric()
    gamma_clifford = tuple(np.asarray(g, dtype=complex) for g in gamma_matrices())
    gamma_ec = tuple(np.asarray(g, dtype=complex) for g in dirac_gamma_matrices_minus_plus())
    if len(gamma_clifford) != 4 or len(gamma_ec) != 4:
        raise ArithmeticError("Both Dirac modules must expose four gamma matrices")

    metric_residual = _max_matrix_residual(eta_clifford, eta_ec)
    gamma_residual = max(_max_matrix_residual(a, b)
                         for a, b in zip(gamma_clifford, gamma_ec))
    gamma5_residual = _max_matrix_residual(gamma5_matrix(), chirality_matrix())

    clifford_residual = 0.0
    identity = np.eye(4, dtype=complex)
    for mu in range(4):
        for nu in range(4):
            anti = gamma_clifford[mu] @ gamma_clifford[nu]
            anti += gamma_clifford[nu] @ gamma_clifford[mu]
            clifford_residual = max(
                clifford_residual,
                float(np.linalg.norm(anti - 2.0 * eta_clifford[mu, nu] * identity,
                                     ord="fro")),
            )

    q6 = build_q6_complex()
    chain = q6.boundary_1.astype(np.int64) @ q6.boundary_2.astype(np.int64)
    chain_residual = int(np.max(np.abs(chain))) if chain.size else 0
    checks = {
        "metric_signature_match": metric_residual <= tolerance,
        "gamma_representation_match": gamma_residual <= tolerance,
        "gamma5_representation_match": gamma5_residual <= tolerance,
        "clifford_relation": clifford_residual <= tolerance,
        "q6_chain_complex": chain_residual == 0,
    }
    return {
        "metric_signature": list(METRIC_SIGNATURE),
        "tolerance": float(tolerance),
        "metric_match_residual": metric_residual,
        "gamma_match_max_residual": gamma_residual,
        "gamma5_match_residual": gamma5_residual,
        "clifford_max_residual": clifford_residual,
        "q6": {
            "vertices": len(q6.vertices), "edges": len(q6.edges),
            "faces": len(q6.faces), "B1_shape": list(q6.boundary_1.shape),
            "B2_shape": list(q6.boundary_2.shape),
            "max_abs_B1B2": chain_residual,
        },
        "checks": checks,
        "passed": bool(all(checks.values())),
    }


def physics_architecture_manifest() -> dict[str, Any]:
    """Return the sector and coupling registry without claiming unification."""
    return {
        "schema": ARCHITECTURE_SCHEMA,
        "metric_signature": list(METRIC_SIGNATURE),
        "sectors": [dict(item) for item in SECTOR_REGISTRY],
        "coupling_gates": [dict(item) for item in COUPLING_GATES],
        "interpretation": (
            "A dependency and evidence registry. It is not a unified action, "
            "a solved field configuration, or a continuum/phenomenological proof."
        ),
    }


def audit_physics_architecture(*, tolerance: float = 2e-12) -> dict[str, Any]:
    """Combine convention checks with explicit open coupling gates."""
    core = audit_core_conventions(tolerance=tolerance)
    manifest = physics_architecture_manifest()
    open_gates = [item["gate"] for item in manifest["coupling_gates"]
                  if item["status"] == "OPEN"]
    return {
        "schema": ARCHITECTURE_SCHEMA,
        "core_conventions": core,
        "sector_count": len(manifest["sectors"]),
        "coupling_gate_count": len(manifest["coupling_gates"]),
        "open_coupling_gates": open_gates,
        "all_core_conventions_compatible": core["passed"],
        "unification_status": (
            "CORE_CONVENTIONS_COMPATIBLE_BUT_PHYSICAL_COUPLINGS_OPEN"
            if core["passed"] else "CORE_CONVENTION_CONFLICT"
        ),
        "manifest": manifest,
    }
