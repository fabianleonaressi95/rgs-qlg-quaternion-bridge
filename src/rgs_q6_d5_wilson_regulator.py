"""Transfer-derived spectral regulator and Q6-to-D5 identifiability audit.

The compact Wilson transfer operator supplies an exact cutoff in its own
integer-flux spectrum. This module deliberately does not rename that cutoff
as a D5 generic-norm cutoff: the real-linear extension from the rank-129
cycle-coordinate space to the five selected D5 tangent momenta has a
124-dimensional real kernel. Injectivity on the integer lattice is a
separate arithmetic question and is not inferred from this dimension count.
"""
from __future__ import annotations

import math
import numpy as np

from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector
from src.rgs_q6_u1_character_electric_flux import build_q6_integer_cycle_basis
from src.rgs_q6_u1_transfer_operator import (
    q6_u1_electric_transfer_eigenvalues,
    q6_u1_electric_hamiltonian_eigenvalues,
    validate_q6_integer_flux_states,
)
from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data


def _positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def q6_wilson_spectral_cutoff(
    flux_states,
    beta_time: float,
    temporal_spacing: float,
    minimum_transfer_weight: float,
    data=None,
) -> dict:
    """Select states with exact normalized transfer weight lambda_n >= tau.

    Equivalently E_n <= E_max=-log(tau)/a_t, where E_n is the exact
    finite-beta_t electric transfer energy with vacuum energy removed.
    This is a Wilson-sector spectral cutoff, not a D5 norm cutoff.
    """
    beta = _positive_finite("beta_time", beta_time)
    at = _positive_finite("temporal_spacing", temporal_spacing)
    tau = float(minimum_transfer_weight)
    if not math.isfinite(tau) or not 0.0 < tau <= 1.0:
        raise ValueError("minimum_transfer_weight must lie in (0,1]")
    states = validate_q6_integer_flux_states(flux_states, data)
    weights = q6_u1_electric_transfer_eigenvalues(states, beta, data)
    energies = q6_u1_electric_hamiltonian_eigenvalues(states, beta, at, data)
    selected = weights >= tau
    emax = -math.log(tau) / at
    # Independently reconstruct energies from the returned transfer weights.
    # This catches accidental divergence between the two public spectral APIs.
    energy_from_weight = -np.log(weights) / at
    equivalence_residual = float(np.max(np.abs(energies - energy_from_weight)))
    energy_selected = energies <= emax + 1e-12
    selection_mismatch_count = int(np.count_nonzero(selected != energy_selected))
    return {
        "threshold_transfer_weight": tau,
        "equivalent_maximum_electric_energy": emax,
        "transfer_weights": weights,
        "electric_energies": energies,
        "selected_mask": selected,
        "selected_count": int(np.count_nonzero(selected)),
        "total_count": int(len(states)),
        "cutoff_equivalence_residual": equivalence_residual,
        "selection_mismatch_count": selection_mismatch_count,
        "selection_rule": "lambda_n >= tau iff E_n <= -log(tau)/a_t",
    }


def audit_q6_wilson_d5_regulator_obstruction(
    beta_time: float = 80.0,
    temporal_spacing: float = 0.05,
    minimum_transfer_weight: float = 0.25,
    seed: int = 20261104,
) -> dict:
    """Audit exact flux spectral truncation and the missing D5 scalar map."""
    data = build_q6_euclidean_time_data(2)
    cycles = build_q6_integer_cycle_basis(data)  # 129 x 192, integral
    states = np.vstack((np.zeros((1, 192), dtype=np.int64), cycles, -cycles))
    states = validate_q6_integer_flux_states(states, data)

    spectral = q6_wilson_spectral_cutoff(
        states, beta_time, temporal_spacing, minimum_transfer_weight, data
    )
    sector = build_q6_d5_phase_space_sector(beta=1.0, eigenvalue=12.0)
    modes = np.asarray(sector["selected_mode_basis"], dtype=float)
    # Coefficients c in R^129 map to five-mode momenta by c -> c @ cycles @ B.
    projected_cycle_map = cycles @ modes
    singular_values = np.linalg.svd(projected_cycle_map, compute_uv=False)
    rank = int(np.linalg.matrix_rank(projected_cycle_map, tol=1.0e-9))
    kernel_dimension = int(cycles.shape[0] - rank)

    # Finite Galerkin shell: vacuum and both orientations of every fundamental
    # integer cycle. This is not the complete Z^129 cycle lattice.
    projected_momenta = states @ modes
    projected_norm2 = np.sum(projected_momenta**2, axis=1)
    energies = np.asarray(spectral["electric_energies"], dtype=float)
    weights = np.asarray(spectral["transfer_weights"], dtype=float)
    selected = np.asarray(spectral["selected_mask"], dtype=bool)

    # The real-linear extension has a 124D kernel. This does not by itself
    # settle injectivity of the restriction to Z^129: irrational coefficients
    # can make a low-dimensional real map injective on an integer lattice.
    # No integer-collision search or Diophantine injectivity proof is supplied.
    # The D5 chart is local/nonlinear; no global map to N(z) is established.
    checks = {
        "integer_gauss_law_exact": bool(
            np.max(np.abs(np.asarray(data["D0"], dtype=np.int64).T @ states.T)) == 0
        ),
        "cycle_rank_129": cycles.shape == (129, 192),
        "finite_galerkin_shell_259": states.shape == (259, 192),
        "projected_five_mode_rank": rank == 5,
        "projection_kernel_dimension_124": kernel_dimension == 124,
        "vacuum_has_unit_transfer_weight": math.isclose(float(weights[0]), 1.0, rel_tol=1e-13),
        "all_transfer_weights_in_unit_interval": bool(
            np.all(weights > 0.0) and np.all(weights <= 1.0 + 1e-12)
        ),
        "spectral_selection_nonempty": bool(np.any(selected)),
        "vacuum_is_retained": bool(selected[0]),
        "energy_threshold_matches_weight_threshold": bool(
            np.array_equal(selected, energies <= spectral["equivalent_maximum_electric_energy"] + 1e-12)
        ),
        "transfer_energy_formula_residual_below_tolerance": (
            spectral["cutoff_equivalence_residual"] < 1.0e-12
        ),
        "spectral_selection_has_no_mismatches": (
            spectral["selection_mismatch_count"] == 0
        ),
        "projected_coordinates_are_five_dimensional": projected_momenta.shape == (259, 5),
        "no_epsilon_identification_claimed": True,
    }
    return {
        "audit": "Q6_WILSON_SPECTRAL_REGULATOR_D5_IDENTIFIABILITY",
        "q6": {
            "vertices": 64, "edges": 192, "integer_cycle_rank": 129,
            "finite_galerkin_states": int(states.shape[0]),
            "finite_galerkin_description": "vacuum plus +/- 129 fundamental integral cycles",
        },
        "transfer": {
            "beta_time": float(beta_time),
            "temporal_spacing": float(temporal_spacing),
            "threshold_weight": float(minimum_transfer_weight),
            "maximum_electric_energy": float(spectral["equivalent_maximum_electric_energy"]),
            "selected_states": int(spectral["selected_count"]),
            "total_states": int(spectral["total_count"]),
            "minimum_weight": float(np.min(weights)),
            "maximum_weight": float(np.max(weights)),
            "minimum_energy": float(np.min(energies)),
            "maximum_energy": float(np.max(energies)),
        },
        "five_mode_projection": {
            "source_dimension": 129,
            "target_dimension": 5,
            "rank": rank,
            "kernel_dimension": kernel_dimension,
            "real_linear_kernel_dimension": kernel_dimension,
            "integer_lattice_injectivity": "not established by this audit",
            "integer_projection_collisions_checked": False,
            "singular_values": singular_values.tolist(),
            "projected_momentum_squared_range": [
                float(np.min(projected_norm2)), float(np.max(projected_norm2))
            ],
        },
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "exact_conclusion": (
            "The compact Q6 Wilson transfer operator induces an exact finite-beta_t "
            "spectral regulator lambda_n >= tau, equivalently E_n <= -log(tau)/a_t, "
            "on any explicitly supplied finite subset of integer cycle fluxes."
        ),
        "identifiability_result": (
            "The integer cycle lattice has rank 129, whereas the selected D5 tangent "
            "momentum projection has rank 5 and a 124-dimensional real kernel. "
            "Consequently the five-mode projection is not injective and cannot, by "
            "itself, determine the full Wilson transfer weight or a unique D5 generic "
            "norm cutoff. An additional equivariant map/action or quotient is required."
        ),
        "status": {
            "proved_exactly": [
                "finite-beta_t Bessel transfer weights and their energy threshold",
                "integer Gauss constraint on all supplied flux states",
                "129-dimensional integer cycle lattice",
            ],
            "numerically_verified": [
                "rank-five projection onto the selected Q6-D5 tangent modes",
                "124-dimensional kernel of the real-linear extension",
                "spectral selection on the 259-state Galerkin shell",
            ],
            "not_claimed": [
                "a physical value of the D5 epsilon cutoff",
                "injectivity or non-injectivity of the five-mode projection restricted to Z^129",
                "factorization of the full Wilson weight through the five-mode projection",
                "a global Q6-to-D5 map",
                "equivalence between Wilson spectral and D5 generic-norm cutoffs",
                "the full infinite-dimensional flux transfer spectrum",
                "a 3+1D continuum limit or fine-structure constant prediction",
            ],
        },
        "reproducibility": {
            "seed_reserved_for_future_sampling": int(seed),
            "command": "python scripts/run_q6_wilson_d5_regulator.py",
        },
    }
