"""Coarse Q(omega), E(omega) branch audit for the PSU–Anti-Twin Q-ball.

This is a finite-radius parameter sweep, not a stability proof. The sign of
dQ/domega is reported as a branch diagnostic only.
"""
from __future__ import annotations

import json
import numpy as np

from src.rgs_psu_antitwin_qball_bvp import psu_antitwin_qball_maxwell_bvp


def run_branch_audit() -> dict:
    central_amplitudes = np.asarray((0.62, 0.63, 0.64, 0.65, 0.66, 0.67))
    rows = []
    for amplitude in central_amplitudes:
        report = psu_antitwin_qball_maxwell_bvp(
            gauge_charge=0.05,
            central_amplitude=float(amplitude),
            radius_max=20.0,
            points=160,
            mass_scale=1.0,
            quartic=1.5,
            sextic=1.0,
            imbalance_penalty=0.2,
            exchange_coupling=0.1,
        )
        rows.append({
            "central_amplitude": float(amplitude),
            "omega": float(report["omega"]),
            "noether_charge": float(report["noether_particle_charge"]),
            "energy": float(report["energy"]),
            "energy_per_charge": float(report["energy_per_particle_charge"]),
            "gauss_law_residual": float(report["gauss_law_residual"]),
            "plus_eom_residual_max": float(report["plus_eom_residual_max"]),
            "minus_eom_residual_max": float(report["minus_eom_residual_max"]),
            "maxwell_eom_residual_max": float(report["maxwell_eom_residual_max"]),
        })

    rows.sort(key=lambda item: item["omega"])
    omega = np.asarray([row["omega"] for row in rows])
    charge = np.asarray([row["noether_charge"] for row in rows])
    derivative = np.gradient(charge, omega, edge_order=2)
    for row, slope in zip(rows, derivative):
        row["dQ_domega_finite_difference"] = float(slope)

    return {
        "schema": "RGS_PSU_ANTITWIN_QBALL_BRANCH_AUDIT_V1",
        "parameters": {
            "mass_scale": 1.0, "quartic": 1.5, "sextic": 1.0,
            "imbalance_penalty": 0.2, "exchange_coupling": 0.1,
            "gauge_charge": 0.05, "radius_max": 20.0,
        },
        "rows_sorted_by_omega": rows,
        "interpretation": (
            "Finite-difference branch diagnostic on a finite radial domain. "
            "A sign change in dQ/domega indicates a candidate turning point; "
            "it is not a proof of linear or nonlinear stability."
        ),
        "full_stability_established": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_branch_audit(), indent=2, sort_keys=True))
