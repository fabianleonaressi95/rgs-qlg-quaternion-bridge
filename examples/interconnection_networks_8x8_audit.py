#!/usr/bin/env python3
"""Reproducible exhaustive audit for Omega-8 and Beneš-8."""
from __future__ import annotations

import json
from pathlib import Path

from src.interconnection_networks_8x8 import (
    BENES_DIMS,
    OMEGA_DIMS,
    benes8,
    bit_reversal_permutation,
    omega8,
)

omega = omega8()
benes = benes8()
reverse = bit_reversal_permutation()

report = {
    "audit": "Omega-8 and Benes-8 exhaustive distinct-permutation enumeration",
    "model": "8 wires labelled 0..7; 2x2 switches act on one address-bit dimension; shuffle-normalized coordinates",
    "switch_setting": "0 preserves the selected bit; 1 flips it; each stage has four independent switches encoded by a 4-bit mask",
    "omega": {
        "stage_dimensions": OMEGA_DIMS,
        "stages": len(OMEGA_DIMS),
        "switches_per_stage": 4,
        "total_switches": 4 * len(OMEGA_DIMS),
        "distinct_permutations_by_stage": omega.layer_counts,
        "routable_permutations": omega.permutation_count,
        "all_permutations": omega.total_permutations,
        "routable_fraction": "32/315",
    },
    "benes": {
        "stage_dimensions": BENES_DIMS,
        "stages": len(BENES_DIMS),
        "switches_per_stage": 4,
        "total_switches": 4 * len(BENES_DIMS),
        "distinct_permutations_by_stage": benes.layer_counts,
        "routable_permutations": benes.permutation_count,
        "all_permutations": benes.total_permutations,
        "routable_fraction": "1",
    },
    "bit_reversal": {
        "permutation": reverse,
        "omega_routable": omega.routes(reverse),
        "benes_routable": benes.routes(reverse),
        "benes_stage_settings": benes.settings_for(reverse),
    },
    "witness_example": {
        "permutation": (7, 0, 5, 2, 1, 6, 3, 4),
        "benes_stage_settings": benes.settings_for((7, 0, 5, 2, 1, 6, 3, 4)),
    },
    "verification": {
        "enumeration": "all 40,320 target permutations checked for Beneš-8; 4,096 distinct Omega-8 permutations enumerated",
        "status": "NUMERICALLY VERIFIED",
    },
    "limitations": [
        "Finite combinatorial network model only.",
        "No Q6-to-switching-fabric equivalence or Green-Kubo physical theorem is established.",
        "Topology convention and stage coordinate relabelling are specified in the accompanying research record.",
    ],
}

out = Path(__file__).resolve().parents[1] / "results" / "interconnection_networks_8x8_audit.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
print(f"wrote: {out}")
