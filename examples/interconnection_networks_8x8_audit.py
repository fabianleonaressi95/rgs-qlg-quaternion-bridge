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
    "scope": "8x8 finite switching networks",
    "wire_labels": "three-bit integers 0..7",
    "switch_model": "four independent 2x2 switches per stage",
    "switch_setting": "each switch preserves or flips its stage bit dimension",
    "omega": {
        "stage_dimensions": OMEGA_DIMS,
        "stages": len(OMEGA_DIMS),
        "switches_per_stage": 4,
        "total_switches": 4 * len(OMEGA_DIMS),
        "distinct_permutations_by_stage": omega.layer_counts,
        "routable_permutations": omega.permutation_count,
        "all_permutations": omega.total_permutations,
        "fraction": omega.permutation_count / omega.total_permutations,
    },
    "benes": {
        "stage_dimensions": BENES_DIMS,
        "stages": len(BENES_DIMS),
        "switches_per_stage": 4,
        "total_switches": 4 * len(BENES_DIMS),
        "distinct_permutations_by_stage": benes.layer_counts,
        "routable_permutations": benes.permutation_count,
        "all_permutations": benes.total_permutations,
        "fraction": benes.permutation_count / benes.total_permutations,
    },
    "bit_reversal": {
        "permutation": reverse,
        "omega_routable": omega.routes(reverse),
        "benes_routable": benes.routes(reverse),
    },
    "evidence_status": "EXHAUSTIVELY_ENUMERATED_BY_COMMITTED_CODE",
    "limitations": [
        "The stage-dimension model is stated explicitly; physical wire labels are in shuffle-normalized coordinates.",
        "The count is of distinct input-output permutations, not switch-setting configurations.",
        "This finite routing result does not establish a physical Q6 transport law or a Green-Kubo relation.",
    ],
}

out = Path(__file__).resolve().parents[1] / "results" / "interconnection_networks_8x8_audit.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
print(f"wrote: {out}")
