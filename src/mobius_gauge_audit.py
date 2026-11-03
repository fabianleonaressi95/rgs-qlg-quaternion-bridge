"""Deterministic randomized SU(2) gauge-covariance audit on the Möbius complex.

The internal SU(2) link field is deliberately kept distinct from the
orientation-reversing transition data of the non-orientable base.
"""
from __future__ import annotations

import json
from typing import Any

import numpy as np

from src.mobius_complex import build_mobius_complex
from src.holonomies import SU2Lattice, dagger


def _random_su2(rng: np.random.Generator) -> np.ndarray:
    """Draw an SU(2) matrix from a normalized four-dimensional Gaussian."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array(
        [[a + 1j * b, c + 1j * d],
         [-c + 1j * d, a - 1j * b]],
        dtype=complex,
    )


def build_random_mobius_lattice(seed: int = 20261103) -> tuple[SU2Lattice, dict[int, np.ndarray]]:
    """Create independent random SU(2) links and vertex gauge matrices."""
    rng = np.random.default_rng(seed)
    complex_ = build_mobius_complex()
    lattice = SU2Lattice(complex_.vertices)
    for source, target in complex_.edges:
        lattice.add_link(source, target, _random_su2(rng))
    gauge = {vertex: _random_su2(rng) for vertex in complex_.vertices}
    return lattice, gauge


def audit_mobius_gauge_covariance(
    seed: int = 20261103,
    atol: float = 1e-10,
) -> dict[str, Any]:
    """Audit open-path covariance, loop conjugation and Wilson invariance."""
    complex_ = build_mobius_complex()
    lattice, gauge = build_random_mobius_lattice(seed)
    transformed = lattice.gauge_transform(gauge)

    open_residuals = []
    conjugation_residuals = []
    wilson_residuals = []
    inverse_residuals = []

    # Test every face boundary and its three-edge open prefix.
    for face in complex_.faces:
        closed_path = face + (face[0],)
        open_path = face[:]

        p = lattice.path_transport(open_path)
        p_g = transformed.path_transport(open_path)
        expected_p_g = gauge[open_path[-1]] @ p @ dagger(gauge[open_path[0]])
        open_residuals.append(float(np.linalg.norm(p_g - expected_p_g, ord="fro")))

        h = lattice.holonomy(closed_path)
        h_g = transformed.holonomy(closed_path)
        expected_h_g = gauge[closed_path[0]] @ h @ dagger(gauge[closed_path[0]])
        conjugation_residuals.append(
            float(np.linalg.norm(h_g - expected_h_g, ord="fro"))
        )
        wilson_residuals.append(
            abs(lattice.wilson_trace(closed_path) -
                transformed.wilson_trace(closed_path))
        )

    for source, target in complex_.edges:
        forward = lattice.link(source, target)
        reverse = lattice.link(target, source)
        inverse_residuals.append(
            float(np.linalg.norm(reverse - dagger(forward), ord="fro"))
        )
        transformed_forward = transformed.link(source, target)
        transformed_reverse = transformed.link(target, source)
        inverse_residuals.append(
            float(np.linalg.norm(
                transformed_reverse - dagger(transformed_forward), ord="fro"
            ))
        )

    maxima = {
        "open_path_covariance_max": max(open_residuals, default=0.0),
        "closed_loop_conjugation_max": max(conjugation_residuals, default=0.0),
        "wilson_trace_invariance_max": max(wilson_residuals, default=0.0),
        "inverse_link_consistency_max": max(inverse_residuals, default=0.0),
    }
    passed = all(value <= atol for value in maxima.values())
    return {
        "audit": "MOBIUS_SU2_GAUGE_COVARIANCE",
        "seed": int(seed),
        "tolerance": float(atol),
        "vertices": len(complex_.vertices),
        "edges": len(complex_.edges),
        "faces_tested": len(complex_.faces),
        "paths_tested": len(complex_.faces),
        **maxima,
        "passed": bool(passed),
    }


if __name__ == "__main__":
    print(json.dumps(audit_mobius_gauge_covariance(), indent=2, sort_keys=True))
