"""Exact Fourier spectrum of the finite A5-root Cayley graph.

For m=3, characters are indexed by k in (Z/3Z)^5 in simple-root
coordinates. Pairing the 30 roots into 15 opposite pairs gives an integer
adjacency eigenvalue: each pair contributes 2 when its character phase is
zero and -1 otherwise.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations, product
from math import comb, cos, pi

from src.gec6_a5_root_cayley_torus import root_step, root_steps


def positive_root_steps() -> tuple[tuple[int, ...], ...]:
    return tuple(root_step(i, j) for i, j in combinations(range(6), 2))


def fourier_laplacian_eigenvalue(character: tuple[int, ...], modulus: int = 3) -> int:
    """Exact Laplacian eigenvalue for a character of (Z/modulus Z)^5.

    The integer formula below is specific to modulus 3.
    """
    if modulus != 3:
        raise ValueError("the exact integer formula is implemented for modulus 3")
    if len(character) != 5 or any(not 0 <= x < 3 for x in character):
        raise ValueError("character must lie in (Z/3Z)^5")
    zero_phase_pairs = sum(
        sum(a * b for a, b in zip(character, root)) % 3 == 0
        for root in positive_root_steps()
    )
    # Adjacency eigenvalue = 2*a - (15-a) = 3*a-15;
    # L = 30I-A, hence lambda_L = 45-3*a.
    return 45 - 3 * zero_phase_pairs


def exact_fourier_spectrum_m3() -> dict[int, int]:
    """Return Laplacian eigenvalue multiplicities from all 3^5 characters."""
    multiplicities = Counter(
        fourier_laplacian_eigenvalue(k)
        for k in product(range(3), repeat=5)
    )
    return dict(sorted(multiplicities.items()))


def normalized_comparison() -> dict:
    """Compare normalized Laplacian spectra/gaps with Q6 and Cartesian PRISM."""
    gamma = exact_fourier_spectrum_m3()
    gamma_normalized = {
        f"{value}/30": multiplicity for value, multiplicity in gamma.items()
    }
    q6 = {f"{k}/3": comb(6, k) for k in range(7)}
    prism_gap = 2.0 - 2.0 * cos(2.0 * pi / 64.0)
    return {
        "GEC6_m3_normalized_laplacian": gamma_normalized,
        "Q6_normalized_laplacian": q6,
        "PRISM_C64xC6xC6xC6xC5_normalized_spectral_gap": prism_gap / 10.0,
        "GEC6_m3_normalized_spectral_gap": 15.0 / 30.0,
        "Q6_normalized_spectral_gap": 2.0 / 6.0,
        "warning": (
            "Different finite graphs and degrees: normalized spectra are "
            "comparable as dimensionless invariants, not evidence of graph, "
            "dynamical, or physical equivalence."
        ),
    }


def dense_adjacency(modulus: int = 3):
    """Build the exact integer adjacency matrix for optional numerical checks."""
    if modulus != 3:
        raise ValueError("the dense reference audit is fixed at modulus 3")
    import numpy as np

    vertices = list(product(range(3), repeat=5))
    index = {v: i for i, v in enumerate(vertices)}
    adjacency = np.zeros((len(vertices), len(vertices)), dtype=np.int8)
    roots = root_steps()
    for i, vertex in enumerate(vertices):
        for root in roots:
            neighbor = tuple((vertex[j] + root[j]) % 3 for j in range(5))
            adjacency[i, index[neighbor]] = 1
    return adjacency


def numerical_matrix_spectrum() -> dict[int, int]:
    """Diagonalize the 243x243 adjacency matrix and round certified integers."""
    import numpy as np

    adjacency = dense_adjacency(3)
    if not np.array_equal(adjacency, adjacency.T):
        raise AssertionError("adjacency matrix is not symmetric")
    if not np.all(adjacency.sum(axis=1) == 30):
        raise AssertionError("graph is not 30-regular")
    eigenvalues = np.linalg.eigvalsh(30 * np.eye(243) - adjacency)
    rounded = np.rint(eigenvalues).astype(int)
    if float(np.max(np.abs(eigenvalues - rounded))) > 1e-9:
        raise AssertionError("numerical eigenvalues are not close to integers")
    return dict(sorted(Counter(int(x) for x in rounded).items()))


def audit() -> dict:
    exact = exact_fourier_spectrum_m3()
    numerical = numerical_matrix_spectrum()
    if exact != numerical:
        raise AssertionError(f"Fourier/matrix mismatch: {exact} != {numerical}")
    return {
        "graph": "Cay((Z/3Z)^5, A5 roots)",
        "vertices": 243,
        "degree": 30,
        "edges": 3645,
        "laplacian_spectrum": [
            {"eigenvalue": value, "multiplicity": mult}
            for value, mult in exact.items()
        ],
        "spectrum_sum_multiplicities": sum(exact.values()),
        "spectral_gap": min(value for value in exact if value > 0),
        "normalized_spectral_gap": 0.5,
        "fourier_method": "exact character sum over all 3^5 characters",
        "matrix_check": "243x243 symmetric integer adjacency; numpy eigvalsh",
        "fourier_equals_matrix": True,
        "normalized_comparison": normalized_comparison(),
        "status": "exact Fourier multiplicities; independent numerical matrix confirmation",
    }


def normalized_fourier_eigenvalue(character: tuple[int, ...], modulus: int) -> float:
    """Normalized-Laplacian eigenvalue from the 15 positive A5 roots."""
    if modulus < 3:
        raise ValueError("modulus must be at least 3")
    if len(character) != 5 or any(not 0 <= x < modulus for x in character):
        raise ValueError("character must lie in (Z/modulus Z)^5")
    return sum(
        1.0 - cos(2.0 * pi * (sum(a * b for a, b in zip(character, root)) % modulus) / modulus)
        for root in positive_root_steps()
    ) / 15.0


def normalized_spectral_gap(modulus: int) -> float:
    """Exhaustively compute the normalized gap over all m^5 characters."""
    if modulus < 3:
        raise ValueError("modulus must be at least 3")
    zero = (0, 0, 0, 0, 0)
    return min(
        normalized_fourier_eigenvalue(k, modulus)
        for k in product(range(modulus), repeat=5)
        if k != zero
    )


def gap_scaling_audit(max_modulus: int = 10) -> dict:
    """Exhaustively compute normalized gaps and m^2-scaled values."""
    if max_modulus < 3:
        raise ValueError("max_modulus must be at least 3")
    rows = []
    for modulus in range(3, max_modulus + 1):
        gap = normalized_spectral_gap(modulus)
        rows.append({
            "modulus": modulus,
            "vertices": modulus ** 5,
            "normalized_gap": gap,
            "m2_times_gap": modulus * modulus * gap,
        })
    return {
        "graph_family": "Cay((Z/mZ)^5, A5 roots)",
        "degree": 30,
        "enumeration": "all m^5 characters; 15 positive roots; floating-point cosine evaluation",
        "rows": rows,
        "asymptotic_constant": 2.0 * pi * pi / 3.0,
        "asymptotic_statement": "m^2 times the normalized gap converges to 2*pi^2/3",
        "status": "exhaustive finite-modulus numerical audit; asymptotic constant proved separately",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2, sort_keys=True))
