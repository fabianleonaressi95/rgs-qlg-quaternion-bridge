"""Gauge-covariant qutrit chain embedded along a Gray path in Q6.

The finite path is a subgraph of the canonical six-cube. Each active Q6
vertex carries C^3. Fibonacci letters modulate nearest-neighbour exchange.
Canonical edge links U_{v<-u} (u<v) live in SU(3); the transported swap is
(U_{v<-u}^dagger tensor U_{v<-u}) P and is covariant under independent
site transformations. Dense exact diagonalisation is restricted to small
active paths.

This is a finite lattice quantum-spin construction, not continuum QCD.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from src.q6_complex import build_q6_complex
from src.qutrit_aperiodic import fibonacci_word


def _validate_su3(matrix: np.ndarray, atol: float = 2e-10) -> np.ndarray:
    value = np.asarray(matrix, dtype=complex)
    if value.shape != (3, 3) or not np.all(np.isfinite(value)):
        raise ValueError("each SU(3) link/gauge matrix must be finite and 3x3")
    if not np.allclose(value.conj().T @ value, np.eye(3), atol=atol, rtol=atol):
        raise ValueError("matrix is not unitary")
    if not np.isclose(np.linalg.det(value), 1.0 + 0.0j, atol=atol, rtol=atol):
        raise ValueError("matrix determinant must equal one")
    return value.copy()


def q6_gray_path(n_sites: int) -> tuple[int, ...]:
    """Return n_sites consecutive vertices of a Hamiltonian Gray path in Q6."""
    if not isinstance(n_sites, int) or not 2 <= n_sites <= 7:
        raise ValueError("n_sites must be an integer in [2, 7] for dense exact diagonalisation")
    return tuple(i ^ (i >> 1) for i in range(n_sites))


def _swap(n_sites: int, bond: int) -> np.ndarray:
    dimension = 3**n_sites
    result = np.zeros((dimension, dimension), dtype=complex)
    for source in range(dimension):
        digits = [0] * n_sites
        value = source
        for k in range(n_sites - 1, -1, -1):
            digits[k] = value % 3
            value //= 3
        digits[bond], digits[bond + 1] = digits[bond + 1], digits[bond]
        target = 0
        for digit in digits:
            target = 3 * target + digit
        result[target, source] = 1.0
    return result


def _oriented_link(
    source: int,
    target: int,
    links: Mapping[tuple[int, int], np.ndarray],
) -> np.ndarray:
    key = (min(source, target), max(source, target))
    canonical = _validate_su3(links[key])
    return canonical if source < target else canonical.conj().T


@dataclass(frozen=True)
class Q6QutritChain:
    n_sites: int
    site_vertices: tuple[int, ...]
    word: str
    couplings: tuple[float, ...]
    links: dict[tuple[int, int], np.ndarray]
    hamiltonian: np.ndarray


def build_q6_qutrit_chain(
    n_sites: int,
    j_a: float = 1.0,
    j_b: float = 0.5,
    links: Mapping[tuple[int, int], np.ndarray] | None = None,
) -> Q6QutritChain:
    """Build Fibonacci-modulated qutrit exchange on a Gray path in Q6.

    For a path edge source->target, U denotes U_{target<-source}. Its
    covariant swap is (U^dagger on source tensor U on target) times P.
    """
    vertices = q6_gray_path(n_sites)
    if not np.isfinite(j_a) or not np.isfinite(j_b) or j_a <= 0 or j_b <= 0:
        raise ValueError("j_a and j_b must be finite and strictly positive")
    q6_edges = set(build_q6_complex().edges)
    path_edges = tuple((min(a, b), max(a, b)) for a, b in zip(vertices[:-1], vertices[1:]))
    if any(edge not in q6_edges for edge in path_edges):
        raise ArithmeticError("Gray path contains a non-edge of Q6")

    word = fibonacci_word(n_sites - 1)
    couplings = tuple(float(j_a if letter == "A" else j_b) for letter in word)
    identity = np.eye(3, dtype=complex)
    supplied = {} if links is None else dict(links)
    if not set(supplied).issubset(set(path_edges)):
        raise ValueError("link keys must be active canonical edges of the embedded Q6 path")
    link_data = {edge: _validate_su3(supplied.get(edge, identity)) for edge in path_edges}

    dimension = 3**n_sites
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for bond, (source, target) in enumerate(zip(vertices[:-1], vertices[1:])):
        transporter = _oriented_link(source, target, link_data)
        factors = [identity] * n_sites
        factors[bond] = transporter.conj().T
        factors[bond + 1] = transporter
        transport_operator = factors[0]
        for factor in factors[1:]:
            transport_operator = np.kron(transport_operator, factor)
        hamiltonian += couplings[bond] * (transport_operator @ _swap(n_sites, bond))

    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=2e-10, rtol=2e-10):
        raise ArithmeticError("gauge-covariant exchange Hamiltonian is not Hermitian")
    return Q6QutritChain(n_sites, vertices, word, couplings, link_data, hamiltonian)


def gauge_transform_q6_links(
    chain: Q6QutritChain,
    gauge: Mapping[int, np.ndarray],
) -> dict[tuple[int, int], np.ndarray]:
    """Transform canonical links U_{v<-u}->g_v U g_u^dagger."""
    if set(gauge) != set(chain.site_vertices):
        raise ValueError("gauge map must contain exactly the active Q6 path vertices")
    matrices = {v: _validate_su3(gauge[v]) for v in chain.site_vertices}
    transformed = {}
    for (u, v), link in chain.links.items():
        transformed[(u, v)] = matrices[v] @ link @ matrices[u].conj().T
    return transformed


def q6_qutrit_gauge_covariance_residual(
    chain: Q6QutritChain,
    gauge: Mapping[int, np.ndarray],
) -> float:
    """Return ||H[U^g] - G H[U] G^dagger||_F for the active path."""
    transformed_links = gauge_transform_q6_links(chain, gauge)
    transformed = build_q6_qutrit_chain(
        chain.n_sites,
        j_a=chain.couplings[0],
        j_b=chain.couplings[1] if len(chain.couplings) > 1 else chain.couplings[0],
        links=transformed_links,
    )
    local = [_validate_su3(gauge[v]) for v in chain.site_vertices]
    global_gauge = local[0]
    for matrix in local[1:]:
        global_gauge = np.kron(global_gauge, matrix)
    expected = global_gauge @ chain.hamiltonian @ global_gauge.conj().T
    return float(np.linalg.norm(transformed.hamiltonian - expected))


def q6_qutrit_spectrum(chain: Q6QutritChain, cut: int | None = None) -> dict[str, object]:
    """Return finite-chain energies and Schmidt entropy of one ground vector."""
    n = chain.n_sites
    cut = n // 2 if cut is None else cut
    if not 1 <= cut < n:
        raise ValueError("cut must split the chain into two non-empty subsystems")
    eigenvalues, eigenvectors = np.linalg.eigh(chain.hamiltonian)
    ground_energy = float(eigenvalues[0])
    tol = 1e-10 * max(1.0, abs(ground_energy))
    degeneracy = int(np.count_nonzero(eigenvalues - eigenvalues[0] <= tol))
    state = eigenvectors[:, 0].reshape(3**cut, 3**(n-cut))
    singular = np.linalg.svd(state, compute_uv=False)
    probabilities = singular**2
    probabilities /= probabilities.sum()
    nonzero = probabilities > 1e-14
    entropy = float(-np.sum(probabilities[nonzero] * np.log2(probabilities[nonzero])))
    gap = float(eigenvalues[degeneracy] - eigenvalues[0]) if degeneracy < len(eigenvalues) else 0.0
    return {
        "ground_energy": ground_energy,
        "ground_degeneracy": degeneracy,
        "gap_above_ground_manifold": gap,
        "schmidt_probabilities": probabilities.tolist(),
        "entanglement_entropy_bits": entropy,
        "schmidt_rank": int(np.count_nonzero(nonzero)),
        "cut": cut,
        "site_vertices": list(chain.site_vertices),
        "fibonacci_word": chain.word,
        "bond_couplings": list(chain.couplings),
        "hermiticity_residual": float(np.linalg.norm(chain.hamiltonian-chain.hamiltonian.conj().T)),
        "entropy_basis_dependent_in_degenerate_ground_space": bool(degeneracy > 1),
        "scope": "finite Q6-embedded qutrit path; selected ground vector; no continuum/QCD claim",
    }
