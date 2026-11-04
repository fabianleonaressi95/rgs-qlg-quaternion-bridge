"""Qutrit (three-level) entanglement on a Fibonacci-aperiodic SU(3) chain.

Each site carries the fundamental C^3 representation. Neighbouring sites are
coupled by the SU(3)-invariant swap operator P_{i,i+1}, with bond strengths
modulated by a finite prefix of the Fibonacci word A->AB, B->A.

This is a finite-dimensional quantum-spin model. It is not, by itself, a
continuum QFT, a QCD colour model, or evidence for physical aperiodicity.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def fibonacci_word(length: int) -> str:
    """Return the requested prefix of the fixed point A->AB, B->A."""
    if length < 0:
        raise ValueError("length must be non-negative")
    word = "A"
    while len(word) < length:
        word = word.replace("B", "a").replace("A", "AB").replace("a", "A")
    return word[:length]


def aperiodic_prefix_audit(word: str) -> dict[str, object]:
    """Check all candidate periods of a finite word; this is not an infinite proof."""
    if not word:
        raise ValueError("word must be non-empty")
    periods = [
        p for p in range(1, len(word))
        if all(word[i] == word[i - p] for i in range(p, len(word)))
    ]
    return {
        "length": len(word),
        "is_periodic_as_finite_word": bool(periods),
        "finite_word_periods": periods,
        "smallest_period": min(periods) if periods else None,
        "scope": "finite-prefix test only; infinite Fibonacci-word aperiodicity is a separate theorem",
    }


def qutrit_bell_state() -> np.ndarray:
    """Return |Phi_3> = (|00>+|11>+|22>)/sqrt(3), shape (3,3)."""
    return np.eye(3, dtype=complex).reshape(9) / np.sqrt(3.0)


def schmidt_spectrum(state: np.ndarray, dims: tuple[int, int]) -> np.ndarray:
    """Return normalized Schmidt probabilities for a bipartite pure state."""
    psi = np.asarray(state, dtype=complex)
    da, db = dims
    if da < 1 or db < 1 or psi.size != da * db:
        raise ValueError("state size must equal the product of positive subsystem dimensions")
    if not np.all(np.isfinite(psi)):
        raise ValueError("state must be finite")
    norm = np.linalg.norm(psi)
    if norm <= 0:
        raise ValueError("state must be nonzero")
    singular = np.linalg.svd(psi.reshape(da, db) / norm, compute_uv=False)
    probs = singular**2
    probs /= probs.sum()
    return probs


def von_neumann_entropy(probabilities: np.ndarray, base: float = 2.0) -> float:
    """Shannon/von-Neumann entropy of a normalized non-negative spectrum."""
    p = np.asarray(probabilities, dtype=float)
    if p.ndim != 1 or np.any(p < -1e-12) or not np.all(np.isfinite(p)):
        raise ValueError("probabilities must be a finite non-negative vector")
    p = np.clip(p, 0.0, None)
    total = p.sum()
    if total <= 0 or base <= 0 or np.isclose(base, 1.0):
        raise ValueError("positive spectrum and valid logarithm base required")
    p = p / total
    nz = p > 0
    return float(-np.sum(p[nz] * np.log(p[nz]) / np.log(base)))


def _swap_operator(n_sites: int, left_site: int) -> np.ndarray:
    """Dense permutation matrix swapping adjacent base-3 tensor factors."""
    dim = 3**n_sites
    out = np.zeros((dim, dim), dtype=float)
    for source in range(dim):
        digits = [0] * n_sites
        value = source
        for k in range(n_sites - 1, -1, -1):
            digits[k] = value % 3
            value //= 3
        digits[left_site], digits[left_site + 1] = digits[left_site + 1], digits[left_site]
        target = 0
        for digit in digits:
            target = 3 * target + digit
        out[target, source] = 1.0
    return out


@dataclass(frozen=True)
class FibonacciQutritChain:
    n_sites: int
    word: str
    couplings: tuple[float, ...]
    hamiltonian: np.ndarray


def build_fibonacci_qutrit_chain(
    n_sites: int,
    j_a: float = 1.0,
    j_b: float = 0.5,
    max_sites: int = 7,
) -> FibonacciQutritChain:
    """Build H=sum_i J_i P_{i,i+1} on an open qutrit chain.

    Positive J_A/J_B give the antiferromagnetic permutation-chain convention.
    Dense exact diagonalisation is intentionally capped at max_sites.
    """
    if not isinstance(n_sites, int) or n_sites < 2 or n_sites > max_sites:
        raise ValueError(f"n_sites must be an integer in [2, {max_sites}]")
    if not np.isfinite(j_a) or not np.isfinite(j_b) or j_a <= 0 or j_b <= 0:
        raise ValueError("j_a and j_b must be finite and strictly positive")
    word = fibonacci_word(n_sites - 1)
    couplings = tuple(float(j_a if letter == "A" else j_b) for letter in word)
    dim = 3**n_sites
    hamiltonian = np.zeros((dim, dim), dtype=float)
    for bond, coupling in enumerate(couplings):
        hamiltonian += coupling * _swap_operator(n_sites, bond)
    return FibonacciQutritChain(n_sites, word, couplings, hamiltonian)


def ground_state_entanglement(
    chain: FibonacciQutritChain,
    cut: int | None = None,
) -> dict[str, object]:
    """Diagonalise the finite chain and measure pure-state entanglement across a cut."""
    n = chain.n_sites
    cut = n // 2 if cut is None else cut
    if not 1 <= cut < n:
        raise ValueError("cut must split the chain into two non-empty subsystems")
    h = np.asarray(chain.hamiltonian, dtype=float)
    if h.shape != (3**n, 3**n) or not np.allclose(h, h.T, atol=1e-12):
        raise ValueError("chain Hamiltonian has invalid shape or is not Hermitian")
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    ground_energy = float(eigenvalues[0])
    degeneracy_tol = 1e-10 * max(1.0, abs(ground_energy))
    ground_degeneracy = int(np.count_nonzero(eigenvalues - eigenvalues[0] <= degeneracy_tol))
    state = eigenvectors[:, 0].astype(complex)
    probs = schmidt_spectrum(state, (3**cut, 3**(n-cut)))
    spectral_gap = (
        float(eigenvalues[ground_degeneracy] - eigenvalues[0])
        if ground_degeneracy < len(eigenvalues) else 0.0
    )
    return {
        "ground_energy": ground_energy,
        "ground_state_degeneracy": ground_degeneracy,
        "gap_above_ground_manifold": spectral_gap,
        "gap": spectral_gap,
        "schmidt_probabilities": probs.tolist(),
        "entanglement_entropy_bits": von_neumann_entropy(probs),
        "schmidt_rank": int(np.count_nonzero(probs > 1e-12)),
        "cut": cut,
        "left_dimension": 3**cut,
        "right_dimension": 3**(n-cut),
        "word": chain.word,
        "bond_couplings": list(chain.couplings),
        "hamiltonian_hermiticity_residual": float(np.linalg.norm(h-h.T)),
        "ground_state_entanglement_basis_dependent": bool(ground_degeneracy > 1),
        "scope": "finite-chain selected ground vector; entropy can depend on basis inside a degenerate ground space; no thermodynamic or continuum claim",
    }


def audit_qutrit_aperiodic_model(n_sites: int = 4, j_a: float = 1.0, j_b: float = 0.5) -> dict:
    """Combined reproducibility report for the qutrit Bell state and Fibonacci chain."""
    bell = qutrit_bell_state()
    bell_probs = schmidt_spectrum(bell, (3, 3))
    chain = build_fibonacci_qutrit_chain(n_sites, j_a, j_b)
    ground = ground_state_entanglement(chain)
    return {
        "local_hilbert_dimension": 3,
        "local_basis": ["|0>", "|1>", "|2>"],
        "bell_state_schmidt_probabilities": bell_probs.tolist(),
        "bell_state_entropy_bits": von_neumann_entropy(bell_probs),
        "fibonacci_word": chain.word,
        "aperiodic_prefix_audit": aperiodic_prefix_audit(chain.word),
        "ground_state": ground,
    }
