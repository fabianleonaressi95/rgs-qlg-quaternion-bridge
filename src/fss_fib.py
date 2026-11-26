"""Finite-size response tools for the Fibonacci Aubry–André chain."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]

def fibonacci_pair(k: int) -> tuple[int, int]:
    """Return (F_k, F_{k-1}) with F_0=0, F_1=1."""
    if not isinstance(k, int) or k < 2:
        raise ValueError("k must be an integer >= 2")
    a, b = 0, 1
    seq = [a, b]
    for _ in range(2, k + 1):
        a, b = b, a + b
        seq.append(b)
    return seq[k], seq[k - 1]

def fibonacci_approximant(k: int) -> tuple[int, int, float]:
    n, p = fibonacci_pair(k)
    return n, p, p / n

def aubry_andre_hamiltonian(n: int, *, hopping: float = 1.0,
    potential: float | None = None, alpha: float, phase: float = 0.0,
    periodic: bool = False) -> FloatArray:
    """H=-t(T+T*)+V cos(2π alpha j+phase); default V=2t."""
    if n < 2 or hopping <= 0 or not np.isfinite(hopping):
        raise ValueError("n >= 2 and finite positive hopping required")
    if not np.isfinite(alpha) or not np.isfinite(phase):
        raise ValueError("alpha and phase must be finite")
    v = 2.0 * hopping if potential is None else float(potential)
    if not np.isfinite(v):
        raise ValueError("potential must be finite")
    j = np.arange(n, dtype=float)
    h = np.diag(v * np.cos(2.0 * np.pi * alpha * j + phase))
    idx = np.arange(n - 1)
    h[idx, idx + 1] = -hopping
    h[idx + 1, idx] = -hopping
    if periodic:
        h[0, n - 1] = h[n - 1, 0] = -hopping
    return h

def diagonalize_chain(h: ArrayLike) -> tuple[FloatArray, FloatArray]:
    matrix = np.asarray(h, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("h must be square")
    if not np.allclose(matrix, matrix.T, atol=1e-12, rtol=1e-12):
        raise ValueError("h must be real symmetric")
    return np.linalg.eigh(matrix)

def density_vertex(q: float, n: int, positions: ArrayLike | None = None) -> ComplexArray:
    if n < 1 or not np.isfinite(q):
        raise ValueError("n must be positive and q finite")
    x = np.arange(n, dtype=float) if positions is None else np.asarray(positions, dtype=float)
    if x.shape != (n,) or not np.all(np.isfinite(x)):
        raise ValueError("positions must be a finite vector of length n")
    return np.diag(np.exp(1j * q * x)).astype(np.complex128)

def occupations_zero_temperature(eigenvalues: ArrayLike, particle_number: int) -> FloatArray:
    e = np.asarray(eigenvalues, dtype=float)
    if e.ndim != 1 or not np.all(np.isfinite(e)):
        raise ValueError("eigenvalues must be a finite vector")
    if not 0 <= particle_number <= e.size:
        raise ValueError("particle_number out of range")
    f = np.zeros(e.size)
    f[np.argsort(e)[:particle_number]] = 1.0
    return f

def occupations_fermi_dirac(eigenvalues: ArrayLike, *, temperature: float,
    chemical_potential: float) -> FloatArray:
    e = np.asarray(eigenvalues, dtype=float)
    if e.ndim != 1 or not np.all(np.isfinite(e)):
        raise ValueError("eigenvalues must be a finite vector")
    if temperature <= 0 or not np.isfinite(temperature) or not np.isfinite(chemical_potential):
        raise ValueError("finite T>0 and finite chemical potential required")
    x = np.clip((e - chemical_potential) / temperature, -700.0, 700.0)
    return 1.0 / (np.exp(x) + 1.0)

def chemical_potential_for_filling(eigenvalues: ArrayLike, *, temperature: float,
    particle_number: float) -> float:
    e = np.asarray(eigenvalues, dtype=float)
    if e.ndim != 1 or e.size == 0 or not np.all(np.isfinite(e)):
        raise ValueError("eigenvalues must be a non-empty finite vector")
    if temperature <= 0 or not np.isfinite(temperature) or not 0 < particle_number < e.size:
        raise ValueError("require T>0 and 0<particle_number<N")
    lo, hi = float(e.min() - 100*temperature - 1), float(e.max() + 100*temperature + 1)
    def count(mu):
        return float(np.sum(occupations_fermi_dirac(e, temperature=temperature, chemical_potential=mu)))
    return float(brentq(lambda mu: count(mu) - particle_number, lo, hi))

def _response_terms(eigenvalues, eigenvectors, q, occupations):
    e = np.asarray(eigenvalues, dtype=float)
    u = np.asarray(eigenvectors, dtype=np.complex128)
    f = np.asarray(occupations, dtype=float)
    n = e.size
    if e.ndim != 1 or u.shape != (n, n) or f.shape != (n,):
        raise ValueError("incompatible eigenvalue, eigenvector, or occupation shapes")
    if not np.all(np.isfinite(e)) or not np.all(np.isfinite(f)) or np.any((f < 0) | (f > 1)):
        raise ValueError("finite eigenvalues and occupations in [0,1] required")
    rho = density_vertex(q, n)
    m = u.conj().T @ rho @ u
    weights = (f[:, None] - f[None, :]) * np.abs(m)**2 / n
    return weights, e[:, None] - e[None, :]

def chi0_retarded(eigenvalues, eigenvectors, q: float, omega: complex,
    occupations, *, eta: float = 0.0) -> complex:
    """χ0^R(q,z)=N^-1 Σ_mn (f_m-f_n)|<m|rho_q|n>|²/(z+E_m-E_n)."""
    if eta < 0 or not np.isfinite(eta):
        raise ValueError("eta must be finite and non-negative")
    z = complex(omega) + 1j*eta
    if z.imag < 0:
        raise ValueError("retarded response requires Im(z)>=0")
    weights, de = _response_terms(eigenvalues, eigenvectors, q, occupations)
    den = z + de
    if np.any(np.abs(den) < 1e-14):
        raise ZeroDivisionError("frequency lies on an unregularized discrete pole")
    return complex(np.sum(weights / den))

def chi0_derivative(eigenvalues, eigenvectors, q: float, omega: complex,
    occupations) -> complex:
    """Analytic derivative dχ0/dz away from discrete poles."""
    z = complex(omega)
    if z.imag < 0:
        raise ValueError("retarded response requires Im(z)>=0")
    weights, de = _response_terms(eigenvalues, eigenvectors, q, occupations)
    den = z + de
    if np.any(np.abs(den) < 1e-14):
        raise ZeroDivisionError("derivative evaluated on a discrete pole")
    return complex(-np.sum(weights / den**2))

def nn_density_kernel(q: float, interaction: float) -> float:
    """Direct density-channel kernel for U sum_j n_j n_(j+1): 2U cos(q)."""
    if not np.isfinite(q) or not np.isfinite(interaction):
        raise ValueError("q and interaction must be finite")
    return float(2.0 * interaction * np.cos(q))

def chi_rpa(chi0: complex, interaction_kernel: float) -> complex:
    den = 1.0 - interaction_kernel*chi0
    if abs(den) < 1e-14:
        raise ZeroDivisionError("RPA denominator is at a pole")
    return complex(chi0 / den)

def rpa_pole_residue(chi0_at_pole: complex, dchi0_dz: complex,
    interaction_kernel: float) -> complex:
    """Residue of a simple pole: -1/[V(q)^2 χ0'(z_L)]."""
    if interaction_kernel == 0 or abs(dchi0_dz) < 1e-14:
        raise ValueError("nonzero interaction and simple-pole derivative required")
    if abs(1.0 - interaction_kernel*chi0_at_pole) > 1e-5:
        raise ValueError("point does not satisfy the RPA pole condition")
    return complex(-1.0 / (interaction_kernel**2 * dchi0_dz))

@dataclass(frozen=True)
class FSSFit:
    limit: float
    amplitude: float
    exponent: float
    residual_norm: float
    n_points: int

def fit_fss_power_law(sizes: Iterable[int], values: Iterable[float], *,
    exponent: float) -> FSSFit:
    """Descriptive fit y(N)=y_inf+a*N^-p with p fixed independently."""
    n = np.asarray(tuple(sizes), dtype=float)
    y = np.asarray(tuple(values), dtype=float)
    if n.ndim != 1 or y.shape != n.shape or n.size < 2:
        raise ValueError("equal-length size/value vectors with >=2 points required")
    if np.any(n <= 0) or not np.all(np.isfinite(n)) or not np.all(np.isfinite(y)):
        raise ValueError("positive sizes and finite values required")
    if exponent <= 0 or not np.isfinite(exponent):
        raise ValueError("exponent must be finite and positive")
    design = np.column_stack((np.ones(n.size), n**(-exponent)))
    coeff, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coeff
    return FSSFit(float(coeff[0]), float(coeff[1]), float(exponent),
                  float(np.linalg.norm(residual)), int(n.size))
