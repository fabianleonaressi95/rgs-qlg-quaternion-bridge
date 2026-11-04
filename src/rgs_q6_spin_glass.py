"""Edwards–Anderson Ising spin glass on the canonical Q6 graph.

The quenched bonds are independent of the SU(2) connection variables used by
other Q6 modules.  This is a finite-graph statistical-mechanics model, not a
claim of a thermodynamic spin-glass phase.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from src.q6_complex import build_q6_complex


@dataclass(frozen=True)
class Q6SpinGlass:
    edges: np.ndarray
    faces: np.ndarray
    couplings: np.ndarray
    n_vertices: int = 64

    @classmethod
    def random(
        cls, seed: int = 20261104, distribution: str = "gaussian",
        scale: float = 1.0,
    ) -> "Q6SpinGlass":
        """Draw one quenched disorder realization and keep it fixed."""
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("scale must be finite and positive")
        q6 = build_q6_complex()
        rng = np.random.default_rng(seed)
        if distribution == "gaussian":
            bonds = rng.normal(0.0, scale, len(q6.edges))
        elif distribution == "bimodal":
            bonds = scale * rng.choice(np.array([-1.0, 1.0]), len(q6.edges))
        else:
            raise ValueError("distribution must be 'gaussian' or 'bimodal'")
        return cls(np.asarray(q6.edges, dtype=np.int16),
                   np.asarray(q6.faces, dtype=np.int16), bonds.astype(float))

    def __post_init__(self) -> None:
        e = np.asarray(self.edges)
        f = np.asarray(self.faces)
        j = np.asarray(self.couplings, dtype=float)
        if e.shape != (192, 2) or f.shape != (240, 4) or j.shape != (192,):
            raise ValueError("Q6 requires edges (192,2), faces (240,4), bonds (192,)")
        if not np.all(np.isfinite(j)):
            raise ValueError("couplings must be finite")

    def energy(self, spins: np.ndarray) -> float:
        s = _validate_spins(spins, self.n_vertices)
        return float(-np.sum(self.couplings * s[self.edges[:, 0]] * s[self.edges[:, 1]]))

    def local_fields(self, spins: np.ndarray) -> np.ndarray:
        s = _validate_spins(spins, self.n_vertices)
        fields = np.zeros(self.n_vertices, dtype=float)
        np.add.at(fields, self.edges[:, 0], self.couplings * s[self.edges[:, 1]])
        np.add.at(fields, self.edges[:, 1], self.couplings * s[self.edges[:, 0]])
        return fields

    def frustration(self) -> np.ndarray:
        """Gauge-invariant Z2 flux product sign(J) around each square."""
        signs = np.sign(self.couplings)
        if np.any(signs == 0):
            raise ValueError("plaquette frustration is undefined for zero bonds")
        edge_index = {tuple(map(int, edge)): k for k, edge in enumerate(self.edges)}
        out = np.empty(len(self.faces), dtype=np.int8)
        for p, (a, b, c, d) in enumerate(self.faces):
            boundary = ((a, b), (b, c), (c, d), (d, a))
            product = 1
            for u, v in boundary:
                product *= int(signs[edge_index[(min(int(u), int(v)), max(int(u), int(v))) ]])
            out[p] = product
        return out

    def gauge_transform(self, vertex_signs: np.ndarray) -> "Q6SpinGlass":
        """Z2 gauge transformation J_ij -> eps_i J_ij eps_j."""
        eps = _validate_spins(vertex_signs, self.n_vertices)
        bonds = self.couplings * eps[self.edges[:, 0]] * eps[self.edges[:, 1]]
        return Q6SpinGlass(self.edges.copy(), self.faces.copy(), bonds)


def _validate_spins(spins: np.ndarray, n: int) -> np.ndarray:
    s = np.asarray(spins, dtype=np.int8)
    if s.shape != (n,) or not np.all(np.isin(s, (-1, 1))):
        raise ValueError(f"spins must be a length-{n} array containing only +/-1")
    return s


def exact_ground_state(model: Q6SpinGlass, max_vertices: int = 24) -> dict:
    """Exact Gray-code enumeration for small induced models only.

    The canonical 64-site Q6 is deliberately rejected: 2**64 enumeration is
    infeasible.  Use parallel tempering for the full graph.
    """
    n = model.n_vertices
    if n > max_vertices:
        raise ValueError(
            f"exact enumeration refused for N={n}; limit is {max_vertices}. "
            "Use parallel_tempering() for Q6."
        )
    best_e = np.inf
    best_s = None
    spins = -np.ones(n, dtype=np.int8)
    current = model.energy(spins)
    best_e, best_s = current, spins.copy()
    previous_gray = 0
    for k in range(1, 1 << n):
        gray = k ^ (k >> 1)
        changed = gray ^ previous_gray
        bit = (changed & -changed).bit_length() - 1
        # Flip one spin; update energy from its local field before the flip.
        fields = model.local_fields(spins)
        current += 2.0 * spins[bit] * fields[bit]
        spins[bit] *= -1
        if current < best_e:
            best_e, best_s = current, spins.copy()
        previous_gray = gray
    return {"energy": float(best_e), "spins": best_s, "enumerated_states": 1 << n}


def _one_parallel_tempering_ladder(
    model: Q6SpinGlass, temperatures: np.ndarray, sweeps: int,
    burn_in: int, measure_every: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    n_t, n = len(temperatures), model.n_vertices
    beta = 1.0 / temperatures
    spins = rng.choice(np.array([-1, 1], dtype=np.int8), size=(n_t, n))
    energies = np.array([model.energy(s) for s in spins], dtype=float)
    n_meas = (sweeps - burn_in) // measure_every
    energy_samples = np.empty((n_t, n_meas), dtype=float)
    spin_samples = np.empty((n_t, n_meas, n), dtype=np.int8)
    measurement = 0
    for sweep in range(sweeps):
        for t in range(n_t):
            for site in rng.permutation(n):
                s = spins[t]
                field = model.local_fields(s)[site]
                delta = 2.0 * s[site] * field
                if delta <= 0.0 or rng.random() < np.exp(-beta[t] * delta):
                    s[site] *= -1
                    energies[t] += delta
        parity = sweep % 2
        for t in range(parity, n_t - 1, 2):
            log_accept = (beta[t] - beta[t + 1]) * (energies[t] - energies[t + 1])
            if log_accept >= 0.0 or np.log(rng.random()) < log_accept:
                spins[[t, t + 1]] = spins[[t + 1, t]]
                energies[[t, t + 1]] = energies[[t + 1, t]]
        if sweep >= burn_in and (sweep - burn_in) % measure_every == 0:
            energy_samples[:, measurement] = energies
            spin_samples[:, measurement, :] = spins
            measurement += 1
    return energy_samples, spin_samples


def parallel_tempering(
    model: Q6SpinGlass, temperatures: Iterable[float] | None = None,
    sweeps: int = 2000, burn_in: int = 1000, measure_every: int = 10,
    seed: int = 20261105,
) -> dict:
    """Run two independent PT ladders for same-disorder replica overlaps."""
    temps = np.asarray(
        list(temperatures) if temperatures is not None else np.geomspace(0.5, 3.0, 8),
        dtype=float,
    )
    if temps.ndim != 1 or len(temps) < 2 or np.any(~np.isfinite(temps)) or np.any(temps <= 0):
        raise ValueError("temperatures must contain at least two finite positive values")
    if np.any(np.diff(temps) <= 0):
        raise ValueError("temperatures must be strictly increasing")
    if sweeps <= burn_in or burn_in < 0 or measure_every < 1:
        raise ValueError("require sweeps > burn_in >= 0 and measure_every >= 1")
    n_meas = (sweeps - burn_in) // measure_every
    if n_meas < 2:
        raise ValueError("at least two measurements are required")
    root = np.random.SeedSequence(seed)
    child = root.spawn(2)
    ladders = [
        _one_parallel_tempering_ladder(
            model, temps, sweeps, burn_in, measure_every,
            np.random.default_rng(child[k]),
        ) for k in range(2)
    ]
    energies = np.stack([x[0] for x in ladders])
    configs = np.stack([x[1] for x in ladders])
    overlaps = np.mean(configs[0] * configs[1], axis=2)  # (T, measurement)
    q2 = np.mean(overlaps**2, axis=1)
    q4 = np.mean(overlaps**4, axis=1)
    binder = np.full(len(temps), np.nan)
    valid = q2 > 0
    binder[valid] = 1.0 - q4[valid] / (3.0 * q2[valid]**2)
    return {
        "audit": "RGS_Q6_EDWARDS_ANDERSON_PARALLEL_TEMPERING",
        "model": "H=-sum_(ij) J_ij s_i s_j; quenched bonds fixed during thermal sampling",
        "n_vertices": model.n_vertices,
        "n_edges": int(len(model.edges)),
        "n_plaquettes": int(len(model.faces)),
        "temperatures": temps.tolist(),
        "sweeps": int(sweeps),
        "burn_in": int(burn_in),
        "measurements": int(n_meas),
        "mean_energy_per_spin": np.mean(energies, axis=(0, 2)).__truediv__(model.n_vertices).tolist(),
        "replica_overlap_mean": np.mean(overlaps, axis=1).tolist(),
        "replica_overlap_second_moment": q2.tolist(),
        "replica_overlap_fourth_moment": q4.tolist(),
        "spin_glass_susceptibility": (model.n_vertices * q2).tolist(),
        "binder_cumulant": binder.tolist(),
        "frustrated_plaquettes": int(np.count_nonzero(model.frustration() < 0)),
        "frustration_fraction": float(np.mean(model.frustration() < 0)),
        "independent_ladders": 2,
        "quenched_disorder_seed_is_external": True,
        "limitations": [
            "finite N=64 graph; no thermodynamic-limit phase claim",
            "parallel tempering output is exploratory; equilibration requires diagnostics",
            "two replicas estimate overlap moments but do not establish replica-symmetry breaking",
            "no continuum or PSU coupling is assumed",
        ],
    }
