"""Koopman bridge for the RGS–PSU–QLG Cantor/Q6 and barycentric dynamics.

Conventions
-----------
* Markov transition matrices are row-stochastic: (K f)(i) = sum_j P[i,j] f(j).
* Probability columns evolve by mu_next = P.T @ mu.
* EDMD uses row-wise feature snapshots and returns K with Phi_next ~= Phi @ K.
* The Cantor–Parry/Q6 chain is a finite symbolic model, not yet the physical
  RGS flow. Its geometry is supplied separately through a seed-derived chart.
"""
from __future__ import annotations

import numpy as np

PSI = 1.465571231876768
PHI = (1.0 + np.sqrt(5.0)) / 2.0


def _parry_data():
    A = np.array([[0.0, 0.0, 1.0],
                  [1.0, 0.0, 0.0],
                  [0.0, 1.0, 1.0]])
    vals, vecs = np.linalg.eig(A)
    idx = int(np.argmax(np.real(vals)))
    psi = float(np.real(vals[idx]))
    r = np.real(vecs[:, idx])
    if r.sum() < 0:
        r = -r
    # Fix the Perron normalization deterministically.
    r = r / r[0]
    P = A * r[None, :] / (psi * r[:, None])
    vals_l, vecs_l = np.linalg.eig(P.T)
    j = int(np.argmin(np.abs(vals_l - 1.0)))
    pi = np.real(vecs_l[:, j])
    if pi.sum() < 0:
        pi = -pi
    pi = pi / pi.sum()
    return A, psi, r, P, pi


def build_cantor_q6_markov_chain():
    """Return the directed 192-state Parry x Q6 chain and its invariant law.

    A state is (i,g), i in the three-state supergolden Parry automaton and
    g in {0,...,63}. For each allowed i->j transition, an emitted bit b is
    drawn with the golden-weighted probability and flips Q6 bit (i+j) mod 6.
    The chain retains orientation; it is not the symmetrized conductance
    Laplacian used by the static spectral operator.
    """
    A, psi, r, P3, pi3 = _parry_data()
    qphi = np.exp(-1.0 / PHI**2) / PHI
    q = qphi * r / np.prod(r)**(1.0 / 3.0)
    p1 = q / (1.0 + q)
    p0 = 1.0 / (1.0 + q)
    T = np.zeros((192, 192), dtype=float)
    for i in range(3):
        for g in range(64):
            source = 64 * i + g
            for j in range(3):
                if A[i, j] == 0:
                    continue
                bit = (i + j) % 6
                for b, pb in ((0, p0[i]), (1, p1[i])):
                    target_g = g ^ (int(b) << bit)
                    target = 64 * j + target_g
                    T[source, target] += P3[i, j] * pb
    invariant = np.repeat(pi3 / 64.0, 64)
    return T, invariant, {
        "adjacency": A,
        "psi": psi,
        "perron_vector": r,
        "parry_transition": P3,
        "parry_stationary": pi3,
        "golden_emission_p0": p0,
        "golden_emission_p1": p1,
        "state_order": "state=64*i+g; g=sum_a 2**a b_a",
    }


def koopman_markov_audit(P, invariant=None, tolerance=1e-11):
    """Audit a finite Markov Koopman operator K=P and its L2(mu) adjoint."""
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or P.shape[0] != P.shape[1]:
        raise ValueError("P must be a square transition matrix")
    if not np.all(np.isfinite(P)) or np.min(P) < -tolerance:
        raise ValueError("P must be finite and nonnegative")
    row_res = float(np.max(np.abs(P.sum(axis=1) - 1.0)))
    if row_res > tolerance:
        raise ValueError("P must be row-stochastic")
    n = P.shape[0]
    if invariant is None:
        vals, vecs = np.linalg.eig(P.T)
        j = int(np.argmin(np.abs(vals - 1.0)))
        mu = np.real(vecs[:, j])
        if mu.sum() < 0:
            mu = -mu
        if np.min(mu) < -tolerance:
            raise ValueError("could not infer a nonnegative invariant distribution")
        mu = np.maximum(mu, 0.0)
        mu /= mu.sum()
    else:
        mu = np.asarray(invariant, dtype=float)
        if mu.shape != (n,) or not np.all(np.isfinite(mu)) or np.min(mu) < 0 or mu.sum() <= 0:
            raise ValueError("invariant must be a finite nonnegative vector")
        mu = mu / mu.sum()
    stationarity = float(np.linalg.norm(P.T @ mu - mu))
    if np.any(mu <= 0):
        raise ValueError("strictly positive invariant weights are required for the L2(mu) adjoint")
    # K*=D_mu^{-1} P.T D_mu under <f,g>_mu = f* D_mu g.
    adjoint = (P.T * mu[None, :]) / mu[:, None]
    adjoint_res = float(np.linalg.norm(
        np.diag(mu) @ P - adjoint.T @ np.diag(mu)
    ))
    eig = np.linalg.eigvals(P)
    return {
        "dimension": int(n),
        "row_stochastic_residual": row_res,
        "invariant_normalization_residual": float(abs(mu.sum() - 1.0)),
        "stationarity_residual": stationarity,
        "weighted_adjoint_identity_residual": adjoint_res,
        "constant_observable_residual": float(np.linalg.norm(P @ np.ones(n) - np.ones(n))),
        "spectral_radius": float(np.max(np.abs(eig))),
        "eigenvalues": [[float(z.real), float(z.imag)] for z in eig],
        "invariant_distribution": mu.tolist(),
        "koopman_convention": "(K f)(i)=sum_j P[i,j] f(j)",
        "adjoint_convention": "K*=D_mu^{-1} P.T D_mu",
    }


def edmd_koopman_audit(features_now, features_next, sample_weights=None, ridge=0.0):
    """Extended DMD: fit Phi_next ~= Phi_now @ K by weighted least squares.

    Rows are snapshots, columns are observables. The returned K is a finite
    Galerkin/EDMD approximation, not automatically an exact Koopman invariant
    subspace. A nonzero residual quantifies failure of closure.
    """
    X = np.asarray(features_now, dtype=float)
    Y = np.asarray(features_next, dtype=float)
    if X.ndim != 2 or Y.shape != X.shape or X.shape[0] < 1:
        raise ValueError("feature snapshots must be equally shaped 2D arrays")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError("feature snapshots must be finite")
    if ridge < 0:
        raise ValueError("ridge must be nonnegative")
    n = X.shape[0]
    if sample_weights is None:
        w = np.full(n, 1.0 / n)
    else:
        w = np.asarray(sample_weights, dtype=float)
        if w.shape != (n,) or not np.all(np.isfinite(w)) or np.min(w) < 0 or w.sum() <= 0:
            raise ValueError("sample_weights must be finite, nonnegative and nonzero")
        w = w / w.sum()
    sw = np.sqrt(w)[:, None]
    Xw, Yw = X * sw, Y * sw
    gram = Xw.T @ Xw
    rhs = Xw.T @ Yw
    reg = ridge * np.eye(gram.shape[0])
    K = np.linalg.solve(gram + reg, rhs) if ridge > 0 else np.linalg.lstsq(Xw, Yw, rcond=None)[0]
    residual = (X @ K - Y) * sw
    singular = np.linalg.svd(Xw, compute_uv=False)
    rank = int(np.linalg.matrix_rank(Xw))
    cond = float(singular[0] / singular[-1]) if singular.size and singular[-1] > 0 else None
    eig = np.linalg.eigvals(K)
    return {
        "koopman_matrix": K.tolist(),
        "rank": rank,
        "feature_dimension": int(X.shape[1]),
        "snapshot_count": int(n),
        "singular_values": singular.tolist(),
        "condition_number": cond,
        "weighted_relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(Yw), np.finfo(float).tiny)),
        "weighted_rms_residual": float(np.sqrt(np.sum(residual**2))),
        "eigenvalues": [[float(z.real), float(z.imag)] for z in eig],
        "ridge": float(ridge),
        "interpretation": "finite EDMD/Galerkin approximation; exact Koopman closure is not assumed",
    }


def barycentric_koopman_features(positions, weights, delta_k, delta_omega=0.0,
                                  times=None, phase0=0.0):
    """Build observables from barycenter, covariance and Moire characteristic.

    positions: (snapshots, particles, dimension); weights: (snapshots, particles).
    Returns a feature matrix suitable for EDMD and a named feature layout.
    """
    X = np.asarray(positions, dtype=float)
    W = np.asarray(weights, dtype=float)
    dk = np.asarray(delta_k, dtype=float)
    if X.ndim != 3 or W.shape != X.shape[:2] or dk.shape != (X.shape[2],):
        raise ValueError("incompatible positions, weights or delta_k shapes")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(W)) or not np.all(np.isfinite(dk)):
        raise ValueError("inputs must be finite")
    if np.any(W < 0) or np.any(W.sum(axis=1) <= 0):
        raise ValueError("each snapshot needs nonnegative weights with positive sum")
    W = W / W.sum(axis=1, keepdims=True)
    t = np.arange(X.shape[0], dtype=float) if times is None else np.asarray(times, dtype=float)
    if t.shape != (X.shape[0],) or not np.all(np.isfinite(t)):
        raise ValueError("times must be finite and match the snapshot count")
    rows, names = [], None
    for s in range(X.shape[0]):
        xb = W[s] @ X[s]
        Y = X[s] - xb
        cov = (Y * W[s, :, None]).T @ Y
        q = Y @ dk
        char = np.sum(W[s] * np.exp(-1j * q))
        theta = float(delta_omega) * t[s] - float(dk @ xb) + float(phase0)
        cov_values = [cov[i, j] for i in range(X.shape[2]) for j in range(i, X.shape[2])]
        row = [1.0, *xb.tolist(), *cov_values,
               float(np.cos(theta)), float(np.sin(theta)),
               float(char.real), float(char.imag),
               float(1.0 - np.real(np.exp(1j * theta) * char))]
        rows.append(row)
    d = X.shape[2]
    names = ["1"] + [f"barycenter_{i}" for i in range(d)]
    names += [f"cov_{i}_{j}" for i in range(d) for j in range(i, d)]
    names += ["cos_phase", "sin_phase", "moire_char_real", "moire_char_imag", "moire_average_exact"]
    return np.asarray(rows), names


def constant_frequency_koopman_generator(omega):
    """Exact generator on observables (cos(theta), sin(theta)), theta_dot=omega."""
    omega = float(omega)
    if not np.isfinite(omega):
        raise ValueError("omega must be finite")
    return np.array([[0.0, omega], [-omega, 0.0]])


def koopman_semigroup_residual(generator, dt, steps=1):
    """Check exp(dt*L) composition against exp(steps*dt*L) by eigendecomposition."""
    L = np.asarray(generator, dtype=float)
    if L.ndim != 2 or L.shape[0] != L.shape[1] or not np.all(np.isfinite(L)):
        raise ValueError("generator must be a finite square matrix")
    if dt < 0 or steps < 1:
        raise ValueError("dt must be nonnegative and steps positive")
    vals, vecs = np.linalg.eig(L)
    inv = np.linalg.inv(vecs)
    E = vecs @ np.diag(np.exp(vals * dt)) @ inv
    En = vecs @ np.diag(np.exp(vals * dt * steps)) @ inv
    return float(np.linalg.norm(np.linalg.matrix_power(E, steps) - En))


def build_q6_rgs_geometry_lookup(X, q6_states):
    """Build the exact finite Q6->RGS lookup from archived 64-site coordinates.

    The canonical Q6 integer is g=sum_a 2**a b_a. The archive stores rows in
    lexicographic bit order, so this routine explicitly reindexes rather than
    assuming row number equals g. This is a lookup chart, not an affine map.
    """
    X = np.asarray(X, dtype=float)
    states = np.asarray(q6_states, dtype=int)
    if X.shape != (64, 3) or states.shape != (64, 6):
        raise ValueError("X must be (64,3) and q6_states must be (64,6)")
    if not np.all(np.isfinite(X)) or not np.all((states == 0) | (states == 1)):
        raise ValueError("coordinates must be finite and states binary")
    ids = states @ (2 ** np.arange(6, dtype=int))
    if not np.array_equal(np.sort(ids), np.arange(64)):
        raise ValueError("q6_states must enumerate every canonical Q6 state once")
    lookup = np.empty((64, 3), dtype=float)
    lookup[ids] = X
    # Affine fit is diagnostic only: a poor fit must not replace the lookup.
    B = states.astype(float)
    Bc, Xc = B - B.mean(axis=0), X - X.mean(axis=0)
    coef = np.linalg.lstsq(Bc, Xc, rcond=None)[0]
    affine_residual = Xc - Bc @ coef
    return {
        "positions_by_q6_integer": lookup,
        "canonical_state_ids": ids.tolist(),
        "centroid": X.mean(axis=0).tolist(),
        "radius_squared_sum": float(np.sum(X * X)),
        "affine_projection_diagnostic": coef.T.tolist(),
        "affine_relative_residual": float(
            np.linalg.norm(affine_residual) /
            max(np.linalg.norm(Xc), np.finfo(float).tiny)
        ),
        "geometry_map": "exact 64-state lookup; no affine assumption",
        "bit_order": "g=sum_a 2**a*b_a; archived state rows explicitly reindexed",
    }
