#!/usr/bin/env python3
"""RGS-PSU: diagnosi dell'incongruenza di base e costruzione esplicita di una coppia compatibile.
Uso: python RGS_PSU_resolve_basis.py  (nella cartella del bundle v0.1)  -> scrive RGS_PSU_basis_resolution.json/.npz
"""
import numpy as np, json
from scipy.linalg import eigh, null_space
d = np.load('RGS_PSU_all_matrices.npz', allow_pickle=False)
K0, PG, R, A, Pi = (d[k] for k in ['K0', 'P_G', 'R_geo_tangent', 'A_geometry_orthonormal', 'Pi_AT'])
K0 = (K0 + K0.T) / 2; PG = (PG + PG.T) / 2; n = 188
nrm = lambda M: float(np.linalg.norm(M))
rep = {}

# ---- 1. DIAGNOSI -------------------------------------------------------------
rng = np.random.default_rng(1); rv = []
for _ in range(500):
    O, _ = np.linalg.qr(rng.standard_normal((n, 64))); P = O @ O.T; rv.append(nrm(R @ P - P @ R))
obs = nrm(R @ PG - PG @ R)
rep['commutator_observed'] = obs
rep['commutator_random_mean_std'] = [float(np.mean(rv)), float(np.std(rv))]
rep['z_score_vs_random'] = float((obs - np.mean(rv)) / np.std(rv))
w, E = eigh(K0)
rep['K0_nonzero_spectrum_max_dev_from_0.5+0.1j'] = float(np.abs(w[64:] - (0.5 + 0.1 * np.arange(124))).max())
rep['K0_nonzero_min_gap'] = float(np.diff(w[64:]).min())
rep['R_signature'] = [int((np.linalg.eigvalsh(R) > 0).sum()), int((np.linalg.eigvalsh(R) < 0).sum())]
# estensione "identita' nel complemento" R_K = U Pi U^T + (I-P_G): segnatura
B = PG @ A; s, V = eigh(B.T @ B); U = B @ V @ np.diag(s ** -.5) @ V.T
RK = U @ Pi @ U.T + np.eye(n) - PG
rep['R_K_signature'] = [int((np.linalg.eigvalsh(RK) > 0).sum()), int((np.linalg.eigvalsh(RK) < 0).sum())]
rep['trace_R_geo_vs_R_K'] = [float(np.trace(R)), float(np.trace(RK))]

# ---- 2. COSTRUZIONE DI W ------------------------------------------------------
Z, Ev, ev = E[:, :64], E[:, 64:], w[64:]                       # kernel e modi non nulli di K0
wp, Vp = eigh((Pi + Pi.T) / 2); Vp_, Vm_ = Vp[:, wp > 0], Vp[:, wp < 0]   # +-1 di Pi (32+32)
Cc = null_space(A.T)                                            # complemento di Im A (124)
wc, Vc = eigh((Cc.T @ R @ Cc + (Cc.T @ R @ Cc).T) / 2)
Cp, Cm = Cc @ Vc[:, wc > 0], Cc @ Vc[:, wc < 0]                 # 63 e 61
# involuzione piu' vicina a R che commuta con K0 (segni nella base propria di K0), vincolo di segnatura
m = np.einsum('ij,ij->j', Ev, R @ Ev); order = np.argsort(-m)
sp, sm = np.sort(order[:Cp.shape[1]]), np.sort(order[Cp.shape[1]:])
Ep, Em = Ev[:, sp], Ev[:, sm]
Rprime = U @ Pi @ U.T + Ep @ Ep.T - Em @ Em.T
W = U @ A.T + Ep @ Cp.T + Em @ Cm.T
Kg = W.T @ K0 @ W                                               # K0 trasportato nella base tangente
rep['W_orthogonality'] = nrm(W.T @ W - np.eye(n))
rep['W R W^T - R_prime'] = nrm(W @ R @ W.T - Rprime)
rep['[R_geo, K_geo]'] = nrm(R @ Kg - Kg @ R)
rep['W^T P_G W - A A^T'] = nrm(W.T @ PG @ W - A @ A.T)
rep['[R_geo, W^T P_G W]'] = nrm(R @ (W.T @ PG @ W) - (W.T @ PG @ W) @ R)
rep['K_geo spectrum = K0 spectrum'] = float(np.abs(np.linalg.eigvalsh(Kg) - w).max())
rep['||R - R_prime|| (distanza dall identita di base)'] = nrm(R - Rprime)
rep['||R - R_K||'] = nrm(R - RK)
# equivarianza trasporto polare nella base trasportata
Bn = (W.T @ PG @ W) @ A; s2, V2 = eigh(Bn.T @ Bn); Un = Bn @ V2 @ np.diag(s2 ** -.5) @ V2.T
rep['R_geo U_new - U_new Pi'] = nrm(R @ Un - Un @ Pi)
json.dump(rep, open('RGS_PSU_basis_resolution.json', 'w'), indent=2)
np.savez('RGS_PSU_basis_resolution.npz', W=W, K_geo=Kg, R_prime=Rprime)
for k, v in rep.items(): print(k, v)
