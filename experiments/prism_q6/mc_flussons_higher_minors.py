#!/usr/bin/env python3
"""
1. Monte Carlo leggero dell'azione flussoni (plaquette Q6)
2. Estensione mappa sorted a minori di ordine superiore (k=3, Gr(3,6))
"""

import numpy as np
from itertools import combinations
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("MONTE CARLO FLUSSONI + MINORI DI ORDINE SUPERIORE")
print("=" * 70)
t0 = time.time()

# ============================================================
# 1. MONTE CARLO AZIONE FLUSSONI
# ============================================================
print("\n[1] Monte Carlo leggero – azione flussoni su 240 plaquette")

N_PLAQ = 240
Z_EFF = 4            # grado medio duale (intero)
N_SWEEP = 400        # sweep Metropolis
N_THERM = 100
BETA = 2.5
KAPPA = 0.3          # punto vicino a m/√σ ≈ 3

# Grafo duale semplificato: ogni plaquette collegata a Z_EFF vicine
# (lista di adiacenza ciclica per velocità)
adj = [[] for _ in range(N_PLAQ)]
for i in range(N_PLAQ):
    for d in range(Z_EFF):
        j = (i + 1 + d) % N_PLAQ
        if j not in adj[i]:
            adj[i].append(j)
            adj[j].append(i)

def action(theta, beta, kappa):
    S_w = beta * np.sum(1.0 - np.cos(theta))
    S_k = 0.0
    for i in range(N_PLAQ):
        for j in adj[i]:
            if j > i:
                S_k += 0.5 * kappa * (theta[i] - theta[j])**2
    return S_w + S_k

def local_delta_S(theta, i, new_th, beta, kappa):
    # variazione solo dai termini che coinvolgono i
    dS = beta * ((1-np.cos(new_th)) - (1-np.cos(theta[i])))
    for j in adj[i]:
        dS += 0.5 * kappa * ((new_th - theta[j])**2 - (theta[i] - theta[j])**2)
    return dS

# Metropolis
theta = np.random.uniform(-np.pi, np.pi, N_PLAQ)
S = action(theta, BETA, KAPPA)
S_hist = []
accept = 0
total = 0

for sweep in range(N_SWEEP):
    for i in range(N_PLAQ):
        prop = theta[i] + np.random.normal(0, 0.4)
        # wrap
        prop = (prop + np.pi) % (2*np.pi) - np.pi
        dS = local_delta_S(theta, i, prop, BETA, KAPPA)
        if dS < 0 or np.random.rand() < np.exp(-dS):
            theta[i] = prop
            S += dS
            accept += 1
        total += 1
    if sweep >= N_THERM:
        S_hist.append(S)

S_hist = np.array(S_hist)
acc_rate = accept / total

# Osservabili
# massa proxy: da correlazione a 1 passo sul duale
# C(1) = <θ_i θ_{i+1}> / <θ²>
corr1 = 0.0
var_th = 0.0
n_meas = 0
# rifai qualche misura
for _ in range(50):
    for i in range(N_PLAQ):
        prop = theta[i] + np.random.normal(0, 0.3)
        prop = (prop + np.pi) % (2*np.pi) - np.pi
        dS = local_delta_S(theta, i, prop, BETA, KAPPA)
        if dS < 0 or np.random.rand() < np.exp(-dS):
            theta[i] = prop
    # misura
    corr1 += np.mean([theta[i]*theta[(i+1)%N_PLAQ] for i in range(N_PLAQ)])
    var_th += np.mean(theta**2)
    n_meas += 1
corr1 /= n_meas
var_th /= n_meas
xi_proxy = -1.0 / np.log(abs(corr1)/(var_th+1e-12) + 1e-12) if var_th > 0 else 0.0
m_proxy = 1.0 / (xi_proxy + 1e-12)

print(f"  β = {BETA},  κ = {KAPPA}")
print(f"  Sweep totali        : {N_SWEEP} (therm {N_THERM})")
print(f"  Acceptance rate     : {acc_rate:.3f}")
print(f"  <S>                 : {S_hist.mean():.2f} ± {S_hist.std():.2f}")
print(f"  <θ²>                : {var_th:.4f}")
print(f"  C(1)/<θ²>           : {corr1/(var_th+1e-12):.4f}")
print(f"  ξ_proxy (duale)     : {xi_proxy:.3f}")
print(f"  m_proxy = 1/ξ       : {m_proxy:.3f}")
print(f"  m²_tree-level       : {2*BETA + Z_EFF*KAPPA:.3f}")
print(f"  m_tree              : {np.sqrt(2*BETA + Z_EFF*KAPPA):.3f}")

# ============================================================
# 2. ESTENSIONE SORTED A MINORI DI ORDINE SUPERIORE (k=3)
# ============================================================
print("\n[2] Estensione mappa sorted a Gr(3,6) e Gr(2,6)")

N_S = 8192
N_EFF = 6

idx = np.random.choice(69120, N_S, replace=False)
e = idx%5; tmp=idx//5; c=tmp%6; tmp//=6; b=tmp%6; tmp//=6; a=tmp%6; i64=tmp//6
coords = np.stack([i64,a,b,c,e],1).astype(float)
u = coords / np.array([64.,6,6,6,5])
u6 = np.zeros((N_S, N_EFF))
u6[:,0] = u[:,0]
u6[:,1] = (u[:,1]+0.1*u[:,0])%1
u6[:,2] = (u[:,2]+0.1*u[:,1])%1
u6[:,3] = (u[:,3]+0.1*u[:,2])%1
u6[:,4] = (u[:,4]+0.1*u[:,3])%1
u6[:,5] = (0.5*u[:,0]+0.3*u[:,4]+0.2*u[:,1])%1
u_sorted = np.sort(u6, axis=1)

def soft_pos_k(C, k):
    """Soft-positivity per minori massimali di ordine k."""
    minors_idx = list(combinations(range(N_EFF), k))
    n_min = len(minors_idx)
    # per k=2 det 2x2; per k=3 det 3x3
    pos_count = np.zeros(C.shape[0])
    for idxs in minors_idx:
        sub = C[:, :, list(idxs)]          # (N, k, k)
        if k == 2:
            det = sub[:,0,0]*sub[:,1,1] - sub[:,0,1]*sub[:,1,0]
        elif k == 3:
            det = (sub[:,0,0]*(sub[:,1,1]*sub[:,2,2]-sub[:,1,2]*sub[:,2,1])
                  -sub[:,0,1]*(sub[:,1,0]*sub[:,2,2]-sub[:,1,2]*sub[:,2,0])
                  +sub[:,0,2]*(sub[:,1,0]*sub[:,2,1]-sub[:,1,1]*sub[:,2,0]))
        else:
            det = np.array([np.linalg.det(sub[s]) for s in range(C.shape[0])])
        pos_count += (det >= -1e-12).astype(float)
    return pos_count / n_min

# k=2 (già fatto, verifica)
C2 = np.zeros((N_S, 2, N_EFF))
C2[:,0,:] = np.cumsum(u_sorted + 0.05, axis=1)
C2[:,1,:] = np.cumsum((u_sorted + 0.05)**2, axis=1)
sp2 = soft_pos_k(C2, 2)

# k=3
C3 = np.zeros((N_S, 3, N_EFF))
C3[:,0,:] = np.cumsum(u_sorted + 0.05, axis=1)
C3[:,1,:] = np.cumsum((u_sorted + 0.05)**2, axis=1)
C3[:,2,:] = np.cumsum((u_sorted + 0.05)**3, axis=1)
sp3 = soft_pos_k(C3, 3)

print(f"  Stati campionati           : {N_S}")
print(f"  Soft-pos Gr(2,6) sorted    : {sp2.mean():.4f}  (frac>0.95: {np.mean(sp2>0.95):.4f})")
print(f"  Soft-pos Gr(3,6) sorted    : {sp3.mean():.4f}  (frac>0.95: {np.mean(sp3>0.95):.4f})")

if sp3.mean() >= 0.95:
    print("  → Estensione a k=3 RIUSCITA (positività quasi totale)")
elif sp3.mean() >= 0.80:
    print("  → Estensione a k=3 parziale; da raffinare")
else:
    print("  → Estensione a k=3 debole")

# ============================================================
# RIEPILOGO
# ============================================================
print("\n" + "=" * 70)
print("RIEPILOGO")
print("=" * 70)
print(f"  MC flussoni: <S>={S_hist.mean():.1f}, m_proxy={m_proxy:.3f}, m_tree={np.sqrt(2*BETA+Z_EFF*KAPPA):.3f}")
print(f"  Soft-pos k=2 (sorted): {sp2.mean():.4f}")
print(f"  Soft-pos k=3 (sorted): {sp3.mean():.4f}")
print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/mc_flussons_higher_minors_results.npz",
         S_mean=S_hist.mean(), S_std=S_hist.std(), m_proxy=m_proxy,
         sp2=sp2.mean(), sp3=sp3.mean(), acc_rate=acc_rate)
print("\nSalvato mc_flussons_higher_minors_results.npz")
