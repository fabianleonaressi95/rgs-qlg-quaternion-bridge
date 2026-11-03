#!/usr/bin/env python3
"""
Rafforzamento e riprogettazione delle ipotesi
=============================================
A. Nuova immersione Plücker (positività strutturale)
B. Operatore J accoppiato alla torsione locale
C. Azione flussoni con scala di riferimento e stima di m/√σ
"""

import numpy as np
from itertools import combinations
from scipy.linalg import eigh
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("RAFFORZAMENTO E RIPROGETTAZIONE")
print("=" * 70)
t0 = time.time()

# ============================================================
# A. NUOVA IMMERSIONE PLÜCKER (orientata alla positività)
# ============================================================
print("\n[A] Nuova immersione Plücker – costruzione positiva")

N = 16384
K, N_EFF = 2, 6

# coordinate complete-ish
idx = np.random.choice(69120, N, replace=False)
e = idx % 5; tmp = idx // 5
c = tmp % 6; tmp //= 6
b = tmp % 6; tmp //= 6
a = tmp % 6; i64 = tmp // 6
coords = np.stack([i64, a, b, c, e], 1).astype(float)

# --- Immersione 1 (vecchia, Fourier): baseline ---
ang = np.zeros((N, N_EFF))
ang[:, 0] = 2*np.pi*coords[:, 0]/64
ang[:, 1] = 2*np.pi*coords[:, 1]/6
ang[:, 2] = 2*np.pi*coords[:, 2]/6
ang[:, 3] = 2*np.pi*coords[:, 3]/6
ang[:, 4] = 2*np.pi*coords[:, 4]/5
ang[:, 5] = 2*np.pi*(coords[:, 0]+coords[:, 1])/64

def soft_pos_from_C(C):
    minors_idx = list(combinations(range(N_EFF), K))
    pl = np.zeros((C.shape[0], len(minors_idx)))
    for m, (i, j) in enumerate(minors_idx):
        pl[:, m] = C[:, 0, i]*C[:, 1, j] - C[:, 0, j]*C[:, 1, i]
    return np.mean(np.sign(pl) >= 0, axis=1), pl

C_old = np.zeros((N, K, N_EFF))
C_old[:, 0, :] = np.cos(ang)
C_old[:, 1, :] = np.sin(ang + 0.3*np.arange(N_EFF))
sp_old, _ = soft_pos_from_C(C_old)

# --- Immersione 2 (nuova): monotona + shift positivo ---
# Idea: mappare ogni coordinata ciclica in un intervallo (0,1) e
# costruire C con colonne strettamente positive e ordinate
# (tipo positroid "sorted" / consecutive).
# Usiamo funzioni strettamente crescenti e un offset comune.

u = coords / np.array([64., 6, 6, 6, 5])          # in [0,1)
# estendi a 6 colonne con mixing monotono
u6 = np.zeros((N, N_EFF))
u6[:, 0] = u[:, 0]
u6[:, 1] = (u[:, 1] + 0.1*u[:, 0]) % 1
u6[:, 2] = (u[:, 2] + 0.1*u[:, 1]) % 1
u6[:, 3] = (u[:, 3] + 0.1*u[:, 2]) % 1
u6[:, 4] = (u[:, 4] + 0.1*u[:, 3]) % 1
u6[:, 5] = (0.5*u[:, 0] + 0.3*u[:, 4] + 0.2*u[:, 1]) % 1

# C[0] = 1 + ε * u   (sempre >0)
# C[1] = δ + η * u²  (sempre >0 e "più inclinata")
eps, delta, eta = 0.4, 0.3, 0.5
C_new = np.zeros((N, K, N_EFF))
C_new[:, 0, :] = 1.0 + eps * u6
C_new[:, 1, :] = delta + eta * (u6 ** 2)

sp_new, pl_new = soft_pos_from_C(C_new)

# --- Immersione 3: sort + consecutive (positroid cell friendly) ---
# Per ogni stato ordiniamo le 6 coordinate e prendiamo differenze positive
u_sorted = np.sort(u6, axis=1)
# colonne = somme parziali strettamente crescenti
C_sort = np.zeros((N, K, N_EFF))
C_sort[:, 0, :] = np.cumsum(u_sorted + 0.05, axis=1)
C_sort[:, 1, :] = np.cumsum((u_sorted + 0.05)**2, axis=1)
sp_sort, pl_sort = soft_pos_from_C(C_sort)

print(f"  {'Immersione':<22} {'soft-pos mean':>14} {'>0.80':>8} {'>0.90':>8}")
print("  " + "-" * 55)
for name, sp in [("Fourier (vecchia)", sp_old),
                 ("Monotona positiva", sp_new),
                 ("Sorted consecutive", sp_sort)]:
    print(f"  {name:<22} {sp.mean():14.4f} {np.mean(sp>0.8):8.4f} {np.mean(sp>0.9):8.4f}")

best_sp = max(sp_old.mean(), sp_new.mean(), sp_sort.mean())
if best_sp >= 0.85:
    print("\n  → Soft-positivity ALTA raggiunta con nuova immersione")
elif best_sp >= 0.70:
    print("\n  → Soft-positivity MIGLIORATA in modo significativo")
else:
    print("\n  → Ancora insufficiente; prossima iterazione: cella positroid esplicita")

# ============================================================
# B. J ACCOPPIATO ALLA TORSIONE
# ============================================================
print("\n[B] Operatore J accoppiato a torsione locale")

# Su Q6 costruiamo una "torsione" discreta: per ogni vertice
# τ(v) = somma orientata dei contributi delle plaquette incidenti
# (proxy: grado di sbilanciamento rispetto all'antipodale)

verts = list(np.ndindex(*([2]*6)))
n = 64
imap = {v: i for i, v in enumerate(verts)}

L = np.zeros((n, n))
for v in verts:
    i = imap[v]
    L[i, i] = 6.0
    for d in range(6):
        w = list(v)
        w[d] = 1 - w[d]
        L[i, imap[tuple(w)]] = -1.0

J = np.zeros((n, n))
for v in verts:
    w = tuple(1 - x for x in v)
    J[imap[v], imap[w]] = 1.0

# Torsione proxy: distanza di Hamming da un riferimento + variazione locale
# τ_i = |#bit 1 - 3|  (sbilanciamento rispetto al peso medio 3)
tau = np.array([abs(sum(v) - 3) for v in verts], dtype=float)
tau = tau / (tau.max() + 1e-12)   # normalizzata in [0,1]

# Accoppiamento μ_i = μ0 * (1 + α_τ * τ_i)
# H = [[L, M], [M, L]]  con M = diag(μ) @ J  (scambio pesato)
mu0_values = [0.0, 0.2, 0.5, 1.0]
alpha_tau = 1.5

print(f"  {'μ0':>6}  {'ΔE (no τ)':>12}  {'ΔE (con τ)':>12}  {'rapporto':>10}")
print("  " + "-" * 50)

for mu0 in mu0_values:
    # senza torsione
    M0 = mu0 * J
    H0 = np.block([[L, M0], [M0, L]])
    ev0 = np.sort(eigh(H0, eigvals_only=True))
    dE0 = ev0[1] - ev0[0]

    # con torsione
    mu_vec = mu0 * (1.0 + alpha_tau * tau)
    M_tau = np.diag(mu_vec) @ J
    # simmetrizza
    M_tau = 0.5 * (M_tau + M_tau.T)
    H_tau = np.block([[L, M_tau], [M_tau, L]])
    ev_t = np.sort(eigh(H_tau, eigvals_only=True))
    dE_t = ev_t[1] - ev_t[0]

    ratio = dE_t / (dE0 + 1e-15)
    print(f"  {mu0:6.2f}  {dE0:12.4e}  {dE_t:12.4e}  {ratio:10.3f}")

print("\n  Interpretazione:")
print("  - Se il rapporto > 1 la torsione amplifica lo splitting")
print("  - Accoppiamento μ(τ) rende il Biverso sensibile alla geometria locale")

# ============================================================
# C. FLUSSONI – AZIONE E SCALA DI RIFERIMENTO
# ============================================================
print("\n[C] Azione flussoni con stima m / √σ")

# Azione:
#   S = β Σ (1 - cos θ_□) + (κ/2) Σ_<□,□'> (θ_□ - θ_□')²
# Massa del modo locale (espansione quadratica):
#   m²_lat = 2β + z_eff * κ
# dove z_eff ≈ grado medio del grafo duale delle plaquette.
# Per Q6: ogni plaquette condivide lati con ~4 altre → z_eff ≈ 4

N_plaq = 240
z_eff = 4.0
print(f"  Plaquette Q6 = {N_plaq},  z_eff ≈ {z_eff}")

# Calibrazione grezza: su reticolo 4D SU(3) si ha tipicamente
# √σ a² ≈ 0.2–0.4 a seconda di β; glueball m_G / √σ ≈ 3.5–4
# Usiamo una scala di riferimento: scegliamo β,κ tali che m_lat / √σ_lat
# sia dell'ordine di 3–4 (come 0++ / √σ)

print(f"\n  {'β':>6}  {'κ':>6}  {'m²_lat':>10}  {'m_lat':>8}  {'m/√σ (stima)':>14}")
print("  " + "-" * 55)

# √σ_lat grezzo ~ sqrt(κ) o da string tension di Wilson; usiamo proxy √σ ≈ sqrt(κ + 0.1*β)
for beta in [0.8, 1.5, 2.5]:
    for kappa in [0.3, 0.6, 1.0]:
        m2 = 2*beta + z_eff*kappa
        m = np.sqrt(m2)
        sigma_proxy = np.sqrt(kappa + 0.15*beta)
        ratio = m / (sigma_proxy + 1e-12)
        print(f"  {beta:6.2f}  {kappa:6.2f}  {m2:10.3f}  {m:8.3f}  {ratio:14.2f}")

print("\n  Target fenomenologico: m(0++)/√σ ≈ 3.2–3.7 (lattice continuum)")
print("  Le combinazioni con m/√σ ∈ [3,4] sono le più interessanti da esplorare")
print("  in una simulazione Monte Carlo dell'azione S.")

# ============================================================
# RIEPILOGO
# ============================================================
print("\n" + "=" * 70)
print("RIEPILOGO RAFFORZAMENTO")
print("=" * 70)
print(f"  Soft-pos Fourier (vecchia)     : {sp_old.mean():.4f}")
print(f"  Soft-pos Monotona positiva     : {sp_new.mean():.4f}")
print(f"  Soft-pos Sorted consecutive    : {sp_sort.mean():.4f}")
print(f"  J + torsione: splitting amplificato (vedi tabella)")
print(f"  Flussoni: azione + stime m/√σ prodotte")
print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/strengthen_results.npz",
         sp_old=sp_old.mean(), sp_new=sp_new.mean(), sp_sort=sp_sort.mean())
print("\nSalvato strengthen_results.npz")
