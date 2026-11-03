#!/usr/bin/env python3
"""
Blocco B — Scaling verso il continuo
====================================
- Sequenza di reticoli effettivi di taglia crescente (n_p plaquette)
- Traiettoria di scaling g(a) tipo asintotico libero (scheletro)
- Δ(g, L), Δ_phys = Δ / a
- Criterio: Δ_phys resta limitato inferiormente mentre a → 0?
"""

import numpy as np
from scipy.linalg import eigh
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("BLOCCO B — Scaling Δ_phys = Δ/a")
print("=" * 70)
t0 = time.time()

# ============================================================
# Modello a plaquette scalabile (stesso di Blocco A, n_p variabile)
# ============================================================
def make_adj_cycle(n_p, degree=4):
    """Grafo duale regolare approssimato (ciclo con degree/2 vicini per lato)."""
    adj = [[] for _ in range(n_p)]
    half = max(1, degree // 2)
    for i in range(n_p):
        for d in range(1, half+1):
            j = (i + d) % n_p
            if j not in adj[i]:
                adj[i].append(j)
                adj[j].append(i)
    return adj

def H_eff(g, n_p, adj, lam_factor=0.25):
    """H su 2^{n_p} (ok fino a n_p≈12–14)."""
    dim = 1 << n_p
    dE = 1.5 * g**2
    t = 2.0 / (g**2 + 1e-15)
    lam = lam_factor * g**2
    # Costruzione densa (per n_p piccoli)
    H = np.zeros((dim, dim))
    for s in range(dim):
        bits = [(s >> k) & 1 for k in range(n_p)]
        e = dE * sum(bits)
        for i in range(n_p):
            for j in adj[i]:
                if j > i:
                    e += lam * bits[i] * bits[j]
        H[s, s] = e
        for i in range(n_p):
            H[s, s ^ (1 << i)] -= t
    return H

def gap_eff(g, n_p, adj):
    H = H_eff(g, n_p, adj)
    ev = np.sort(eigh(H, eigvals_only=True))
    return ev[1] - ev[0], ev[0]

# ============================================================
# Sequenza di “reticoli”: n_p = 4, 6, 8, 10  (proxy di L crescente)
# a = 1/L_proxy con L_proxy = sqrt(n_p)  (2D-like) oppure n_p^{1/4} (4D-like)
# ============================================================
print("\n[1] Gap Δ(g, n_p) su reticoli effettivi crescenti")

n_p_list = [4, 6, 8, 10]
g_list = [1.0, 1.5, 2.0, 2.5, 3.0]

print(f"  {'n_p':>4}  {'g':>6}  {'Δ':>10}  {'E0':>10}")
print("  " + "-" * 36)

results = {}  # (n_p, g) -> delta
for n_p in n_p_list:
    adj = make_adj_cycle(n_p, degree=4)
    for g in g_list:
        dlt, E0 = gap_eff(g, n_p, adj)
        results[(n_p, g)] = dlt
        print(f"  {n_p:4d}  {g:6.2f}  {dlt:10.4f}  {E0:10.4f}")

# ============================================================
# Traiettoria di scaling (scheletro asintotico libero SU(2))
# β = 4/g²,  aΛ ≈ exp(-2π² β / β0)  con β0 = 22/3  (1-loop SU(2))
# Invertiamo: fissiamo una scala fisica e variamo a
# ============================================================
print("\n[2] Traiettoria di scaling e Δ_phys = Δ/a")

# Proxy: a(n_p) = 1 / n_p^{1/4}   (come se volume ~ n_p in 4D, L = n_p^{1/4})
# g(a) lungo una traiettoria in cui β aumenta (g diminuisce) al diminuire di a
# Forma schematica: g²(a) = g0² / (1 + b0 log(1/(a Λ)))
# con parametri di esempio

def g_of_a(a, g0=2.5, b0=0.3, Lambda=1.0):
    """Traiettoria schematica: g diminuisce quando a diminuisce."""
    return g0 / np.sqrt(1.0 + b0 * np.log(1.0/(a*Lambda) + 1.0))

print(f"  {'n_p':>4}  {'a':>8}  {'g(a)':>8}  {'Δ':>10}  {'Δ/a':>12}")
print("  " + "-" * 50)

a_list = []
dphys_list = []
for n_p in n_p_list:
    a = n_p ** (-0.25)
    g = g_of_a(a)
    # interpola Δ sul g più vicino nella griglia, o ricalcola
    # ricalcola esatto
    adj = make_adj_cycle(n_p, degree=4)
    dlt, _ = gap_eff(g, n_p, adj)
    dphys = dlt / a
    a_list.append(a)
    dphys_list.append(dphys)
    print(f"  {n_p:4d}  {a:8.4f}  {g:8.4f}  {dlt:10.4f}  {dphys:12.4f}")

a_list = np.array(a_list)
dphys_list = np.array(dphys_list)

# Trend: se Δ_phys cresce o stabilizza mentre a→0 → supporto;
# se collassa → ipotesi di gap continuo indebolita su questo modello.
print("\n[3] Criterio sul continuo (questo modello)")
slope = (dphys_list[-1] - dphys_list[0]) / (a_list[-1] - a_list[0] + 1e-15)
print(f"  a: {a_list[0]:.4f} → {a_list[-1]:.4f}")
print(f"  Δ_phys: {dphys_list[0]:.4f} → {dphys_list[-1]:.4f}")
print(f"  d(Δ_phys)/da (grezzo) ≈ {slope:.4f}")
if dphys_list[-1] >= dphys_list[0] * 0.5 and dphys_list[-1] > 0:
    print("  → Su questa traiettoria Δ_phys non collassa: compatibile con gap")
    print("    (modello a plaquette, non ancora continuum YM)")
else:
    print("  → Δ_phys in calo significativo: segnale di chiusura del gap sul modello")

# ============================================================
# Bound α lungo la traiettoria
# ============================================================
print("\n[4] α_est lungo la traiettoria di scaling")
print(f"  {'n_p':>4}  {'g(a)':>8}  {'α_est':>10}  {'α<1?':>6}")
print("  " + "-" * 36)
for n_p, a in zip(n_p_list, a_list):
    g = g_of_a(a)
    dE = 1.5 * g**2
    HBnorm = 2.0 * n_p / (g**2 + 1e-15)
    al = HBnorm / (dE + 1e-15)
    print(f"  {n_p:4d}  {g:8.4f}  {al:10.4f}  {'sì' if al<1 else 'no':>6}")

print("""
  Nota epistemica:
  - Il modello a plaquette non è Yang-Mills continuo.
  - Scaling e Δ/a qui sono un PROTOTIPO del Blocco B sul formalismo
    che useremo su reticoli KS veri (Q3, Q4, tori, …).
  - Il criterio falsificabile resta: Δ_phys limitato inferiormente
    su una sequenza a→0 con g(a) lungo il flusso di rinormalizzazione.
""")

print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/ks_block_B_results.npz",
         n_p_list=np.array(n_p_list), g_list=np.array(g_list),
         results={str(k):v for k,v in results.items()},
         a_list=a_list, dphys_list=dphys_list)
print("\nSalvato ks_block_B_results.npz")
