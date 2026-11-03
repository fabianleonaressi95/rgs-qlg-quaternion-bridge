#!/usr/bin/env python3
"""
Ipotesi: la rinormalizzazione contiene la sezione superaurea ψ
=============================================================
ψ = radice reale di x³ - x² - 1 = 0  ≈ 1.465571231876768
Ricorrenza: K³ - K² - I = 0  ⇒ autovalori in {ψ, ω±}

Test:
1. Identità algebrica di ψ e scaling discreto a_{n+1} = a_n / ψ
2. Inserire ψ in una traiettoria di RG schematica e confrontare Δ_phys
3. Verificare se un operatore di “blocco” 3→1 (tipo K) riproduce la ricorrenza
4. Criterio falsificabile: predizione quantitativa su Δ(g) o su g(a)
"""

import numpy as np
from scipy.linalg import eigh
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("IPOTESI: RG CONTIENE LA SEZIONE SUPERAUREA ψ")
print("=" * 70)
t0 = time.time()

# ---------- Costante ψ ----------
# radice reale di x³ - x² - 1 = 0
psi = (1 + ((29+3*np.sqrt(93))/2)**(1/3) + ((29-3*np.sqrt(93))/2)**(1/3)) / 3
# verifica
res = psi**3 - psi**2 - 1
print(f"\n[0] ψ = {psi:.15f}")
print(f"    ψ³ - ψ² - 1 = {res:.3e}")

# ---------- 1. Scaling discreto a_{n+1} = a_n / ψ ----------
print("\n[1] Traiettoria di reticolo con fattore ψ")
a0 = 1.0
n_steps = 8
a_psi = [a0]
for _ in range(n_steps):
    a_psi.append(a_psi[-1] / psi)
a_psi = np.array(a_psi)
print(f"  a_n = a0 / ψ^n")
print(f"  n=0..{n_steps}: a = {np.array2string(a_psi, precision=4)}")

# Confronto con scaling 2 (blocco 2^4 in 4D) e con φ (golden)
phi = (1+np.sqrt(5))/2
print(f"  ψ = {psi:.6f}  |  φ = {phi:.6f}  |  2^{1/4} = {2**0.25:.6f}")

# ---------- 2. RG schematica con β-function e inserzione di ψ ----------
print("\n[2] Due traiettorie g(a): standard 1-loop vs ‘ψ-step’")

# SU(2) 1-loop: β(g) = -b0 g³/(16π²) con b0 = 22/3
# a dg/da = β(g)  ⇒  1/g²(a) - 1/g²(a0) = (b0/(8π²)) log(a0/a)
b0 = 22.0/3.0
def g_1loop(a, a0=1.0, g0=2.5):
    inv_g2 = 1/g0**2 + (b0/(8*np.pi**2)) * np.log(a0/(a+1e-30))
    return 1.0/np.sqrt(max(inv_g2, 1e-8))

# Traiettoria ‘ψ-step’: a ogni step discreto a → a/ψ, g si aggiorna
# con la stessa β ma Δlog = log(ψ)
def g_psi_step(n, g0=2.5):
    inv_g2 = 1/g0**2 + (b0/(8*np.pi**2)) * n * np.log(psi)
    return 1.0/np.sqrt(max(inv_g2, 1e-8))

print(f"  {'n':>3}  {'a=1/ψ^n':>10}  {'g_1loop':>10}  {'g_ψ-step':>10}")
print("  " + "-" * 40)
for n in range(n_steps+1):
    a = a_psi[n]
    print(f"  {n:3d}  {a:10.5f}  {g_1loop(a):10.5f}  {g_psi_step(n):10.5f}")

# ---------- 3. Operatore di blocco tipo Krist-All ----------
print("\n[3] Operatore di blocco e ricorrenza K³−K²−I=0")

# Costruiamo una matrice di trasferimento / blocco 3×3
# che soddisfa la ricorrenza per costruzione (companion matrix)
# K companion di x³ - x² - 1:
# [[0,0,1],[1,0,0],[0,1,1]]  o analoga
K = np.array([
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 1.0]
], dtype=float)
# Verifica K³ - K² - I
I = np.eye(3)
resK = K@K@K - K@K - I
print(f"  ||K³ - K² - I||_F = {np.linalg.norm(resK):.3e}")
evK = np.sort(np.linalg.eigvals(K))
print(f"  autovalori di K   = {evK}")
print(f"  |λ_max|           = {np.abs(evK).max():.15f}  (deve ≈ ψ)")

# Interpretazione RG: uno step di blocco dilata le lunghezze di ψ
# e scala gli operatori con potenze di ψ (dimensione anomala discreta)
print("""
  Lettura RG:
  - Uno step K: a → a/ψ
  - Un operatore locale O di dimensione Δ_O scala come ψ^{Δ_O}
  - Il gap di massa (dimensione 1) scala come Δ → ψ Δ  (in unità di reticolo)
    quindi Δ_phys = Δ_lat / a resta invariante se Δ_lat → Δ_lat / ψ quando a → a/ψ
""")

# ---------- 4. Test numerico sul modello a plaquette ----------
print("\n[4] Test sul modello a plaquette: Δ_lat(a) vs predizione ψ")

def H_eff_gap(g, n_p=6):
    dim = 1 << n_p
    dE = 1.5 * g**2
    t = 2.0 / (g**2 + 1e-15)
    lam = 0.25 * g**2
    adj = [((i+1)%n_p, (i+2)%n_p) for i in range(n_p)]
    H = np.zeros((dim, dim))
    for s in range(dim):
        bits = [(s>>k)&1 for k in range(n_p)]
        e = dE*sum(bits)
        for i in range(n_p):
            for j in adj[i]:
                if j>i: e += lam*bits[i]*bits[j]
        H[s,s] = e
        for i in range(n_p):
            H[s, s^(1<<i)] -= t
    ev = np.sort(eigh(H, eigvals_only=True))
    return ev[1]-ev[0]

# Lungo la traiettoria ψ: a_n = 1/ψ^n, g = g_1loop(a_n)
# Predizione: se RG è “ψ-covariante”, Δ_lat(n+1) ≈ Δ_lat(n)/ψ
print(f"  {'n':>3}  {'a':>8}  {'g':>8}  {'Δ_lat':>10}  {'Δ_lat·ψ':>10}  {'Δ_{n-1}':>10}  {'ratio':>8}")
print("  " + "-" * 70)
prev = None
ratios = []
for n in range(0, 6):
    a = 1.0 / (psi**n)
    g = g_1loop(a, a0=1.0, g0=2.5)
    # n_p cresce come volume ~ 1/a^4 proxy
    n_p = min(10, max(4, int(6 * (1.0/a)**0.5)))
    # per stabilità usiamo n_p fisso=6 e solo g(a)
    n_p = 6
    dlt = H_eff_gap(g, n_p)
    if prev is not None:
        ratio = prev / (dlt + 1e-30)
        ratios.append(ratio)
        print(f"  {n:3d}  {a:8.4f}  {g:8.4f}  {dlt:10.4f}  {dlt*psi:10.4f}  {prev:10.4f}  {ratio:8.4f}")
    else:
        print(f"  {n:3d}  {a:8.4f}  {g:8.4f}  {dlt:10.4f}")
    prev = dlt

if ratios:
    print(f"\n  ratio medio Δ_n/Δ_{n+1} = {np.mean(ratios):.4f}  (target ψ={psi:.4f})")
    print(f"  scarto relativo |ratio/ψ - 1| = {abs(np.mean(ratios)/psi - 1):.4f}")
    if abs(np.mean(ratios)/psi - 1) < 0.15:
        print("  → ratio vicino a ψ: segnale compatibile con scaling superaureo")
    else:
        print("  → ratio NON vicino a ψ su questo modello (ipotesi non supportata qui)")

# ---------- 5. Predizione falsificabile ----------
print("\n[5] Predizione falsificabile")
print("""
  Se la RG del modello contiene ψ come fattore di blocco:
  1. Δ_lat(a/ψ) / Δ_lat(a) → 1/ψ
  2. Δ_phys = Δ_lat/a è stazionario lungo a_n = a0/ψ^n
  3. L'operatore di trasferimento di blocco ha raggio spettrale ψ

  Su KS vero (non solo plaquette):
  - costruire un blocco 3³→1 (o 3→1 in 1D ridotto)
  - stimare l'autovalore dominante del trasferimento
  - confrontare con ψ

  Se λ_max ≠ ψ e i ratio Δ non convergono a ψ,
  l'ipotesi “la rinormalizzazione contiene la sezione superaurea”
  è falsificata in quel setup.
""")

print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/psi_rg_results.npz",
         psi=psi, a_psi=a_psi, ratios=np.array(ratios) if ratios else np.array([]))
print("\nSalvato psi_rg_results.npz")
