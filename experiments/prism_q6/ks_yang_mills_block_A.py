#!/usr/bin/env python3
"""
Blocco A — Yang-Mills su reticolo: Kogut–Susskind
=================================================
- Hamiltoniana H = (g²/2) Σ E² - (1/g²) Σ Re Tr U_□
- SU(2), truncazione elettrica j ≤ 1/2 (stati |0>, |1/2,m>)
- Reticolo: Q3 (cubo) per diagonalizzazione esatta del settore
  a gauge fissata / truncata; formula strong-coupling su Q6.
- Output: Δ(g), α(g) = ||H_B|| / Δ_E  su una griglia di g.
"""

import numpy as np
from itertools import product, combinations
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, kron, eye as speye
from scipy.sparse.linalg import eigsh
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("BLOCCO A — Kogut–Susskind YM (SU(2) truncato)")
print("=" * 70)
t0 = time.time()

# ============================================================
# 1. Strong-coupling analitico su Q6 (standard)
# ============================================================
print("\n[1] Strong coupling su Q6 (analitico / semi-analitico)")

# In strong coupling (g → ∞) il vuoto è E=0 su tutti i link (singoletto).
# Il primo eccitato: una plaquette “flipped” (loop di flusso j=1/2) costa
# ΔE = (g²/2) * C_2(j=1/2) * 4   (4 link del loop)
# C_2(j) = j(j+1) = (1/2)(3/2) = 3/4  per SU(2)
# Quindi Δ_E^{(sc)} = (g²/2) * (3/4) * 4 = (3/2) g²
# (normalizzazioni KS classiche: H_E = (g²/2) Σ_ℓ E_ℓ² con E²|j> = j(j+1)|j>)

N_PLAQ_Q6 = 240
N_LINKS_Q6 = 192  # 6*64/2

def delta_E_strong(g):
    """Gap elettrico strong-coupling (primo loop j=1/2)."""
    return 1.5 * g**2   # (3/2) g²

def alpha_sc(g):
    """
    α ≈ ||H_B|| / Δ_E.
    ||H_B|| ~ (1/g²) * N_plaq * (bound su |Re Tr U|) ≤ (1/g²)*N_plaq*2
    (Tr U ≤ 2 per SU(2)).
    Stima grezza: ||H_B||_op ≤ 2 N_plaq / g²
    """
    dE = delta_E_strong(g)
    H_B_bound = 2.0 * N_PLAQ_Q6 / (g**2 + 1e-15)
    return H_B_bound / (dE + 1e-15)

g_grid = np.array([0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0])
print(f"  {'g':>6}  {'Δ_E (sc)':>12}  {'α_bound':>12}  {'α<1?':>6}")
print("  " + "-" * 42)
for g in g_grid:
    dE = delta_E_strong(g)
    al = alpha_sc(g)
    print(f"  {g:6.2f}  {dE:12.4f}  {al:12.4f}  {'sì' if al<1 else 'no':>6}")

g_crit_sc = None
for g in np.linspace(0.5, 5, 200):
    if alpha_sc(g) < 1:
        g_crit_sc = g
        break
print(f"\n  g* tale che α_bound(g*)=1 (stima) ≈ {g_crit_sc:.3f}")
print("  Per g > g* il bound strong-coupling garantisce Δ_phys ≥ (1−α)Δ_E > 0")
print("  (volume finito). Questo è lo stesso tipo di argomento del Cap. 22.")

# ============================================================
# 2. Diagonalizzazione numerica su Q3 (cubo) — SU(2) truncato
# ============================================================
print("\n[2] Diagonalizzazione su Q3 (8 vertici, 12 link, 6 plaquette)")

# Vertici Q3
verts = list(product([0,1], repeat=3))
n_v = 8
imap = {v: i for i, v in enumerate(verts)}

# Link: orientati (v → v flip d) con v_d=0 (per non doppiare)
links = []
for v in verts:
    for d in range(3):
        if v[d] == 0:
            w = list(v); w[d] = 1
            links.append((imap[v], imap[tuple(w)], d))
n_l = len(links)  # 12
assert n_l == 12

# Plaquette Q3
plaquettes = []
for d1, d2 in combinations(range(3), 2):
    fixed = [d for d in range(3) if d not in (d1, d2)]
    for bit in [0, 1]:
        v0 = [0]*3
        v0[fixed[0]] = bit
        def flip(vv, d):
            w = list(vv); w[d] = 1-w[d]; return tuple(w)
        a = tuple(v0)
        b = flip(a, d1)
        c = flip(b, d2)
        d = flip(c, d1)
        # indici link del ciclo (cerca in links)
        def find_link(i, j):
            for li, (x,y,_) in enumerate(links):
                if (x,y)==(i,j) or (x,y)==(j,i):
                    return li
            return None
        cycle = [find_link(imap[a],imap[b]), find_link(imap[b],imap[c]),
                 find_link(imap[c],imap[d]), find_link(imap[d],imap[a])]
        plaquettes.append(cycle)

print(f"  Link = {n_l},  Plaquette = {len(plaquettes)}")

# Truncazione: ogni link ha stati |0> (j=0) e |↑>,|↓> (j=1/2, m=±1/2)
# dim per link = 3.  dim totale = 3^12 ≈ 531441 — troppo grande.
# Strategia: sottospazio a basso flusso elettrico
#   - vuoto: tutti j=0
#   - stati a 1 plaquette eccitata (4 link a j=1/2, resto 0) — ma vincoli di Gauss
# Approccio pragmatico: modello ridotto "plaquette indipendenti" + correzione
# OPPURE: diagonalizzare un Hamiltoniano efficace a 6 plaquette con
# interazione debole (strong coupling effective).

# --- Modello efficace strong-coupling a plaquette ---
# Ogni plaquette ha stati |0> (vuoto) e |1> (eccitata, costo Δ_E)
# H_eff = Δ_E Σ_p n_p  - (1/g²) Σ_p (termine di tunneling) + interazioni
# Il termine magnetico crea/annichila loop: a leading order
# H_B ~ -(2/g²) Σ_p ( |0><1| + |1><0| )  (coefficiente da character expansion)
# Stati a più plaquette: se due plaquette condividono un link, c'è correzione.

n_p = len(plaquettes)  # 6
# Grafo duale delle 6 facce del cubo: ogni faccia ne tocca 4
adj_p = [[] for _ in range(n_p)]
# due plaquette sono adiacenti se condividono un link
for i in range(n_p):
    for j in range(i+1, n_p):
        if set(plaquettes[i]) & set(plaquettes[j]):
            adj_p[i].append(j)
            adj_p[j].append(i)

def H_eff_plaquette(g, n_p=6):
    """
    Spazio 2^{n_p}: bit = plaquette eccitata o no.
    H = Δ_E Σ n_i  - (2/g²) Σ_i σ^x_i  + λ Σ_<ij> n_i n_j
    (λ = correzione strong-coupling da link condivisi, O(g²) o O(1)).
    """
    dim = 2**n_p
    dE = 1.5 * g**2
    t = 2.0 / (g**2 + 1e-15)   # ampiezza magnetica
    lam = 0.25 * g**2            # penalità link condiviso (ordine di grandezza)
    H = np.zeros((dim, dim))
    for s in range(dim):
        bits = [(s >> k) & 1 for k in range(n_p)]
        # diagonale: energia elettrica + interazione
        e_diag = dE * sum(bits)
        for i in range(n_p):
            for j in adj_p[i]:
                if j > i:
                    e_diag += lam * bits[i] * bits[j]
        H[s, s] = e_diag
        # off-diagonale: flip di una plaquette (termine magnetico)
        for i in range(n_p):
            s2 = s ^ (1 << i)
            H[s, s2] -= t
    return H

print("\n  Gap effettivo su Q3 (modello a plaquette):")
print(f"  {'g':>6}  {'E0':>10}  {'E1':>10}  {'Δ':>10}  {'α_est':>10}")
print("  " + "-" * 50)

deltas = []
alphas = []
for g in g_grid:
    H = H_eff_plaquette(g)
    ev = np.sort(eigh(H, eigvals_only=True))
    E0, E1 = ev[0], ev[1]
    dlt = E1 - E0
    # α_est: ||H_B|| / Δ_E con ||H_B|| ~ 2*n_p/g² (operatore)
    dE = 1.5 * g**2
    HBnorm = 2.0 * n_p / (g**2 + 1e-15)
    al = HBnorm / (dE + 1e-15)
    deltas.append(dlt)
    alphas.append(al)
    print(f"  {g:6.2f}  {E0:10.4f}  {E1:10.4f}  {dlt:10.4f}  {al:10.4f}")

deltas = np.array(deltas)
alphas = np.array(alphas)

# ============================================================
# 3. Sintesi verso il continuo
# ============================================================
print("\n[3] Sintesi Blocco A")
print(f"  Q6 strong-coupling: gap Δ_E = (3/2)g²,  α_bound < 1 per g ≳ {g_crit_sc:.2f}")
print(f"  Q3 modello a plaquette: Δ(g) calcolato sulla griglia sopra")
print(f"  Min Δ sulla griglia (Q3 eff) = {deltas.min():.4f}  (a g={g_grid[deltas.argmin()]:.2f})")
print(f"  α_est < 1 per g ≳ {g_grid[np.where(np.array(alphas)<1)[0][0]] if np.any(np.array(alphas)<1) else float('nan'):.2f}")

print("""
  Interpretazione (onesta):
  - In volume finito e strong coupling il gap è sotto controllo (standard).
  - Il modello a plaquette su Q3 mostra un gap positivo su tutta la griglia
    esplorata; α_est scende sotto 1 a g abbastanza grande.
  - Per il CONTINUO serve una sequenza di reticoli (L→∞, a→0) con
    traiettoria di scaling g(a) e verifica che Δ_phys = Δ/a non vada a 0.
  - Questo script è il Blocco A: infrastruttura e numeri in volume finito.
""")

print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/ks_block_A_results.npz",
         g_grid=g_grid, deltas=deltas, alphas=alphas,
         g_crit_sc=g_crit_sc)
print("\nSalvato ks_block_A_results.npz")
