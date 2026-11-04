#!/usr/bin/env python3
"""
Cerca definizioni di H_E, H_B e blocco con s_E ≈ s_B
====================================================
1. H_E = laplaciano di link (1-forme) / laplaciano di grafo
2. H_B = laplaciano di plaquette (2-forme) = d1 d1^T  (coboundary)
   → entrambi dalla stessa successione di cobordo del complesso
3. Blocco: proiezione di Hodge / restrizione coerente su cocatene
4. Misura s_E, s_B, α su una torre di coarsening uniforme
"""

import numpy as np
from itertools import product, combinations
from scipy.sparse import lil_matrix, csr_matrix, diags
from scipy.sparse.linalg import norm as snorm, eigsh
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("RICERCA DI SCALING OMOGENEO s_E = s_B")
print("=" * 70)
t0 = time.time()

psi = float((1 + ((29+3*np.sqrt(93))/2)**(1/3) + ((29-3*np.sqrt(93))/2)**(1/3)) / 3)

# ============================================================
# Complesso cellulare di un prodotto di cicli
# ============================================================
def build_complex(cycles):
    """Vertici, link orientati, plaquette orientate."""
    verts = list(product(*[range(c) for c in cycles]))
    vmap = {v: i for i, v in enumerate(verts)}
    n_v = len(verts)
    D = len(cycles)

    # Link orientati: (v, d) con passo +1 lungo d
    links = []  # (id_link, v_from, v_to, direction)
    for v in verts:
        for d in range(D):
            w = list(v)
            w[d] = (w[d] + 1) % cycles[d]
            links.append((vmap[v], vmap[tuple(w)], d))
    n_l = len(links)

    # Plaquette orientate: (v, d1, d2) d1<d2
    plaq = []
    for v in verts:
        for d1, d2 in combinations(range(D), 2):
            def stp(vv, d, s=+1):
                ww = list(vv)
                ww[d] = (ww[d] + s) % cycles[d]
                return tuple(ww)
            # 4 edge del bordo orientato
            # +e(v,d1), +e(v+d1,d2), -e(v+d2,d1), -e(v,d2)
            v1 = stp(v, d1)
            v2 = stp(v, d2)
            v12 = stp(v1, d2)
            edges = [
                (vmap[v], vmap[v1], d1, +1),
                (vmap[v1], vmap[v12], d2, +1),
                (vmap[v2], vmap[v12], d1, -1),
                (vmap[v], vmap[v2], d2, -1),
            ]
            plaq.append(edges)
    n_p = len(plaq)
    return verts, vmap, links, plaq, n_v, n_l, n_p, cycles

def coboundary_matrices(links, plaq, n_v, n_l, n_p):
    """
    d0: 0-cochains → 1-cochains  (n_l × n_v)
    d1: 1-cochains → 2-cochains  (n_p × n_l)
    H_E ~ d0 d0^T  (laplaciano su 1-forme) oppure d0^T d0 (su 0-forme)
    H_B ~ d1 d1^T  (laplaciano su 2-forme)
    Usiamo:
      H_E = d0^T d0   (laplaciano di grafo, 0-forme)  [elettrico / hop]
      H_B = d1 d1^T   (laplaciano di plaquette, 2-forme) [magnetico]
    Entrambi dalla stessa successione di cobordo → stessa geometria.
    """
    # mappa (i,j,d) → id link: ogni link orientato è univoco (v,d)
    # links[k] = (i, j, d) con j = i + e_d
    d0 = lil_matrix((n_l, n_v), dtype=float)
    for k, (i, j, d) in enumerate(links):
        d0[k, i] = -1.0
        d0[k, j] = +1.0
    d0 = d0.tocsr()

    # per d1 serve id del link (v,d)
    link_id = {}
    for k, (i, j, d) in enumerate(links):
        link_id[(i, d)] = k  # partenza i, direzione d

    d1 = lil_matrix((n_p, n_l), dtype=float)
    for p, edges in enumerate(plaq):
        for (i, j, d, sgn) in edges:
            # link id: se sgn=+1, link parte da i dir d; se -1, link parte da i dir d
            # edges già hanno (from, to, d, sgn) con from --d--> to quando sgn=+1
            if sgn == +1:
                key = (i, d)
            else:
                # -e: bordo opposto, stesso link geometrico
                key = (i, d)
            if key not in link_id:
                # prova from = j, direzione d (orientazione inversa)
                key = (j, d)
            if key in link_id:
                d1[p, link_id[key]] += sgn
    d1 = d1.tocsr()

    HE = (d0.T @ d0).tocsr()          # n_v × n_v
    HB = (d1 @ d1.T).tocsr()          # n_p × n_p
    # Per confrontare norme su spazi diversi, usiamo anche
    # H_B su 1-forme: d1^T d1  (n_l × n_l)  — stesso grado "link"
    HB_link = (d1.T @ d1).tocsr()     # n_l × n_l
    HE_link = (d0 @ d0.T).tocsr()     # n_l × n_l
    return HE, HB, HE_link, HB_link, d0, d1

def fro(M):
    return snorm(M)

# ============================================================
# Torre di coarsening UNIFORME (tutti i cicli / factor)
# ============================================================
print("\n[1] Torre di complessi con coarsening uniforme")

# Partenza: C6^5  (piccolo ma uniforme)
base_cycles = [6, 6, 6, 6, 6]
factors = [1, 2, 3]  # cycles → cycles/factor

results = []
print(f"  {'factor':>6}  {'cycles':>20}  {'n_v':>6}  {'n_l':>6}  {'n_p':>6}  "
      f"{'||HE||':>10}  {'||HB_l||':>10}  {'||HE_l||':>10}  {'||HB||':>10}")
print("  " + "-" * 100)

for f in factors:
    cyc = [max(2, c // f) for c in base_cycles]  # min 2 per restare cicli
    verts, vmap, links, plaq, n_v, n_l, n_p, cyc = build_complex(cyc)
    HE, HB, HE_l, HB_l, d0, d1 = coboundary_matrices(links, plaq, n_v, n_l, n_p)
    nHE, nHB = fro(HE), fro(HB)
    nHEl, nHBl = fro(HE_l), fro(HB_l)
    results.append(dict(f=f, cyc=cyc, n_v=n_v, n_l=n_l, n_p=n_p,
                        nHE=nHE, nHB=nHB, nHEl=nHEl, nHBl=nHBl))
    print(f"  {f:6d}  {str(cyc):>20}  {n_v:6d}  {n_l:6d}  {n_p:6d}  "
          f"{nHE:10.4f}  {nHBl:10.4f}  {nHEl:10.4f}  {nHB:10.4f}")

# ============================================================
# Stima s su coppie di livelli
# ============================================================
print("\n[2] Stima s_E, s_B (norme sullo stesso grado: 1-forme)")

def s_from_ratio(r, base=psi):
    return np.log(r) / np.log(base)

print(f"  {'step':>10}  {'r_E(link)':>12}  {'r_B(link)':>12}  "
      f"{'s_E(ψ)':>10}  {'s_B(ψ)':>10}  {'Δs':>10}  {'r_E/r_B':>10}")
print("  " + "-" * 80)

ds_list = []
alpha_link = []
for i in range(len(results)-1):
    a, b = results[i], results[i+1]
    rE = a['nHEl'] / (b['nHEl'] + 1e-30)
    rB = a['nHBl'] / (b['nHBl'] + 1e-30)
    sE = s_from_ratio(rE)
    sB = s_from_ratio(rB)
    ds = sE - sB
    ds_list.append(abs(ds))
    print(f"  f{a['f']}→f{b['f']:}  {rE:12.6f}  {rB:12.6f}  "
          f"{sE:10.4f}  {sB:10.4f}  {ds:10.4f}  {rE/(rB+1e-30):10.6f}")

# α su 1-forme
print("\n[3] α su 1-forme (stesso spazio)")
print(f"  {'factor':>6}  {'α = ||HB_l||/||HE_l||':>22}")
for r in results:
    al = r['nHBl'] / (r['nHEl'] + 1e-30)
    alpha_link.append(al)
    print(f"  {r['f']:6d}  {al:22.8f}")

var_alpha = (max(alpha_link)-min(alpha_link))/(np.mean(alpha_link)+1e-15)
print(f"  variazione relativa α: {var_alpha:.6f}")

# ============================================================
# Normalizzazione per volume / densità
# ============================================================
print("\n[4] Norme densità ( / sqrt(dim) ) — scaling intensivo")
print(f"  {'factor':>6}  {'||HE_l||/√n_l':>14}  {'||HB_l||/√n_l':>14}  {'α_dens':>10}")
dens_E, dens_B = [], []
for r in results:
    dE = r['nHEl'] / np.sqrt(r['n_l'])
    dB = r['nHBl'] / np.sqrt(r['n_l'])
    dens_E.append(dE); dens_B.append(dB)
    print(f"  {r['f']:6d}  {dE:14.6f}  {dB:14.6f}  {dB/(dE+1e-30):10.6f}")

print("\n  s da densità:")
for i in range(len(results)-1):
    rE = dens_E[i]/(dens_E[i+1]+1e-30)
    rB = dens_B[i]/(dens_B[i+1]+1e-30)
    print(f"  f{results[i]['f']}→f{results[i+1]['f']}:  "
          f"s_E={s_from_ratio(rE):.4f}  s_B={s_from_ratio(rB):.4f}  "
          f"Δs={s_from_ratio(rE)-s_from_ratio(rB):.4f}  rE/rB={rE/(rB+1e-30):.6f}")

# ============================================================
# Verdetto
# ============================================================
print("\n" + "=" * 70)
print("VERDETTO")
print("=" * 70)
mean_ds = np.mean(ds_list) if ds_list else float('nan')
print(f"""
Definizione trovata (successione di cobordo del prodotto di cicli):
  H_E = d0 d0^T   (laplaciano su 1-forme / link)
  H_B = d1^T d1   (laplaciano su 1-forme da bordo plaquette)

Stesso spazio (link), stessa torre di coarsening uniforme.

  |s_E - s_B| medio (log_ψ) ≈ {mean_ds:.4f}
  variazione relativa di α   ≈ {var_alpha:.4f}
""")

if mean_ds < 0.15 and var_alpha < 0.15:
    print("  → s_E ≈ s_B e α quasi invariante: (A2) SUPPORTATA su questo complesso.")
elif mean_ds < 0.4:
    print("  → s_E ≈ s_B in modo approssimato: (A2) parzialmente supportata.")
else:
    print("  → s_E e s_B ancora separati: serve altro blocco o normalizzazione.")

print(f"\n  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/homogeneous_scaling_results.npz",
         alpha_link=np.array(alpha_link), dens_E=np.array(dens_E), dens_B=np.array(dens_B),
         psi=psi)
print("\nSalvato homogeneous_scaling_results.npz")
