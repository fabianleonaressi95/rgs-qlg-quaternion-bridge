#!/usr/bin/env python3
"""
Test genuino Kogut-Susskind SU(2) sul cubo Q3 (8 vertici, 12 link, 6 plaquette)
================================================================================
Hamiltoniana:  H = (g^2/2) sum_l E_l^2  -  (1/g^2) sum_p Tr U_p      (+ costante)
Troncazione:   j <= 1/2 su ogni link.
Spazio fisico: stati gauge-invarianti ESATTI.  Il cubo ha vertici di grado 3, quindi
               con j in {0,1/2} ogni vertice ha 0 o 2 linee occupate: gli stati sono
               i sottografi pari (spazio dei cicli, dim 2^5 = 32), |S> = prod_cicli Tr U_c,
               ortonormali (integrale di Haar).
Elementi di matrice di Tr U_p: integrali di Haar ESATTI (contrazione di tensori con
               int U_ab U_cd = eps_ac eps_bd / 2, int U_ab conj(U_cd) = delta_ac delta_bd / 2),
               VALIDATI con un campionamento Monte Carlo indipendente su SU(2)^12.

Nessuna etichetta predefinita: tutto e' calcolato.
"""
import numpy as np
from itertools import product, combinations
from scipy.linalg import eigh
import time, sys

rng = np.random.default_rng(1)

# ------------------------------------------------------------------ geometria
verts = list(product([0, 1], repeat=3))
vid = {v: i for i, v in enumerate(verts)}
links = []                                   # (da, a, direzione)
for v in verts:
    for d in range(3):
        if v[d] == 0:
            w = list(v); w[d] = 1
            links.append((vid[v], vid[tuple(w)], d))
NL = len(links)                              # 12
lid = {}
for k, (a, b, d) in enumerate(links):
    lid[(a, b)] = (k, +1)                    # percorso lungo l'orientazione
    lid[(b, a)] = (k, -1)                    # percorso contro l'orientazione

def path_cycle(vs):
    return [lid[(vs[i], vs[(i + 1) % len(vs)])] for i in range(len(vs))]

plaqs = []
for d1, d2 in combinations(range(3), 2):
    d3 = 3 - d1 - d2
    for b in (0, 1):
        v0 = [0, 0, 0]; v0[d3] = b
        v1 = v0.copy(); v1[d1] = 1
        v12 = v1.copy(); v12[d2] = 1
        v2 = v0.copy(); v2[d2] = 1
        plaqs.append(path_cycle([vid[tuple(x)] for x in (v0, v1, v12, v2)]))
pmask = [sum(1 << k for k, _ in cy) for cy in plaqs]
assert len(plaqs) == 6

def is_even(mask):
    deg = [0] * 8
    for k, (a, b, d) in enumerate(links):
        if mask >> k & 1:
            deg[a] += 1; deg[b] += 1
    return all(x % 2 == 0 for x in deg)

states = [m for m in range(1 << NL) if is_even(m)]
sidx = {m: i for i, m in enumerate(states)}
NS = len(states)
assert NS == 32, NS

def cycles_of(mask):
    inc = {v: [] for v in range(8)}
    for k, (a, b, d) in enumerate(links):
        if mask >> k & 1:
            inc[a].append(k); inc[b].append(k)
    used, out = set(), []
    for k0 in range(NL):
        if (mask >> k0 & 1) and k0 not in used:
            a, b, _ = links[k0]
            vs, cur, e = [a], a, k0
            while True:
                used.add(e)
                ea, eb, _ = links[e]
                nxt = eb if cur == ea else ea
                if nxt == vs[0]:
                    break
                vs.append(nxt); cur = nxt
                e = [kk for kk in inc[cur] if kk != e][0]
            out.append(path_cycle(vs))
    return out

# --------------------------------------------------- integrale di Haar esatto
eps = np.array([[0., 1.], [-1., 0.]])
dl = np.eye(2)
P_same = 0.5 * np.einsum('ac,bd->abcd', eps, eps)   # (U,U) o (U^dag,U^dag)
P_opp = 0.5 * np.einsum('ad,bc->abcd', dl, dl)      # (U,U^dag)

def haar(loops):
    slots, labels, c = {}, [], 0
    for m, cy in enumerate(loops):
        labels.append(list(range(c, c + len(cy)))); c += len(cy)
        for t, (k, s) in enumerate(cy):
            slots.setdefault(k, []).append((m, t, s))
    ops = []
    for k, sl in slots.items():
        if len(sl) % 2:
            return 0.0                               # simmetria del centro Z2
        assert len(sl) == 2
        (m1, t1, s1), (m2, t2, s2) = sl
        L1, L2 = labels[m1], labels[m2]
        r1, c1 = L1[(t1 - 1) % len(L1)], L1[t1]
        r2, c2 = L2[(t2 - 1) % len(L2)], L2[t2]
        ops += [P_same if s1 == s2 else P_opp, [r1, c1, r2, c2]]
    return float(np.einsum(*ops, [], optimize='greedy'))

t0 = time.time()
Bp = np.zeros((6, NS, NS))                           # Tr U_p proiettato su j<=1/2
for p in range(6):
    for j, S in enumerate(states):
        S2 = S ^ pmask[p]
        i = sidx[S2]
        Bp[p, i, j] = haar(cycles_of(S2) + [plaqs[p]] + cycles_of(S))
B = Bp.sum(0)
size = np.array([bin(S).count('1') for S in states], float)

print("=" * 72)
print("KS SU(2) SUL CUBO Q3  -  j<=1/2, spazio gauge-invariante esatto")
print("=" * 72)
print(f"stati gauge-invarianti: {NS}   (attesi 2^5=32)   [{time.time()-t0:.1f}s]")
print(f"B simmetrica: {np.allclose(B, B.T)}   <p|TrU_p|0> = {Bp[0, sidx[pmask[0]], sidx[0]]:.6f} (atteso 1)")

# ---------------------------------------- validazione Monte Carlo indipendente
def rand_su2(n):
    q = rng.normal(size=(n, 4)); q /= np.linalg.norm(q, axis=1, keepdims=True)
    a, b, c, d = q.T
    U = np.zeros((n, 2, 2), complex)
    U[:, 0, 0] = a + 1j * b; U[:, 0, 1] = c + 1j * d
    U[:, 1, 0] = -c + 1j * d; U[:, 1, 1] = a - 1j * b
    return U

NMC = 200000
Ul = [rand_su2(NMC) for _ in range(NL)]
def tr_cycle(cy):
    M = None
    for k, s in cy:
        W = Ul[k] if s > 0 else np.conj(np.transpose(Ul[k], (0, 2, 1)))
        M = W if M is None else M @ W
    return np.real(np.trace(M, axis1=1, axis2=2))
fS = np.array([np.prod([tr_cycle(c) for c in cycles_of(S)], axis=0) if S else np.ones(NMC)
               for S in states])
trp = np.array([tr_cycle(c) for c in plaqs])
Bmc = np.zeros((NS, NS)); SE = np.zeros((NS, NS)); Bfull_norm2 = np.zeros(NS)
for j in range(NS):
    gS = trp.sum(0) * fS[j]
    for i in range(NS):
        x = fS[i] * gS
        Bmc[i, j] = x.mean(); SE[i, j] = x.std() / np.sqrt(NMC)
    Bfull_norm2[j] = np.mean(gS ** 2)
err = np.max(np.abs(Bmc - B))
z = np.abs(Bmc - B) / np.maximum(SE, 1e-12)
print(f"validazione MC (N={NMC}):  max|B_exact - B_MC| = {err:.4f};  max z-score = {z.max():.2f}"
      f"  (frazione z<3: {np.mean(z < 3):.4f}, elementi nulli esatti riprodotti: "
      f"{np.mean(np.abs(Bmc[B == 0]) < 4 * SE[B == 0].clip(1e-12)):.4f})")
leak = 1.0 - (B ** 2).sum(0) / Bfull_norm2           # peso fuori dallo spazio j<=1/2
print(f"peso di (sum_p TrU_p)|S> FUORI dalla troncazione j<=1/2: media {leak.mean():.3f}, max {leak.max():.3f}")

# --------------------------------------------------------- modello 'giocattolo'
adj = [[j for j in range(6) if j != i and (pmask[i] & pmask[j])] for i in range(6)]
def toy_gap(g):
    dim = 64; H = np.zeros((dim, dim))
    for s in range(dim):
        bits = [(s >> k) & 1 for k in range(6)]
        e = 1.5 * g**2 * sum(bits)
        for i in range(6):
            for j in adj[i]:
                if j > i: e += 0.25 * g**2 * bits[i] * bits[j]
        H[s, s] = e
        for i in range(6): H[s, s ^ (1 << i)] -= 2.0 / g**2
    ev = np.linalg.eigvalsh(H)
    return ev[1] - ev[0]

# ------------------------------------------------------------- spettro vs g
normB = np.max(np.abs(np.linalg.eigvalsh(B)))
print(f"\n||B||_op (g-indipendente) = {normB:.6f}")
dE_unit = 0.5 * 0.75 * 4                              # Delta_E = 1.5 g^2 (loop di 4 link)
print(f"\n{'g':>5} {'E0':>10} {'Delta_real':>11} {'Delta_toy':>10} {'Delta_E':>9} {'||H_B||':>9}"
      f" {'alpha':>9} {'Weyl':>9} {'Weyl ok':>8}")
rows = []
for g in [0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
    H = 0.5 * g**2 * 0.75 * np.diag(size) - B / g**2
    ev = np.linalg.eigvalsh(H)
    D = ev[1] - ev[0]
    dE = dE_unit * g**2
    hb = normB / g**2
    alpha = hb / dE
    weyl = dE - 2 * hb
    rows.append((g, alpha))
    print(f"{g:5.2f} {ev[0]:10.4f} {D:11.4f} {toy_gap(g):10.4f} {dE:9.3f} {hb:9.4f}"
          f" {alpha:9.4f} {weyl:9.3f} {'si' if D >= weyl - 1e-9 else 'NO':>8}")

gs = np.array([r[0] for r in rows]); al = np.array([r[1] for r in rows])
slope = np.polyfit(np.log(gs), np.log(al), 1)[0]
g_w = (normB / 0.75) ** 0.25
print(f"\npendenza log(alpha) vs log(g) = {slope:.3f}   (alpha = ||B||/(1.5 g^4): atteso -4)")
print(f"bound di Weyl  Delta >= Delta_E - 2||H_B||  positivo per g > {g_w:.3f}")
print(f"(il vecchio bound con (1-alpha) e coefficiente 2N_p/g^2 dava g > 4.23)")
print("""
LETTURA
 - H_E = g^2 * (matrice indipendente da g),  H_B = (1/g^2) * B  =>  alpha(g) = c/g^4 ESATTAMENTE.
   Lungo una traiettoria di continuo g(a) -> 0, alpha -> infinito: alpha NON e' invariante
   di scala e la condizione alpha<1 vale solo nel regime di accoppiamento forte.
 - La troncazione j<=1/2 perde peso O(1) per g piccolo (vedi 'peso fuori troncazione'):
   per g <~ 1.5 i valori di Delta_real NON sono affidabili.
 - Un volume (1 cubo) e una taglia: nessuna informazione sul limite a -> 0.
""")
