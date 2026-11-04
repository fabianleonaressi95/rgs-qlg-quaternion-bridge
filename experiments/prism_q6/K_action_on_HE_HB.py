#!/usr/bin/env python3
"""
Azione di K su operatori di tipo elettrico e magnetico
======================================================
Obiettivo: stimare gli esponenti di scaling s_E, s_B sotto lo step
a → a/ψ  (equivalente a coniugazione / blocco tramite K).

Costruzione:
1. K = companion di x³−x²−1 sul fibrato di scale (3 livelli)
2. H_E, H_B come operatori su (scala ⊗ spazio plaquette/link)
3. Step di blocco: proiezione scala n → scala n+1 con fattore ψ
4. Stima s da  ||H^{(n+1)}|| / ||H^{(n)}|| ≈ ψ^{-s}
"""

import numpy as np
from scipy.linalg import eigh, expm
import warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("AZIONE DI K SU H_E E H_B — stima di s_E, s_B")
print("=" * 70)

psi = float((1 + ((29+3*np.sqrt(93))/2)**(1/3) + ((29-3*np.sqrt(93))/2)**(1/3)) / 3)
print(f"ψ = {psi:.15f}")

# ============================================================
# 1. Operatore K (companion) e potenze
# ============================================================
K = np.array([
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 1.0]
], dtype=float)
assert np.linalg.norm(K@K@K - K@K - np.eye(3)) < 1e-12

# Autodecomposizione: K = V diag(ψ, ω, ω̄) V^{-1}
evals, V = np.linalg.eig(K)
# ordina: ψ reale prima
idx = np.argsort(-np.real(evals))
evals, V = evals[idx], V[:, idx]
print(f"autovalori K: {evals}")

# ============================================================
# 2. Spazio: 3 scale ⊗ N modi (link / plaquette)
# ============================================================
N_MODES = 12  # proxy: 12 link di Q3 o 6 plaquette×2
N_SCALE = 3

# H_E diagonale in base di modi: energie elettriche
# Al livello di scala 0: E² ~ 0.75 su modi "eccitati", 0 sul vuoto
# Distribuiamo autovalori elettrici tipici
e_electric = np.array([0.0] + [0.75]*4 + [1.5]*4 + [2.25]*3)  # 12 valori
e_electric = e_electric[:N_MODES]

# H_B: matrice densa che accoppia i modi (magnetico / plaquette)
# Costruita come laplaciano di interazione + termine off-diagonale
rng = np.random.RandomState(7)
A = rng.randn(N_MODES, N_MODES)
A = 0.5*(A+A.T)
# spettro controllato
evA, UA = eigh(A)
A = UA @ np.diag(evA - evA.min() + 0.1) @ UA.T  # positivo

def build_H_on_scales(g, a_scales):
    """
    H_E e H_B come operatori su (scale ⊗ modi).
    Scaling naive di prova:
      H_E ~ (g²/2) / a   ×  diag(e_electric)   su ogni scala
      H_B ~ (1/g²) / a   ×  A
    (poi confronteremo con lo scaling indotto da K)
    """
    dim = N_SCALE * N_MODES
    HE = np.zeros((dim, dim))
    HB = np.zeros((dim, dim))
    for s, a in enumerate(a_scales):
        sl = slice(s*N_MODES, (s+1)*N_MODES)
        HE[sl, sl] = (g**2 / 2.0) / (a + 1e-15) * np.diag(e_electric)
        HB[sl, sl] = (1.0 / (g**2 + 1e-15)) / (a + 1e-15) * A
    return HE, HB

# Scale gerarchiche legate a ψ
a_scales = np.array([1.0, 1/psi, 1/psi**2])
g = 2.0
HE, HB = build_H_on_scales(g, a_scales)

print(f"\n[1] Operatori su {N_SCALE} scale × {N_MODES} modi")
print(f"  ||H_E||_F = {np.linalg.norm(HE):.6f}")
print(f"  ||H_B||_F = {np.linalg.norm(HB):.6f}")
print(f"  α_naive   = ||H_B||/||H_E|| = {np.linalg.norm(HB)/(np.linalg.norm(HE)+1e-15):.6f}")

# ============================================================
# 3. Azione di K: coniugazione sul settore scale
# ============================================================
print("\n[2] Azione di K (blocco sulle scale)")

# K agisce sulle scale: O ↦ (K ⊗ I) O (K ⊗ I)^{-1}
# oppure, per operatori diagonali a blocchi, rimescola i livelli
K_big = np.kron(K, np.eye(N_MODES))
K_big_inv = np.linalg.inv(K_big)

def conjugate(O):
    return K_big @ O @ K_big_inv

HE1 = conjugate(HE)
HB1 = conjugate(HB)

def block_norms(O):
    """Norme di Frobenius su ciascun blocco di scala."""
    norms = []
    for s in range(N_SCALE):
        sl = slice(s*N_MODES, (s+1)*N_MODES)
        norms.append(np.linalg.norm(O[sl, sl]))
    return np.array(norms)

nE0 = block_norms(HE)
nE1 = block_norms(HE1)
nB0 = block_norms(HB)
nB1 = block_norms(HB1)

print(f"  ||H_E|| per scala PRIMA di K: {nE0}")
print(f"  ||H_E|| per scala DOPO  di K: {nE1}")
print(f"  ||H_B|| per scala PRIMA di K: {nB0}")
print(f"  ||H_B|| per scala DOPO  di K: {nB1}")

# ============================================================
# 4. Stima di s_E, s_B dal rapporto di norme globali
# ============================================================
print("\n[3] Stima s da rapporti di norme")

# Sotto a → a/ψ ci aspettiamo ||H|| → ψ^{-s} ||H||
# Confrontiamo norme totali prima/dopo un step di gerarchia
# Metodo A: rapporto tra scala 0 e scala 1 NELLO STESSO operatore
# (dove a1/a0 = 1/ψ per costruzione)
sE_from_hierarchy = np.log(nE0[0]/(nE0[1]+1e-30)) / np.log(psi)
sB_from_hierarchy = np.log(nB0[0]/(nB0[1]+1e-30)) / np.log(psi)
print(f"  Dalla gerarchia a_n (naive 1/a):")
print(f"    s_E ≈ {sE_from_hierarchy:.6f}")
print(f"    s_B ≈ {sB_from_hierarchy:.6f}")
print(f"    s_E - s_B ≈ {sE_from_hierarchy - sB_from_hierarchy:.6e}")

# Metodo B: dopo coniugazione con K
# Se K mescola le scale "alla Krist-All", le norme possono redistribuirsi
rE = np.linalg.norm(HE1) / (np.linalg.norm(HE)+1e-30)
rB = np.linalg.norm(HB1) / (np.linalg.norm(HB)+1e-30)
print(f"\n  Dopo coniugazione K O K^{{-1}}:")
print(f"    ||H_E'||/||H_E|| = {rE:.6f}")
print(f"    ||H_B'||/||H_B|| = {rB:.6f}")

# ============================================================
# 5. Versione Krist-All: scaling OMOGENEO forzato da K
# ============================================================
print("\n[4] Versione Krist-All: entrambi gli operatori trasformano come K")

# Idea: H_E e H_B non sono indipendenti; sono sezioni dello stesso
# fibrato di operatori su PRISM, e lo step di blocco è
#   H ↦ ψ^{-1} Ad_K(H)   (stesso fattore per tutti)
# Allora s_E = s_B = 1 per definizione, e α è invariante.

def krist_step(O, s=1.0):
    """Step di blocco Krist-All: coniugazione + riscala ψ^{-s}."""
    return (psi**(-s)) * conjugate(O)

HE_k = krist_step(HE, s=1.0)
HB_k = krist_step(HB, s=1.0)
alpha0 = np.linalg.norm(HB) / (np.linalg.norm(HE)+1e-15)
alpha1 = np.linalg.norm(HB_k) / (np.linalg.norm(HE_k)+1e-15)
print(f"  α prima dello step K: {alpha0:.8f}")
print(f"  α dopo  dello step K: {alpha1:.8f}")
print(f"  |α' - α| = {abs(alpha1-alpha0):.3e}")
if abs(alpha1-alpha0) < 1e-10:
    print("  → α ESATTAMENTE invariante sotto lo step Krist-All (s_E=s_B=1)")
else:
    print("  → α non invariante (inatteso se s uguali)")

# Iterazione di più step
alpha_hist = [alpha0]
O_E, O_B = HE.copy(), HB.copy()
for step in range(6):
    O_E = krist_step(O_E, s=1.0)
    O_B = krist_step(O_B, s=1.0)
    alpha_hist.append(np.linalg.norm(O_B)/(np.linalg.norm(O_E)+1e-15))
print(f"  α lungo 6 step Krist-All: {[f'{a:.6f}' for a in alpha_hist]}")

# ============================================================
# 6. Confronto con scaling naive YM (s_E e s_B diversi via g(a))
# ============================================================
print("\n[5] Confronto: scaling YM naive (g dipende da a)")

def alpha_naive_YM(a, g0=2.5, b0=22/3):
    inv = 1/g0**2 + (b0/(8*np.pi**2))*np.log(1/(a+1e-30))
    g2 = 1/max(inv, 1e-8)
    # H_E ~ g²/a, H_B ~ 1/(g² a)  → α ~ 1/g^4
    return 1.0 / (g2**2 + 1e-15)

print(f"  {'a':>10}  {'α_YM_naive':>12}  {'α_Krist':>12}")
print("  " + "-" * 38)
for n in range(6):
    a = 1.0/psi**n
    print(f"  {a:10.5f}  {alpha_naive_YM(a):12.6f}  {alpha0:12.6f}")

# ============================================================
# 7. Verdetto
# ============================================================
print("\n" + "=" * 70)
print("VERDETTO")
print("=" * 70)
print("""
1. Con scaling naive 1/a uguale per H_E e H_B:
   s_E = s_B = 1  e α è invariante per costruzione.

2. Sotto lo step Krist-All H ↦ ψ^{-1} Ad_K(H) (stesso s):
   α resta ESATTAMENTE costante lungo la torre a_n = a0/ψ^n.

3. In YM continuum standard α_YM ~ 1/g(a)^4 VARIA lungo il flusso.
   La differenza è tutta in (A2): s_E = s_B.

4. Quindi:
   - SE si assume (o si dimostra) che K agisce allo stesso modo
     su H_E e H_B, ALLORA α è invariante e m_phys>0 segue
     dalla serie convergente + α<1.
   - La dimostrazione di (A2) deve venire dalla definizione di
     H_E, H_B come operatori sul complesso PRISM e dall'azione
     di K su quel complesso — non dalla companion matrix da sola.

5. Prossimo passo matematico:
   definire H_E, H_B come elementi di un'algebra di operatori
   su PRISM e mostrare che Ad_K li scala con lo stesso peso.
""")
print("=" * 70)

np.savez("/home/workdir/artifacts/K_action_sE_sB_results.npz",
         psi=psi, alpha0=alpha0, alpha_hist=np.array(alpha_hist),
         sE=sE_from_hierarchy, sB=sB_from_hierarchy)
print("\nSalvato K_action_sE_sB_results.npz")
