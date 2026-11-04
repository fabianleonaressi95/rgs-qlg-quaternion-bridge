#!/usr/bin/env python3
"""
Quattro passi falsificabili
===========================
1. Immersione Plücker che deve FALLIRE su casuali e (eventualmente) passare su PRISM
2. QFIM: derivate lungo i cicli, fit a convergenza, confronto omogeneo, split strutturale
3. Torsione tipo Weitzenböck e splitting non riducibile a [L,J]=0
4. Adiacenza reale plaquette Q6 + azione di gauge + stima σ da area law
"""

import numpy as np
from itertools import combinations, product
from scipy.linalg import eigh
from scipy.optimize import minimize
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("PASSI FALSIFICABILI 1–2–3–4")
print("=" * 70)
t0 = time.time()

# ============================================================
# 1. IMMERSIONE PLÜCKER NON BANALE (deve fallire sui casuali)
# ============================================================
print("\n[1] Immersione Plücker strutturale (non curva dei momenti)")

N, N_EFF, K = 8192, 6, 2

def prism_coords(n):
    idx = np.random.choice(69120, n, replace=False)
    e = idx % 5; tmp = idx // 5
    c = tmp % 6; tmp //= 6
    b = tmp % 6; tmp //= 6
    a = tmp % 6; i64 = tmp // 6
    return np.stack([i64, a, b, c, e], 1).astype(float)

def immersion_structured(coords):
    """
    Immersione che usa relazioni discrete tra i cicli
    (differenze, prodotti mod, parity) — non monotona per costruzione.
    Se è PRISM-specifica, soft-pos su casuali deve essere ~0.5.
    """
    x = coords  # (n, 5)
    # 6 colonne costruite da relazioni non monotone
    cols = np.zeros((coords.shape[0], N_EFF))
    cols[:, 0] = np.sin(2*np.pi*x[:, 0]/64) * np.cos(2*np.pi*x[:, 1]/6)
    cols[:, 1] = np.cos(2*np.pi*x[:, 0]/64) * np.sin(2*np.pi*x[:, 2]/6)
    cols[:, 2] = np.sin(2*np.pi*(x[:, 1]+x[:, 2])/6) * np.cos(2*np.pi*x[:, 3]/6)
    cols[:, 3] = np.cos(2*np.pi*(x[:, 3]-x[:, 1])/6) * np.sin(2*np.pi*x[:, 4]/5)
    cols[:, 4] = np.sin(2*np.pi*x[:, 4]/5) * np.cos(2*np.pi*x[:, 0]/64)
    cols[:, 5] = np.cos(2*np.pi*(x[:, 0]+x[:, 4])/64) * np.sin(2*np.pi*(x[:, 2]-x[:, 3])/6)
    # C: due righe = cols e una rotazione/shift non monotona
    C = np.zeros((coords.shape[0], 2, N_EFF))
    C[:, 0, :] = cols
    C[:, 1, :] = np.roll(cols, 1, axis=1) * np.cos(np.arange(N_EFF)*0.7)
    return C

def soft_pos(C):
    mins = list(combinations(range(N_EFF), 2))
    pos = np.zeros(C.shape[0])
    for i, j in mins:
        det = C[:,0,i]*C[:,1,j] - C[:,0,j]*C[:,1,i]
        pos += (det >= -1e-12).astype(float)
    return pos / len(mins)

coords_p = prism_coords(N)
C_p = immersion_structured(coords_p)
sp_p = soft_pos(C_p)

# NULL: stesse formule ma coordinate uniformi casuali in [0,64)x[0,6)^3x[0,5)
coords_r = np.column_stack([
    np.random.uniform(0, 64, N),
    np.random.uniform(0, 6, N),
    np.random.uniform(0, 6, N),
    np.random.uniform(0, 6, N),
    np.random.uniform(0, 5, N),
])
C_r = immersion_structured(coords_r)
sp_r = soft_pos(C_r)

print(f"  Soft-pos PRISM (immersione strutturale): {sp_p.mean():.4f} ± {sp_p.std():.4f}")
print(f"  Soft-pos RANDOM (stessa immersione)    : {sp_r.mean():.4f} ± {sp_r.std():.4f}")
print(f"  Δ (PRISM − RANDOM)                     : {sp_p.mean()-sp_r.mean():.4f}")
if abs(sp_p.mean() - sp_r.mean()) < 0.05 and sp_p.mean() < 0.7:
    print("  → Nessuna positività PRISM-specifica rilevata (atteso per immersione generica)")
elif sp_p.mean() > sp_r.mean() + 0.15:
    print("  → Segnale di positività PRISM-specifica (da approfondire)")
else:
    print("  → Differenza debole o nulla; immersione non certifica positroid")

# ============================================================
# 2. QFIM CON DERIVATE LUNGO I CICLI + CONFRONTO OMOGENEO
# ============================================================
print("\n[2] QFIM: derivate lungo i 5 cicli, fit, confronto omogeneo, split strutturale")

N2 = 4096
coords2 = prism_coords(N2)
# Parametri = le 5 coordinate cicliche (derivate finite lungo i cicli)
# Proiettori da Fourier strutturato
th = 2*np.pi*coords2/np.array([64.,6,6,6,5])
def make_P(th):
    v1 = np.stack([np.cos(th[:,0])+1j*np.sin(th[:,1]),
                   np.cos(th[:,2])+1j*np.sin(th[:,3]),
                   np.cos(th[:,4])+1j*np.sin(th[:,0]+th[:,4])],1)
    v2 = np.stack([-np.sin(th[:,0])+1j*np.cos(th[:,1]),
                    np.sin(th[:,2])-1j*np.cos(th[:,3]),
                    np.sin(th[:,4])+1j*np.cos(th[:,0]-th[:,2])],1)
    v1 /= np.linalg.norm(v1,1,keepdims=True)+1e-15
    v2 -= np.sum(v1.conj()*v2,1,keepdims=True)*v1
    v2 /= np.linalg.norm(v2,1,keepdims=True)+1e-15
    P = v1[:,:,None]*v1.conj()[:,None,:] + v2[:,:,None]*v2.conj()[:,None,:]
    return 0.5*(P+np.conjugate(np.transpose(P,(0,2,1))))

P0 = make_P(th)
# Derivate finite lungo ciascun ciclo
eps = 1e-4
dP = []
periods = np.array([64.,6,6,6,5])
for mu in range(5):
    th_p = th.copy()
    th_p[:, mu] += eps * 2*np.pi / periods[mu]
    Pp = make_P(th_p)
    dP.append((Pp - P0)/eps)
dP = np.stack(dP, axis=1)  # (N, 5, 3, 3)

# g_μν = (1/2) Re Tr(∂μP ∂νP)
g = np.zeros((N2, 5, 5))
for mu in range(5):
    for nu in range(mu, 5):
        val = 0.5*np.real(np.trace(dP[:,mu]@dP[:,nu], axis1=1, axis2=2))
        g[:,mu,nu] = g[:,nu,mu] = val
F = 4*g

# embedding z dalle componenti di g e F (senza proiezione random eccessiva)
triu = np.triu_indices(5)
feat = np.hstack([g[:,triu[0],triu[1]], F[:,triu[0],triu[1]]])
# PCA-like: SVD
feat -= feat.mean(0)
U,S,Vt = np.linalg.svd(feat, full_matrices=False)
D_Z = min(10, feat.shape[1])
z = U[:,:D_Z] * S[:D_Z]
z -= z.mean(0)

def fit_Q(z, maxiter=400):
    n,d = z.shape
    def pack(Q): return Q[np.triu_indices(d)]
    def unpack(x):
        Q=np.zeros((d,d)); Q[np.triu_indices(d)]=x
        return Q+Q.T-np.diag(np.diag(Q))
    def loss(x):
        Q=unpack(x); ev,V=eigh(Q)
        Q=V@np.diag(np.maximum(ev,1e-10))@V.T
        return np.mean((np.einsum('ni,ij,nj->n',z,Q,z)-1)**2)
    r = minimize(loss, pack(np.eye(d)), method='L-BFGS-B',
                 options={'maxiter':maxiter,'ftol':1e-15})
    Q=unpack(r.x); ev,V=eigh(Q)
    Q=V@np.diag(np.maximum(ev,1e-10))@V.T
    tr=np.trace(Q)
    if abs(tr)>1e-12: Q*=(d/tr)
    return Q, r.success, r.nit, r.fun

Q_ell, ok, nit, fun = fit_Q(z)
print(f"  Fit ellissoide: success={ok}, nit={nit}, loss={fun:.4e}")

def R_norm(z, Q):
    v = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n',z,Q,z),0))
    v = v/(np.mean(v)+1e-30)
    return np.var(v)/(np.mean(v)**2+1e-30)

Q_I = np.eye(z.shape[1])
R1_I = R_norm(z, Q_I)
R1_ell = R_norm(z, Q_ell)

# Split STRUTTURALE: foglio + = i64 pari, foglio − = i64 dispari (proxy Anti-Twin)
mask_plus = (coords2[:,0] % 2 == 0)
z_plus, z_minus = z[mask_plus], z[~mask_plus]
if len(z_plus)>100 and len(z_minus)>100:
    Qp,_, _,_ = fit_Q(z_plus, maxiter=200)
    Qm,_, _,_ = fit_Q(z_minus, maxiter=200)
    vp = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n',z,Qp,z),0))
    vm = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n',z,Qm,z),0))
    norms = (vp+vm)/(np.mean(vp+vm)+1e-30)
    R1_ov = np.var(norms)/(np.mean(norms)**2+1e-30)
else:
    R1_ov = np.nan

print(f"  R grado-1 sfera      : {R1_I:.6e}")
print(f"  R grado-1 ellissoide : {R1_ell:.6e}")
print(f"  R grado-1 ovaoloide (split strutturale i64 pari/dispari): {R1_ov:.6e}")
if R1_ell <= R1_I * 1.01:
    print("  → Fit ellissoide coerente (R_ell ≤ R_I)")
else:
    print("  → Fit ancora non ottimale (R_ell > R_I)")
if not np.isnan(R1_ov):
    best = np.nanmin([R1_I, R1_ell, R1_ov])
    print(f"  → Migliore a confronto omogeneo: ", end="")
    if best==R1_I: print("sfera")
    elif best==R1_ell: print("ellissoide")
    else: print("ovaoloide strutturale")

# ============================================================
# 3. TORSIONE TIPO WEITZENBÖCK + SPLITTING NON TRIVIALE
# ============================================================
print("\n[3] Torsione discreta tipo Weitzenböck e H non commutante")

# Vierbein discreto grezzo su Q6: e_i(v) = direzione del bit i
# Connessione Weitzenböck: Γ ~ e^{-1} de
# Torsione T ~ Γ - Γ^T (proxy: antisimmetrizzazione del trasporto)

verts = list(np.ndindex(*([2]*6)))
n=64
imap={v:i for i,v in enumerate(verts)}

# Operatore di trasporto lungo direzione d: T_d |v> = |v flip d>
T = np.zeros((6, n, n))
for d in range(6):
    for v in verts:
        w=list(v); w[d]=1-w[d]
        T[d, imap[v], imap[tuple(w)]] = 1.0

# Torsione proxy: commutatore di trasporti [T_a, T_b] (curvatura di plaquette = 0 su ipercubo,
# ma possiamo definire un "twist" pesato)
# Per avere [H, qualcosa] ≠ 0 introduciamo un potenziale di torsione diagonale
# che NON commuta con J in modo banale.

# τ_W(v) = somma_d (-1)^{v_d} * (v·ω)  con ω irrazionale → rompe la simmetria di J
omega = np.array([1.0, np.sqrt(2), np.sqrt(3), np.sqrt(5), np.sqrt(7), np.sqrt(11)])
tau_W = np.array([np.dot(v, omega)*((-1)**sum(v)) for v in verts])
tau_W = (tau_W - tau_W.mean())/(tau_W.std()+1e-12)

L = np.zeros((n,n))
for v in verts:
    i=imap[v]; L[i,i]=6
    for d in range(6):
        w=list(v); w[d]=1-w[d]; L[i,imap[tuple(w)]]=-1

J = np.zeros((n,n))
for v in verts:
    J[imap[v], imap[tuple(1-x for x in v)]] = 1.0

# H = L⊗I + μ (τ_W-weighted J)  ma con τ sul foglio, non solo scala globale
# Costruiamo M = diag(τ_W) @ J  — questo in generale NON dà ΔE=2μ
mu = 0.5
M = np.diag(tau_W) @ J
M = 0.5*(M+M.T)
H = np.block([[L, mu*M],[mu*M, L]])
ev = np.sort(eigh(H, eigvals_only=True))
dE = ev[1]-ev[0]

# Confronto con il caso banale μ J
H_triv = np.block([[L, mu*J],[mu*J, L]])
ev_triv = np.sort(eigh(H_triv, eigvals_only=True))
dE_triv = ev_triv[1]-ev_triv[0]

print(f"  ||[L,J]||_F           = {np.linalg.norm(L@J-J@L):.3e}")
print(f"  ΔE (J banale, μ={mu}) = {dE_triv:.6f}  (atteso 2μ={2*mu})")
print(f"  ΔE (J pesato τ_W)     = {dE:.6f}")
print(f"  |ΔE − 2μ| (pesato)    = {abs(dE-2*mu):.6f}")
if abs(dE - 2*mu) > 0.05:
    print("  → Splitting NON riducibile a 2μ: la torsione modifica lo spettro")
else:
    print("  → Splitting ancora vicino a 2μ; rafforzare la torsione")

# ============================================================
# 4. ADIACENZA REALE PLAQUETTE Q6 + AZIONE DI GAUGE + σ
# ============================================================
print("\n[4] Adiacenza reale plaquette Q6 + Wilson gauge + area law")

# Plaquette di Q6: scegli 2 direzioni distinte (d1,d2), e un vertice di partenza
# con i bit fuori da {d1,d2} fissi. Numero = C(6,2)*2^4 = 15*16 = 240.
plaquettes = []  # lista di 4-cicli (tuple di 4 indici vertice)
for d1, d2 in combinations(range(6), 2):
    fixed_dirs = [d for d in range(6) if d not in (d1,d2)]
    for bits in product([0,1], repeat=4):
        v0 = [0]*6
        for k,d in enumerate(fixed_dirs):
            v0[d] = bits[k]
        # ciclo: v0 --d1--> v1 --d2--> v2 --d1--> v3 --d2--> v0
        def flip(v, d):
            w=list(v); w[d]=1-w[d]; return tuple(w)
        v0t = tuple(v0)
        v1 = flip(v0t, d1)
        v2 = flip(v1, d2)
        v3 = flip(v2, d1)
        # v0 = flip(v3, d2)
        plaquettes.append((imap[v0t], imap[v1], imap[v2], imap[v3]))

assert len(plaquettes)==240
print(f"  Plaquette enumerate: {len(plaquettes)}")

# Adiacenza duale: due plaquette sono adiacenti se condividono un edge
from collections import defaultdict
edge_to_plaq = defaultdict(list)
for p_idx, (a,b,c,d) in enumerate(plaquettes):
    for e in [(a,b),(b,c),(c,d),(d,a)]:
        edge_to_plaq[tuple(sorted(e))].append(p_idx)

adj_plaq = [set() for _ in range(240)]
for plist in edge_to_plaq.values():
    for i in range(len(plist)):
        for j in range(i+1, len(plist)):
            adj_plaq[plist[i]].add(plist[j])
            adj_plaq[plist[j]].add(plist[i])
deg = [len(s) for s in adj_plaq]
print(f"  Grado medio duale (reale): {np.mean(deg):.2f}  (min={min(deg)}, max={max(deg)})")

# Azione di gauge U(1) sulle plaquette: θ_p ∈ [−π,π]
# S = β Σ_p (1−cos θ_p)   (Wilson)
# Stima grezza di σ: in U(1) 2D/4D proxy, σ ≈ −log(<W(C)>)/Area
# Qui misuriamo <cos θ> e usiamo σ_proxy = −log( max(<cos θ>,1e-6) ) come ordine di grandezza

BETA_G = 1.5
N_SWEEP = 300
theta = np.random.uniform(-np.pi, np.pi, 240)
def S_wilson(th, beta):
    return beta * np.sum(1.0 - np.cos(th))

S = S_wilson(theta, BETA_G)
hist_cos = []
for sweep in range(N_SWEEP):
    for i in range(240):
        prop = theta[i] + np.random.normal(0, 0.5)
        prop = (prop+np.pi)%(2*np.pi)-np.pi
        dS = BETA_G * ((1-np.cos(prop))-(1-np.cos(theta[i])))
        if dS<0 or np.random.rand()<np.exp(-dS):
            theta[i] = prop
            S += dS
    if sweep > 80:
        hist_cos.append(np.mean(np.cos(theta)))

hist_cos = np.array(hist_cos)
W = hist_cos.mean()
sigma_proxy = -np.log(max(abs(W), 1e-6))
print(f"  β_gauge = {BETA_G}")
print(f"  <cos θ_p> (Wilson loop 1×1) = {W:.4f} ± {hist_cos.std():.4f}")
print(f"  σ_proxy = −log|<W|>        = {sigma_proxy:.4f}")
print(f"  (Stima grezza; per area law vera servono loop grandi e fit lineare)")

# ============================================================
# RIEPILOGO
# ============================================================
print("\n" + "=" * 70)
print("RIEPILOGO PASSI FALSIFICABILI")
print("=" * 70)
print(f"  1. Plücker strutturale: SP_PRISM={sp_p.mean():.3f}  SP_RAND={sp_r.mean():.3f}")
print(f"  2. QFIM omogeneo: R_I={R1_I:.4e}  R_ell={R1_ell:.4e}  R_ov_strutt={R1_ov:.4e}")
print(f"  3. Splitting con τ_W: ΔE={dE:.4f}  (banale 2μ={2*mu})")
print(f"  4. Plaquette reali=240, grado duale medio={np.mean(deg):.1f}, σ_proxy={sigma_proxy:.3f}")
print(f"  Tempo: {time.time()-t0:.1f}s")
print("=" * 70)

np.savez("/home/workdir/artifacts/falsifiable_1234_results.npz",
         sp_prism=sp_p.mean(), sp_rand=sp_r.mean(),
         R1_I=R1_I, R1_ell=R1_ell, R1_ov=R1_ov,
         dE=dE, dE_triv=dE_triv, sigma_proxy=sigma_proxy,
         deg_mean=np.mean(deg))
print("\nSalvato falsifiable_1234_results.npz")
