#!/usr/bin/env python3
"""
Controlli nulli e correzioni oneste
==================================
1. Plücker sorted: test nullo su dati casuali (deve dare 1.0 → fallisce il test)
2. QFIM: confronto omogeneo grado-1 vs grado-1; verifica convergenza fit
3. J: verifica che ΔE = 2μ è identità algebrica
4. Nessun claim di evidenza dove non c'è
"""

import numpy as np
from itertools import combinations
from scipy.linalg import eigh
from scipy.optimize import minimize
import warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

print("=" * 70)
print("CONTROLLI NULLI E CORREZIONI ONESTE")
print("=" * 70)

# ============================================================
# 1. PLÜCKER SORTED – TEST NULLO
# ============================================================
print("\n[1] Plücker sorted – test nullo (dati casuali vs PRISM)")

N, N_EFF, K = 4096, 6, 2

def sorted_immersion(u6):
    u_sorted = np.sort(u6, axis=1)
    C = np.zeros((u6.shape[0], 2, N_EFF))
    C[:, 0, :] = np.cumsum(u_sorted + 0.05, axis=1)
    C[:, 1, :] = np.cumsum((u_sorted + 0.05)**2, axis=1)
    return C

def soft_pos(C):
    mins = list(combinations(range(N_EFF), 2))
    pos = np.zeros(C.shape[0])
    for i, j in mins:
        det = C[:,0,i]*C[:,1,j] - C[:,0,j]*C[:,1,i]
        pos += (det >= -1e-12).astype(float)
    return pos / len(mins)

# PRISM-like
idx = np.random.choice(69120, N, replace=False)
e=idx%5; tmp=idx//5; c=tmp%6; tmp//=6; b=tmp%6; tmp//=6; a=tmp%6; i64=tmp//6
coords = np.stack([i64,a,b,c,e],1).astype(float)
u = coords/np.array([64.,6,6,6,5])
u6 = np.column_stack([u, (0.5*u[:,0]+0.3*u[:,4])%1])
sp_prism = soft_pos(sorted_immersion(u6))

# NULL: dati completamente casuali in [0,1)
u_rand = np.random.rand(N, N_EFF)
sp_rand = soft_pos(sorted_immersion(u_rand))

print(f"  Soft-pos sorted su PRISM-like : {sp_prism.mean():.6f}")
print(f"  Soft-pos sorted su RANDOM     : {sp_rand.mean():.6f}")
print(f"  Differenza                    : {abs(sp_prism.mean()-sp_rand.mean()):.6e}")
if abs(sp_prism.mean() - sp_rand.mean()) < 1e-6:
    print("  → TEST NULLO FALLITO: la positività è della costruzione, non di PRISM")
else:
    print("  → Differenza rilevata (inaspettato)")

# ============================================================
# 2. QFIM – CONFRONTO OMOGENEO (tutti come norme di grado 1)
# ============================================================
print("\n[2] QFIM – confronto omogeneo grado-1 e verifica fit")

N2, D_PARAM, D_Z = 4096, 6, 8
# proiettori strutturati minimi
th = 2*np.pi*coords[:N2]/np.array([64.,6,6,6,5])
v1 = np.stack([np.cos(th[:,0])+1j*np.sin(th[:,1]),
               np.cos(th[:,2])+1j*np.sin(th[:,3]),
               np.cos(th[:,4])+1j*np.sin(th[:,0]+th[:,4])],1)
v2 = np.stack([-np.sin(th[:,0])+1j*np.cos(th[:,1]),
                np.sin(th[:,2])-1j*np.cos(th[:,3]),
                np.sin(th[:,4])+1j*np.cos(th[:,0]-th[:,2])],1)
v1 /= np.linalg.norm(v1,1,keepdims=True)+1e-15
v2 -= np.sum(v1.conj()*v2,1,keepdims=True)*v1
v2 /= np.linalg.norm(v2,1,keepdims=True)+1e-15
P = 0.5*((v1[:,:,None]*v1.conj()[:,None,:]+v2[:,:,None]*v2.conj()[:,None,:])
     +np.conjugate(np.transpose(v1[:,:,None]*v1.conj()[:,None,:]+v2[:,:,None]*v2.conj()[:,None,:],(0,2,1))))

rng = np.random.RandomState(7)
base = 0.5*(rng.randn(D_PARAM,3,3)+1j*rng.randn(D_PARAM,3,3)
     +np.conjugate(np.transpose(rng.randn(D_PARAM,3,3)+1j*rng.randn(D_PARAM,3,3),(0,2,1))))
g = np.zeros((N2,D_PARAM,D_PARAM))
for mu in range(D_PARAM):
    dmu = P@(0.03*base[mu])+(0.03*base[mu])@P-2*P@(0.03*base[mu])@P
    for nu in range(mu,D_PARAM):
        dnu = P@(0.03*base[nu])+(0.03*base[nu])@P-2*P@(0.03*base[nu])@P
        val = 0.5*np.real(np.trace(dmu@dnu,1,2))
        g[:,mu,nu]=g[:,nu,mu]=val
F=4*g
triu=np.triu_indices(D_PARAM)
feat=np.hstack([g[:,triu[0],triu[1]],F[:,triu[0],triu[1]]])
proj,_=np.linalg.qr(np.random.RandomState(13).randn(feat.shape[1],D_Z))
z=feat@proj; z-=z.mean(0,keepdims=True)

def fit_Q(z, maxiter=200):
    n,d=z.shape
    def pack(Q): return Q[np.triu_indices(d)]
    def unpack(x):
        Q=np.zeros((d,d)); Q[np.triu_indices(d)]=x
        return Q+Q.T-np.diag(np.diag(Q))
    def loss(x):
        Q=unpack(x); ev,V=eigh(Q)
        Q=V@np.diag(np.maximum(ev,1e-8))@V.T
        return np.mean((np.einsum('ni,ij,nj->n',z,Q,z)-1)**2)
    r=minimize(loss, pack(np.eye(d)), method='L-BFGS-B',
               options={'maxiter':maxiter, 'ftol':1e-14})
    Q=unpack(r.x); ev,V=eigh(Q)
    Q=V@np.diag(np.maximum(ev,1e-8))@V.T
    tr=np.trace(Q)
    if abs(tr)>1e-12: Q*=(d/tr)
    return Q, r.success, r.nit, r.fun

Q_ell, ok, nit, fun = fit_Q(z, maxiter=300)
print(f"  Fit ellissoide: success={ok}, nit={nit}, final loss={fun:.6e}")

# Statistiche grado 2
def R_quad(z, Q):
    v = np.einsum('ni,ij,nj->n', z, Q, z)
    return np.var(v)/(np.mean(v)**2+1e-30)

# Statistiche grado 1 (norme)
def R_norm(z, Q):
    v = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n', z, Q, z), 0))
    v = v/(np.mean(v)+1e-30)
    return np.var(v)/(np.mean(v)**2+1e-30)

Q_I = np.eye(D_Z)
# normalizza Tr
Q_I = Q_I * (D_Z/np.trace(Q_I))

R2_I   = R_quad(z, Q_I)
R2_ell = R_quad(z, Q_ell)
R1_I   = R_norm(z, Q_I)
R1_ell = R_norm(z, Q_ell)

# ovaoloide: due fit su split casuale, poi norma somma
ix = np.random.permutation(N2)
Qp,_ ,_,_ = fit_Q(z[ix[:N2//2]], maxiter=150)
Qm,_ ,_,_ = fit_Q(z[ix[N2//2:]], maxiter=150)
vp = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n',z,Qp,z),0))
vm = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n',z,Qm,z),0))
norms = (vp+vm)/(np.mean(vp+vm)+1e-30)
R1_ov = np.var(norms)/(np.mean(norms)**2+1e-30)

print(f"\n  Confronto OMOGENEO:")
print(f"  {'Modello':<14} {'R grado-2':>12} {'R grado-1 (norma)':>18}")
print("  " + "-"*48)
print(f"  {'Sfera':<14} {R2_I:12.4e} {R1_I:18.4e}")
print(f"  {'Ellissoide':<14} {R2_ell:12.4e} {R1_ell:18.4e}")
print(f"  {'Ovaoloide':<14} {'—':>12} {R1_ov:18.4e}")

if R2_ell > R2_I * 1.01:
    print("\n  → R_ell (grado-2) > R_I: il fit NON ha raggiunto un ottimo globale")
    print("    (la sfera è un sottospazio dell'ellissoide; un ottimo non può peggiorare)")
else:
    print("\n  → Fit ellissoide coerente con R_ell ≤ R_I")

print(f"\n  A parità di grado-1: sfera R={R1_I:.4e}, ell R={R1_ell:.4e}, ov R={R1_ov:.4e}")
best = min(R1_I, R1_ell, R1_ov)
if best == R1_I:
    print("  → A confronto omogeneo vince la SFERA (nessun vantaggio ovoidale netto)")
elif best == R1_ov:
    print("  → A confronto omogeneo l'ovaoloide resta migliore")
else:
    print("  → A confronto omogeneo vince l'ellissoide")

# ============================================================
# 3. J – VERIFICA IDENTITÀ ALGEBRICA ΔE = 2μ
# ============================================================
print("\n[3] Operatore J – verifica che ΔE = 2μ è identità algebrica")

verts = list(np.ndindex(*([2]*6)))
n=64
imap={v:i for i,v in enumerate(verts)}
L=np.zeros((n,n))
for v in verts:
    i=imap[v]; L[i,i]=6
    for d in range(6):
        w=list(v); w[d]=1-w[d]; L[i,imap[tuple(w)]]=-1
J=np.zeros((n,n))
for v in verts:
    J[imap[v], imap[tuple(1-x for x in v)]]=1.0

# Verifica [L,J]=0
comm = L@J - J@L
print(f"  ||[L,J]||_F = {np.linalg.norm(comm):.3e}  (deve essere ~0)")

for mu in [0.1, 0.5, 1.0]:
    H = np.block([[L, mu*J],[mu*J, L]])
    ev = np.sort(eigh(H, eigvals_only=True))
    dE = ev[1]-ev[0]
    print(f"  μ={mu:.1f}  ΔE={dE:.6f}  2μ={2*mu:.6f}  |ΔE-2μ|={abs(dE-2*mu):.3e}")

print("  → ΔE = 2μ esattamente (modo costante di J=+1). Non è evidenza dinamica.")

# ============================================================
# SINTESI ONESTA
# ============================================================
print("\n" + "=" * 70)
print("SINTESI ONESTA")
print("=" * 70)
print("""
1. Plücker sorted soft-pos=1 è un artefatto della curva dei momenti
   (totally positive). Fallisce il test nullo su dati casuali.

2. Confronto R_ov vs R_ell era non omogeneo (grado 1 vs 2).
   A parità di statistica il vantaggio ovoidale scompare o si riduce.
   Il fit ellissoide spesso non converge (R_ell > R_I impossibile all'ottimo).

3. Splitting di J è l'identità algebrica ΔE=2μ perché [L,J]=0.
   Non dimostra un Biverso dinamico.

4. Flussoni: grafo duale inventato, formula m_tree non derivata,
   tuning sul target. Nessuna previsione.

STATO REALE:
- Solido: spettro Q6, cardinalità 69120, ψ algebrico.
- Esplorazioni preliminari senza evidenza ancora:
  positività PRISM-specifica, ovaoloide vs sfera, Biverso dinamico, flussoni.
""")
print("=" * 70)
