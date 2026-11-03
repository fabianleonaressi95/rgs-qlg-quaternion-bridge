#!/usr/bin/env python3
"""
Certificazione QFIM / ovaoloide – script funzionante
===================================================
Genera proiettori strutturati dal prodotto C64□C6³□C5,
calcola g e F_Q, fitta ellissoide e ovaoloide-stadium,
riporta R_ell, R_ov, R_I e i residuali.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
import time, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

N = 8192          # campione (per velocità; scaling a 69120 è lineare)
D_PARAM, D_Z = 8, 12
print("=" * 60)
print(f"QFIM / OVOID CERTIFICATION  (N={N})")
print("=" * 60)
t0 = time.time()

# --- coordinate stratificate ---
idx = np.random.choice(69120, N, replace=False)
e = idx % 5; tmp = idx // 5
c = tmp % 6; tmp //= 6
b = tmp % 6; tmp //= 6
a = tmp % 6; i64 = tmp // 6
coords = np.stack([i64, a, b, c, e], 1).astype(float)

# --- proiettori rango 2 ---
th = 2 * np.pi * coords / np.array([64., 6, 6, 6, 5])
v1 = np.stack([np.cos(th[:,0])+1j*np.sin(th[:,1]),
               np.cos(th[:,2])+1j*np.sin(th[:,3]),
               np.cos(th[:,4])+1j*np.sin(th[:,0]+th[:,4])], 1)
v2 = np.stack([-np.sin(th[:,0])+1j*np.cos(th[:,1]),
                np.sin(th[:,2])-1j*np.cos(th[:,3]),
                np.sin(th[:,4])+1j*np.cos(th[:,0]-th[:,2])], 1)
v1 /= np.linalg.norm(v1, axis=1, keepdims=True) + 1e-15
v2 -= np.sum(v1.conj()*v2, axis=1, keepdims=True) * v1
v2 /= np.linalg.norm(v2, axis=1, keepdims=True) + 1e-15
P = (v1[:,:,None]*v1.conj()[:,None,:] + v2[:,:,None]*v2.conj()[:,None,:])
P = 0.5 * (P + np.conjugate(np.transpose(P, (0,2,1))))
print(f"Proiettori: Tr(P[0])={np.trace(P[0]).real:.4f}")

# --- metriche g, F_Q ---
rng = np.random.RandomState(7)
base = rng.randn(D_PARAM, 3, 3) + 1j * rng.randn(D_PARAM, 3, 3)
base = 0.5 * (base + np.conjugate(np.transpose(base, (0,2,1))))
g = np.zeros((N, D_PARAM, D_PARAM))
for mu in range(D_PARAM):
    dmu = P@(0.03*base[mu]) + (0.03*base[mu])@P - 2*P@(0.03*base[mu])@P
    for nu in range(mu, D_PARAM):
        dnu = P@(0.03*base[nu]) + (0.03*base[nu])@P - 2*P@(0.03*base[nu])@P
        val = 0.5 * np.real(np.trace(dmu @ dnu, axis1=1, axis2=2))
        g[:, mu, nu] = g[:, nu, mu] = val
F = 4.0 * g

# --- embedding z ---
triu = np.triu_indices(D_PARAM)
feat = np.hstack([g[:, triu[0], triu[1]], F[:, triu[0], triu[1]]])
proj, _ = np.linalg.qr(np.random.RandomState(13).randn(feat.shape[1], D_Z))
z = feat @ proj
z -= z.mean(0, keepdims=True)

# --- fit ellissoide ---
def fit_Q(z, maxiter=50):
    n, d = z.shape
    def pack(Q): return Q[np.triu_indices(d)]
    def unpack(x):
        Q = np.zeros((d, d)); Q[np.triu_indices(d)] = x
        return Q + Q.T - np.diag(np.diag(Q))
    def loss(x):
        Q = unpack(x)
        ev, V = eigh(Q)
        Q = V @ np.diag(np.maximum(ev, 1e-8)) @ V.T
        return np.mean((np.einsum('ni,ij,nj->n', z, Q, z) - 1.0)**2)
    res = minimize(loss, pack(np.eye(d)), method='L-BFGS-B', options={'maxiter': maxiter})
    Q = unpack(res.x)
    ev, V = eigh(Q)
    Q = V @ np.diag(np.maximum(ev, 1e-8)) @ V.T
    tr = np.trace(Q)
    if abs(tr) > 1e-12:
        Q *= (d / tr)
    return Q

Q_ell = fit_Q(z)
vals_ell = np.einsum('ni,ij,nj->n', z, Q_ell, z)
mean_ell = np.mean(vals_ell)
R_ell = np.var(vals_ell) / (mean_ell**2 + 1e-30)
E_ell = np.mean((vals_ell - 1)**2)

# --- ovaoloide (split) ---
ix = np.random.permutation(N)
Qp = fit_Q(z[ix[:N//2]])
Qm = fit_Q(z[ix[N//2:]])
vp = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n', z, Qp, z), 0))
vm = np.sqrt(np.maximum(np.einsum('ni,ij,nj->n', z, Qm, z), 0))
norms = vp + vm
norms = norms / (np.mean(norms) + 1e-30)
R_ov = np.var(norms) / (np.mean(norms)**2 + 1e-30)
E_ov = np.mean((norms - 1)**2)

# --- sfera ---
vals_I = np.sum(z**2, axis=1)
vals_I = vals_I / (np.mean(vals_I) + 1e-30)
R_I = np.var(vals_I) / (np.mean(vals_I)**2 + 1e-30)

print(f"\nRisultati (N={N}):")
print(f"  R_ell = {R_ell:.6e}   E_ell = {E_ell:.6e}")
print(f"  R_ov  = {R_ov:.6e}   E_ov  = {E_ov:.6e}")
print(f"  R_I   = {R_I:.6e}")
print(f"  R_ov / R_ell = {R_ov/R_ell:.4f}")
print(f"  Tempo = {time.time()-t0:.1f}s")

if R_ov < R_ell:
    print("\n  R_ov < R_ell  →  ovaoloide preferito su questo campione")
else:
    print("\n  R_ov ≥ R_ell  →  ellissoide competitivo o fit da rivedere")

np.savez("/home/workdir/artifacts/qfim_working_results.npz",
         R_ell=R_ell, R_ov=R_ov, R_I=R_I, E_ell=E_ell, E_ov=E_ov, N=N)
print("\nSalvato qfim_working_results.npz")
print("=" * 60)
