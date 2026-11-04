#!/usr/bin/env python3
"""RGS-PSU Fase 2: matching diretto vs Riccati, audit, soglia, risonanze.
Modello: -Psi'' + V(u) Psi = w^2 Psi,  V(u) = diag(mu^2) - 2 sech^2(u) I - 2 sech^2(u-s) Coff
(Coff: accoppiamenti alpha (primi vicini), beta (secondi vicini), centrati in u=s: rompe la parita').
Nessuna proiezione/correzione per forzare l'unitarieta'.
"""
import os, sys, time, csv, itertools, copy, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import newton
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
LOGF = open(os.path.join(OUT, "audit_log.txt"), "w")
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()

# ------------------------------------------------------------------ modello
class Model:
    def __init__(self, mu, alpha, beta, shift, flip=False):
        self.mu = np.array(mu, float); self.N = len(mu)
        self.D = np.diag(self.mu**2)
        C = np.zeros((self.N, self.N))
        for i in range(self.N - 1): C[i, i+1] = C[i+1, i] = alpha
        for i in range(self.N - 2): C[i, i+2] = C[i+2, i] = beta
        self.C = C; self.shift = shift; self.flip = flip
    def V(self, u):
        if self.flip: u = -u
        return self.D - 2/np.cosh(u)**2*np.eye(self.N) - 2/np.cosh(u - self.shift)**2*self.C
    def flipped(self):
        m = copy.copy(self); m.flip = not self.flip; return m

def chan(om, mu):
    op = om > mu
    k = np.sqrt(np.maximum(om**2 - mu**2, 0.0)); ka = np.sqrt(np.maximum(mu**2 - om**2, 0.0))
    return op, k, ka

# ------------------------------------------------------------ matching diretto
def fundamental(model, om2, a, b, rtol, atol):
    N = model.N; cplx = np.iscomplexobj(om2)
    def rhs(u, y):
        P = y.reshape(2*N, 2*N)
        A = np.zeros((2*N, 2*N), dtype=y.dtype)
        A[:N, N:] = np.eye(N); A[N:, :N] = model.V(u) - om2*np.eye(N)
        return (A @ P).ravel()
    y0 = np.eye(2*N, dtype=complex if cplx else float).ravel()
    sol = solve_ivp(rhs, (a, b), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:, -1].reshape(2*N, 2*N)

def stk(cols, N):
    return np.array(cols, dtype=complex).reshape(len(cols), 2*N).T

def S_direct(model, om, L, rtol, atol):
    N = model.N; mu = model.mu
    op, k, ka = chan(om, mu); io = np.where(op)[0]; ic = np.where(~op)[0]
    no, nc = len(io), len(ic)
    Phi = fundamental(model, om**2, -L, L, rtol, atol)
    J = np.block([[np.zeros((N, N)), np.eye(N)], [-np.eye(N), np.zeros((N, N))]])
    symp = np.linalg.norm(Phi.T @ J @ Phi - J, 2) / max(1.0, np.linalg.norm(Phi, 2))
    def copen(j, sgn, pos):
        ch = io[j]; v = np.zeros(2*N, complex); ph = np.exp(1j*sgn*k[ch]*pos)/np.sqrt(k[ch])
        v[ch] = ph; v[N+ch] = 1j*sgn*k[ch]*ph; return v
    def cclosed(j, sgn):
        ch = ic[j]; v = np.zeros(2*N, complex); v[ch] = 1.0; v[N+ch] = sgn*ka[ch]; return v
    BLin  = stk([copen(j, +1, -L) for j in range(no)], N)
    BLout = stk([copen(j, -1, -L) for j in range(no)], N)
    BLc   = stk([cclosed(j, +1) for j in range(nc)], N)      # e^{kappa(u+L)}
    BRout = stk([copen(j, +1, +L) for j in range(no)], N)
    BRin  = stk([copen(j, -1, +L) for j in range(no)], N)
    BRc   = stk([cclosed(j, -1) for j in range(nc)], N)      # e^{-kappa(u-L)}
    M = np.hstack([Phi@BLout, Phi@BLc, -BRout, -BRc])
    rhs = np.hstack([-Phi@BLin, BRin])
    cond = np.linalg.cond(M)
    X = np.linalg.solve(M, rhs)
    rL = X[:no, :no]; tR = X[:no, no:]
    tL = X[no+nc:2*no+nc, :no]; rR = X[no+nc:2*no+nc, no:]
    S = np.block([[rL, tR], [tL, rR]])
    return S, dict(cond=cond, symp=symp, no=no)

# ------------------------------------------------------------------- Riccati
def riccati_half(model, om, L, rtol, atol):
    """Incidenza da sinistra: restituisce rL, tL (flusso-normalizzate) e diagnostica."""
    N = model.N; mu = model.mu
    op, k, ka = chan(om, mu); io = np.where(op)[0]; ic = np.where(~op)[0]
    no, nc = len(io), len(ic)
    Y0 = np.diag(np.where(op, 1j*k, -ka)).astype(complex)   # uscente a destra / decadente
    def rY(u, y):
        Y = y.reshape(N, N)
        return (model.V(u) - om**2*np.eye(N) - Y@Y).ravel()
    sy = solve_ivp(rY, (L, -L), Y0.ravel(), method="DOP853", rtol=rtol, atol=atol, dense_output=True)
    if not sy.success: raise RuntimeError("Riccati: " + sy.message)
    Ymax = float(np.max(np.abs(sy.y)))
    Ym = sy.y[:, -1].reshape(N, N)
    Vin = np.zeros((N, no), complex); Din = np.zeros((N, no), complex)
    Vout = np.zeros((N, no), complex); Dout = np.zeros((N, no), complex)
    Vc = np.zeros((N, nc), complex); Dc = np.zeros((N, nc), complex)
    for j, ch in enumerate(io):
        a = np.exp(-1j*k[ch]*L)/np.sqrt(k[ch]); b = np.exp(1j*k[ch]*L)/np.sqrt(k[ch])
        Vin[ch, j] = a; Din[ch, j] = 1j*k[ch]*a
        Vout[ch, j] = b; Dout[ch, j] = -1j*k[ch]*b
    for j, ch in enumerate(ic):
        Vc[ch, j] = 1.0; Dc[ch, j] = ka[ch]
    A = np.hstack([Ym@Vout - Dout, Ym@Vc - Dc]); rhs = -(Ym@Vin - Din)
    condA = np.linalg.cond(A)
    z = np.linalg.solve(A, rhs)
    P0 = Vin + np.hstack([Vout, Vc]) @ z
    def rP(u, p):
        return (sy.sol(u).reshape(N, N) @ p.reshape(N, no)).ravel()
    sp = solve_ivp(rP, (-L, L), P0.ravel().astype(complex), method="DOP853", rtol=rtol, atol=atol)
    if not sp.success: raise RuntimeError("Psi: " + sp.message)
    PL = sp.y[:, -1].reshape(N, no)
    r = z[:no, :]
    t = np.array([PL[io[j], :]*np.sqrt(k[io[j]])*np.exp(-1j*k[io[j]]*L) for j in range(no)])
    return r, t, dict(Ymax=Ymax, cond=condA, no=no)

def S_riccati(model, om, L, rtol, atol):
    rL, tL, d1 = riccati_half(model, om, L, rtol, atol)
    rR, tR, d2 = riccati_half(model.flipped(), om, L, rtol, atol)
    S = np.block([[rL, tR], [tL, rR]])
    return S, dict(cond=max(d1["cond"], d2["cond"]), Ymax=max(d1["Ymax"], d2["Ymax"]), no=d1["no"])

# ------------------------------------------------------------------ audit
def metrics(S):
    n = S.shape[0]
    unit = np.linalg.norm(S.conj().T@S - np.eye(n), 2)
    rec = np.linalg.norm(S - S.T, 2)
    flux_in = np.abs(np.sum(np.abs(S)**2, axis=0) - 1).max()     # per canale incidente
    flux_out = np.abs(np.sum(np.abs(S)**2, axis=1) - 1).max()
    return unit, rec, flux_in, flux_out

def eig_tracks(Ss):
    th = [np.sort(np.angle(np.linalg.eigvals(S))) for S in Ss]
    n = len(th[0]); raw = np.zeros((len(th), n)); raw[0] = th[0]
    for i in range(1, len(th)):
        best = None
        for p in itertools.permutations(range(n)):
            d = np.sum(np.abs(np.angle(np.exp(1j*(th[i][list(p)] - raw[i-1])))))
            if best is None or d < best[0]: best = (d, list(p))
        raw[i] = th[i][best[1]]
    return np.unwrap(raw, axis=0)

def detphase(Ss):
    return np.unwrap(np.array([np.angle(np.linalg.det(S)) for S in Ss]))

# ================================================================== MAIN
t_start = time.time()
RT, AT = 1e-11, 1e-13
MU3 = (1.0, 1.5, 2.0)
ALPHA, BETA, SHIFT = 0.2, 0.1, 0.8
M3 = Model(MU3, ALPHA, BETA, SHIFT)
log("# RGS-PSU fase 2 -- avvio", time.ctime())
log(f"Modello 3 canali: mu={MU3}, alpha={ALPHA}, beta={BETA}, shift accoppiamento={SHIFT}")
log("V(u)=diag(mu^2)-2sech^2(u) I-2sech^2(u-s)*Coff; V(inf)=diag(1,2.25,4)")

# ---------------------------------------------------- 0. validazione codice
log("\n## 0. Validazioni del codice")
# 0a modello storico 2 canali a w=1.6
M2 = Model((1.0, 1.5), 0.2, 0.0, 0.0)
om = 1.6; k = np.sqrt(om**2 - M2.mu**2); K = np.diag(k)
Sd, dd = S_direct(M2, om, 12, RT, AT); Sr, dr = S_riccati(M2, om, 12, RT, AT)
def unnorm(S, k):
    no = len(k); r = S[:no, :no]; t = S[no:, :no]
    Kh, Kmh = np.diag(np.sqrt(k)), np.diag(1/np.sqrt(k))
    return Kmh@r@Kh, Kmh@t@Kh
R_un, T_un = unnorm(Sd, k)
log("2 canali, w=1.6, L=12 -- R non normalizzata (matching):\n", np.array2string(R_un, precision=6))
log("T non normalizzata (matching):\n", np.array2string(T_un, precision=6))
R_doc = np.array([[0.003214+0.015416j, 0.039634+0.005600j], [0.088909+0.012563j, -0.015348+0.009249j]])
T_doc = np.array([[0.203309+0.943670j, -0.167201-0.026689j], [-0.375077-0.059870j, -0.488848+0.832347j]])
log("max|R - R_doc| =", np.abs(R_un - R_doc).max(), " max|T - T_doc| =", np.abs(T_un - T_doc).max())
res_mine = np.linalg.norm(R_un.conj().T@K@R_un + T_un.conj().T@K@T_un - K, 2)
res_doc = np.linalg.norm(R_doc.conj().T@K@R_doc + T_doc.conj().T@K@T_doc - K, 2)
log(f"residuo flusso pesato: calcolato qui = {res_mine:.3e};  con le cifre del testo (6 decimali) = {res_doc:.3e}")
Rf = np.diag(np.sqrt(k))@R_un@np.diag(1/np.sqrt(k)); Tf = np.diag(np.sqrt(k))@T_un@np.diag(1/np.sqrt(k))
log(f"reciprocita' blocchi flusso-norm.: ||R-R^T||={np.linalg.norm(Rf-Rf.T,2):.3e}  ||T-T^T||={np.linalg.norm(Tf-Tf.T,2):.3e}")
log("S globale 2 canali (matching): unit=%.3e recip=%.3e | Riccati: unit=%.3e recip=%.3e | ||Sd-Sr||=%.3e" %
    (*metrics(Sd)[:2], *metrics(Sr)[:2], np.linalg.norm(Sd-Sr, 2)))
# 0b disaccoppiato: sech^2 con l=1 e' riflessionless
M3dec = Model(MU3, 0.0, 0.0, 0.0)
rmax = 0
for om in (1.2, 1.4, 1.6, 1.8):
    S, _ = S_riccati(M3dec, om, 16, RT, AT); no = S.shape[0]//2
    rmax = max(rmax, np.abs(S[:no, :no]).max())
log(f"disaccoppiato (Poschl-Teller l=1, atteso riflessionless): max|r| = {rmax:.3e}")
VALID_DEC = rmax

# ---------------------------------------------------- 1. scansione
log("\n## 1. Scansione in omega (3 canali)")
g = np.concatenate([np.linspace(1.05, 1.95, 181),
                    1.118 + np.linspace(-0.03, 0.03, 61), 1.732 + np.linspace(-0.03, 0.03, 61),
                    1.5 + np.concatenate([-np.geomspace(2e-5, 0.1, 30), np.geomspace(2e-5, 0.1, 30)])])
g = np.unique(np.round(g, 9)); g = g[np.abs(g - 1.5) > 1e-9]
log(f"{len(g)} frequenze; direct+Riccati a L=8, Riccati a L=16 (riferimento); rtol={RT}, atol={AT}")
rows = []; Sd8, Sr8, Sr16 = [], [], []
fails = []
for i, om in enumerate(g):
    try:
        a, da = S_direct(M3, om, 8, RT, AT)
        b, db = S_riccati(M3, om, 8, RT, AT)
        c, dc = S_riccati(M3, om, 16, RT, AT)
    except Exception as e:
        fails.append((om, str(e))); log("ERRORE a w=%.9f: %s" % (om, e)); Sd8.append(None); Sr8.append(None); Sr16.append(None); continue
    Sd8.append(a); Sr8.append(b); Sr16.append(c)
    no = da["no"]; ma = metrics(a); mb = metrics(b); mc = metrics(c)
    row = dict(omega=om, n_open=no, k1=np.sqrt(max(om**2-1, 0)), k2=np.sqrt(max(om**2-2.25, 0)),
               unit_dir8=ma[0], recip_dir8=ma[1], fluxin_dir8=ma[2], fluxout_dir8=ma[3], cond_dir8=da["cond"], symp_dir8=da["symp"],
               unit_ric8=mb[0], recip_ric8=mb[1], fluxin_ric8=mb[2], fluxout_ric8=mb[3], cond_ric8=db["cond"], Ymax_ric8=db["Ymax"],
               unit_ric16=mc[0], recip_ric16=mc[1], fluxin_ric16=mc[2], fluxout_ric16=mc[3], cond_ric16=dc["cond"], Ymax_ric16=dc["Ymax"],
               dS_dir8_ric8=np.linalg.norm(a-b, 2), dS_ric8_ric16=np.linalg.norm(b-c, 2), dS_dir8_ric16=np.linalg.norm(a-c, 2))
    for j in range(no):
        row[f"T_L{j+1}"] = np.sum(np.abs(c[no:, j])**2); row[f"R_L{j+1}"] = np.sum(np.abs(c[:no, j])**2)
    for p in range(2*no):
        for q in range(2*no):
            row[f"S16_{p}{q}_re"] = c[p, q].real; row[f"S16_{p}{q}_im"] = c[p, q].imag
    rows.append(row)
    if i % 40 == 0: log(f"  [{i}/{len(g)}] w={om:.5f} unit_ric16={mc[0]:.2e} unit_dir8={ma[0]:.2e} dS(d8,r8)={row['dS_dir8_ric8']:.2e} cond_dir8={da['cond']:.1e}")
keys = []
for r in rows:
    for kk in r:
        if kk not in keys: keys.append(kk)
with open(os.path.join(OUT, "results_scan.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow(r)
okidx = [i for i in range(len(g)) if Sr16[i] is not None]
oms = g[okidx]; S16 = [Sr16[i] for i in okidx]; SD8 = [Sd8[i] for i in okidx]; SR8 = [Sr8[i] for i in okidx]
nopen = np.array([S.shape[0]//2 for S in S16])
def mx(key, mask=None):
    v = np.array([r[key] for r in rows]); return v.max() if mask is None else v[mask].max()
log("\nRiassunto scansione:")
for key in ["unit_dir8", "unit_ric8", "unit_ric16", "recip_dir8", "recip_ric8", "recip_ric16", "fluxin_ric16", "fluxout_ric16",
            "dS_dir8_ric8", "dS_ric8_ric16", "dS_dir8_ric16", "cond_dir8", "cond_ric16", "Ymax_ric16", "symp_dir8"]:
    log(f"  max {key:16s} = {mx(key):.3e}")
log("  errori/eccezioni:", len(fails))

# eigenfasi e fase del determinante per segmento
seg = {}
for name, nn in (("A_sotto_soglia(1 canale aperto)", 1), ("B_sopra_soglia(2 canali aperti)", 2)):
    idx = np.where(nopen == nn)[0]
    seg[name] = (oms[idx], eig_tracks([S16[i] for i in idx]), detphase([S16[i] for i in idx]))
    log(f"segmento {name}: {len(idx)} punti, eigenfasi tracciate con continuita' di ramo")

# indicatore di risoluzione in energia
Tarr = np.array([r["T_L1"] for r in rows])
curv = np.abs(Tarr[1:-1] - 0.5*(Tarr[:-2] + Tarr[2:]))
# ---------------------------------------------------- 2. convergenza
log("\n## 2. Tabella di convergenza (dominio x tolleranza)")
wl = [1.12, 1.4, 1.499, 1.501, 1.6, 1.74, 1.8]
Ls = [10, 12, 16, 20]; tols = [(1e-8, 1e-10), (1e-10, 1e-12), (1e-12, 1e-14)]
conv = []
for om in wl:
    ref, _ = S_riccati(M3, om, 20, 1e-12, 1e-14)
    for L in Ls:
        for rt, at in tols:
            rec = dict(omega=om, L=L, rtol=rt, atol=at)
            for nm, fn in (("dir", S_direct), ("ric", S_riccati)):
                try:
                    S, d = fn(M3, om, L, rt, at); mt = metrics(S)
                    rec.update({f"unit_{nm}": mt[0], f"recip_{nm}": mt[1], f"fluxin_{nm}": mt[2], f"cond_{nm}": d["cond"],
                                f"dref_{nm}": np.linalg.norm(S-ref, 2)})
                except Exception as e:
                    rec.update({f"unit_{nm}": np.nan, f"recip_{nm}": np.nan, f"fluxin_{nm}": np.nan, f"cond_{nm}": np.nan, f"dref_{nm}": np.nan})
                    log(f"ERRORE conv {nm} w={om} L={L} tol={rt}: {e}")
            try:
                S1, _ = S_direct(M3, om, L, rt, at); S2, _ = S_riccati(M3, om, L, rt, at); rec["dS_dir_ric"] = np.linalg.norm(S1-S2, 2)
            except Exception:
                rec["dS_dir_ric"] = np.nan
            conv.append(rec)
    log(f"  w={om} fatto")
ck = list(conv[0].keys())
with open(os.path.join(OUT, "convergence_table.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=ck); w.writeheader(); [w.writerow(r) for r in conv]
log("(riferimento: Riccati L=20 rtol=1e-12 atol=1e-14)")
for om in wl:
    sub = [r for r in conv if r["omega"] == om and r["rtol"] == 1e-12]
    log(f" w={om}: " + " | ".join(f"L={r['L']}: ric dref={r['dref_ric']:.1e}, dir dref={r['dref_dir']:.1e} (cond {r['cond_dir']:.0e})" for r in sub))
conv_ric = max(r["dref_ric"] for r in conv if r["L"] >= 12 and r["rtol"] <= 1e-10)
log(f"max ||S_ric(L>=12, rtol<=1e-10) - S_ref|| = {conv_ric:.3e}")

# ---------------------------------------------------- 3. soglia, cusp, tunnelling
log("\n## 3. Soglia w=1.5 e limite k2->0")
cusp = {}
for L in (12, 16, 20):
    d = np.geomspace(1e-4, 5e-3, 16)
    up = []; dn = []
    for x in d:
        S, _ = S_riccati(M3, 1.5 + x, L, RT, AT); up.append(np.sum(np.abs(S[2:, 0])**2))
        S, _ = S_riccati(M3, 1.5 - x, L, RT, AT); dn.append(np.sum(np.abs(S[1:, 0])**2))
    k2 = np.sqrt((1.5+d)**2 - 2.25)
    cu, cov_u = np.polyfit(k2, up, 2, cov=True)
    cd, cov_d = np.polyfit(-d, dn, 2, cov=True)
    su = np.sqrt(np.diag(cov_u)); sd = np.sqrt(np.diag(cov_d))
    cusp[L] = dict(c0=cu[2], c0s=su[2], c1=cu[1], c1s=su[1], a0=cd[2], a0s=sd[2], up=np.array(up), dn=np.array(dn), d=d, k2=k2)
    log(f" L={L}: sopra P=c0+c1*k2+c2*k2^2: c0={cu[2]:.8f}+-{su[2]:.1e}, c1={cu[1]:.5f}+-{su[1]:.1e}; sotto P(w=1.5-)={cd[2]:.8f}+-{sd[2]:.1e}")
c1s_ = [cusp[L]["c1"] for L in (12, 16, 20)]
c1_stab = (max(c1s_) - min(c1s_)) / max(abs(np.mean(c1s_)), 1e-300)
sig = abs(cusp[20]["c1"]) / cusp[20]["c1s"]
cont = abs(cusp[20]["c0"] - cusp[20]["a0"]) / np.hypot(cusp[20]["c0s"], cusp[20]["a0s"])
log(f" significativita' c1 = {sig:.1f} sigma; variazione relativa di c1 fra L=12,16,20 = {c1_stab:.2e}; |c0-a0|/sigma = {cont:.2f}")
CUSP_OK = (sig > 10) and (c1_stab < 0.1) and (cont < 3)
log(" criterio cuspide (c1>10 sigma, stabilita' <10%, continuita' in soglia <3 sigma):", "SODDISFATTO" if CUSP_OK else "NON soddisfatto")

# tunnelling vs disaccoppiato: disaccoppiato e' riflessionless (T=1)
wsup = np.array([r["omega"] for r in rows]); Tt = np.array([r["T_L1"] for r in rows])
supp = 1 - Tt
log(f" soppressione di trasmissione canale 1 (1-T rispetto al disaccoppiato T=1): max={supp.max():.4f} a w={wsup[np.argmax(supp)]:.4f}, min={supp.min():.2e}")
Tdec = []
for om in (1.2, 1.4, 1.6, 1.8):
    S, _ = S_riccati(M3dec, om, 16, RT, AT); no = S.shape[0]//2; Tdec.append(np.sum(np.abs(S[no:, 0])**2))
log(f" controllo numerico disaccoppiato T(canale 1) a w=1.2,1.4,1.6,1.8: {np.array2string(np.array(Tdec), precision=10)}")

# ---------------------------------------------------- 4. poli
log("\n## 4. Poli di S(w): det(Y_R(0)-Y_L(0)) con Riccati da entrambi i lati")
def sheet_k(omc, mu):
    s = np.where(omc.real > mu, -1.0, 1.0)
    return s*1j*np.sqrt(mu**2 - omc**2 + 0j)
def Fjost(omc, model, L, rtol=1e-11, atol=1e-13):
    N = model.N; k = sheet_k(complex(omc), model.mu); om2 = complex(omc)**2
    def rY(u, y):
        Y = y.reshape(N, N); return (model.V(u) - om2*np.eye(N) - Y@Y).ravel()
    sR = solve_ivp(rY, (L, 0), np.diag(1j*k).ravel(), method="DOP853", rtol=rtol, atol=atol)
    sL = solve_ivp(rY, (-L, 0), np.diag(-1j*k).ravel(), method="DOP853", rtol=rtol, atol=atol)
    return np.linalg.det(sR.y[:, -1].reshape(N, N) - sL.y[:, -1].reshape(N, N))
def Fjost_fund(omc, model, L):
    N = model.N; k = sheet_k(complex(omc), model.mu)
    Phi = fundamental(model, complex(omc)**2, -L, L, 1e-11, 1e-13)
    BL = np.vstack([np.eye(N), -1j*np.diag(k)]); BR = np.vstack([np.eye(N), 1j*np.diag(k)])
    return np.linalg.det(np.hstack([Phi@BL, BR]))
poles = []
guesses = [complex(a, b) for a in np.arange(1.03, 1.96, 0.03) for b in (-0.003, -0.015, -0.05)]
guesses = [z for z in guesses if abs(z.real - 1.5) > 0.006 and abs(z.real - 1.0) > 0.006]
t0 = time.time()
for z0 in guesses:
    try:
        z = newton(lambda x: Fjost(x, M3, 10), z0, x1=z0 + 1e-4 - 1e-4j, tol=1e-13, maxiter=60, disp=False)
        if isinstance(z, complex) or np.iscomplexobj(z):
            z = complex(z)
            if 1.0 < z.real < 1.98 and abs(z.imag) < 0.5 and abs(z.real - 1.5) > 0.004:
                if not any(abs(z - p) < 1e-6 for p in poles): poles.append(z)
    except Exception:
        pass
log(f" ricerca poli: {len(guesses)} punti di partenza, {time.time()-t0:.0f}s, candidati grezzi: {len(poles)}")
ptab = []
for z in sorted(poles, key=lambda x: x.real):
    try:
        zs = {}
        for L in (8, 10, 12):
            zs[L] = complex(newton(lambda x: Fjost(x, M3, L), z, x1=z + 1e-4, tol=1e-13, maxiter=60, disp=False))
        zf = complex(newton(lambda x: Fjost_fund(x, M3, 6), z, x1=z + 1e-4, tol=1e-12, maxiter=60, disp=False))
        spread = max(abs(zs[a] - zs[b]) for a in zs for b in zs)
        lhs = abs(Fjost(zs[10], M3, 10)); rhs_ = abs(Fjost(zs[10] + 1e-3, M3, 10))
        sheet_ok = (sheet_k(zs[10], M3.mu).real > -1e-12).all() or True
        # fase del determinante di S attraverso Re(w)
        wr, wi = zs[10].real, zs[10].imag; dphi = np.nan
        for name, (ww, trk, dph) in seg.items():
            if ww.min() <= wr <= ww.max():
                half = max(3*abs(wi), 0.005); lo, hi = max(wr-half, ww.min()), min(wr+half, ww.max())
                dphi = (np.interp(hi, ww, dph) - np.interp(lo, ww, dph))/np.pi
        ptab.append(dict(re=zs[10].real, im=zs[10].imag, spread_L8_10_12=spread, d_fund_L6=abs(zf - zs[10]), absF_root=lhs,
                         absF_offroot=rhs_, ratio=rhs_/max(lhs, 1e-300), dphase_over_pi=dphi, width=-2*zs[10].imag))
    except Exception as e:
        log("  polo scartato:", z, e)
# decoupled reference positions
log(" posizioni attese del disaccoppiato (stato legato canale j: w^2 = mu_j^2 - 1): ch2 -> 1.1180, ch3 -> 1.7321 (nel continuo dei canali piu' bassi)")
with open(os.path.join(OUT, "poles.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["re_w", "im_w", "width_Gamma=-2Im", "spread_L8_10_12", "diff_vs_fundamental_L6", "absF_root", "absF_off", "ratio", "dphase_over_pi"])
    for p in ptab: w.writerow([p["re"], p["im"], p["width"], p["spread_L8_10_12"], p["d_fund_L6"], p["absF_root"], p["absF_offroot"], p["ratio"], p["dphase_over_pi"]])
for p in ptab:
    log(f"  w = {p['re']:.8f} {p['im']:+.3e}i | Gamma={p['width']:.3e} | spread(L=8,10,12)={p['spread_L8_10_12']:.1e} | vs fund.L=6: {p['d_fund_L6']:.1e} | |F|off/|F|root={p['ratio']:.1e} | d(arg det S)/pi attraverso Re w = {p['dphase_over_pi']:.3f}")
res_ok = [p for p in ptab if p["spread_L8_10_12"] < 1e-6 and p["ratio"] > 1e4 and p["im"] < -1e-12]
log(f" poli che soddisfano (spread<1e-6, ratio>1e4, Im<0): {len(res_ok)} su {len(ptab)}")

# ---------------------------------------------------- 5. grafici
log("\n## 5. Grafici")
cs = ["C0", "C1", "C2", "C3"]
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
for j, lab in enumerate(["T_L1", "R_L1", "T_L2", "R_L2"]):
    v = np.array([r.get(lab, np.nan) for r in rows]); ax.flat[0].plot(wsup, v, ".", ms=3, label=lab)
ax.flat[0].axvline(1.5, c="k", ls=":"); ax.flat[0].set_title("Trasmissioni/riflessioni totali per canale incidente (Riccati L=16)"); ax.flat[0].legend(); ax.flat[0].set_xlabel("w")
for (p, q) in [(0, 0), (1, 0), (2, 0), (3, 0)]:
    v = np.array([np.hypot(r.get(f"S16_{p}{q}_re", np.nan), r.get(f"S16_{p}{q}_im", np.nan)) for r in rows])
    ax.flat[1].plot(wsup, v, ".", ms=3, label=f"|S[{p},{q}]|")
ax.flat[1].axvline(1.5, c="k", ls=":"); ax.flat[1].set_title("Ampiezze |S_p0| (incidenza canale 1 da sinistra)"); ax.flat[1].legend(); ax.flat[1].set_xlabel("w")
ax.flat[2].semilogy(wsup, supp.clip(1e-16), ".", ms=3); ax.flat[2].axvline(1.5, c="k", ls=":"); ax.flat[2].set_title("1 - T(canale 1) rispetto al disaccoppiato (T=1)"); ax.flat[2].set_xlabel("w")
for name, key, ls_ in (("sotto", "dn", "o"), ("sopra", "up", "s")):
    pass
cc = cusp[20]
ax.flat[3].plot(cc["k2"], cc["up"], "s-", label="sopra: P vs k2"); ax.flat[3].plot(-cc["d"], cc["dn"], "o-", label="sotto: P vs (w-1.5)")
ax.flat[3].axvline(0, c="k", ls=":"); ax.flat[3].set_title("Soglia: P_tot(canale 1) (L=20)"); ax.flat[3].legend(); ax.flat[3].set_xlabel("k2 (sopra) / w-1.5 (sotto)")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_amplitudes.png"), dpi=130); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a_, (name, (ww, trk, dph)) in zip(ax, seg.items()):
    for j in range(trk.shape[1]): a_.plot(ww, trk[:, j], ".-", ms=3, label=f"ramo {j}")
    a_.plot(ww, dph/ trk.shape[1] if False else dph, "k--", lw=0.8, label="arg det S")
    for p in ptab:
        if ww.min() <= p["re"] <= ww.max(): a_.axvline(p["re"], c="r", ls=":")
    a_.set_title(name + " -- eigenfasi (rad)"); a_.set_xlabel("w"); a_.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_eigenphases.png"), dpi=130); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for key in ["unit_dir8", "unit_ric8", "unit_ric16", "recip_ric16", "dS_dir8_ric8", "dS_ric8_ric16"]:
    ax[0].semilogy(wsup, np.array([r[key] for r in rows]).clip(1e-18), ".", ms=3, label=key)
ax[0].axvline(1.5, c="k", ls=":"); ax[0].legend(fontsize=7); ax[0].set_title("Residui (unitarieta', reciprocita', differenze fra metodi)"); ax[0].set_xlabel("w")
ax[1].semilogy(wsup, [r["cond_dir8"] for r in rows], ".", ms=3, label="cond M (matching L=8)")
ax[1].semilogy(wsup, [r["cond_ric16"] for r in rows], ".", ms=3, label="cond A (Riccati L=16)")
ax[1].semilogy(wsup, [r["Ymax_ric16"] for r in rows], ".", ms=3, label="max|Y| Riccati L=16"); ax[1].legend(fontsize=7); ax[1].set_title("Condizionamento / polo Riccati"); ax[1].set_xlabel("w")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_residuals.png"), dpi=130); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for om in (1.4, 1.499, 1.6, 1.8):
    for nm, ls_ in (("dir", "--"), ("ric", "-")):
        xs = [r for r in conv if r["omega"] == om and r["rtol"] == 1e-12]
        ax[0].semilogy([r["L"] for r in xs], [max(r[f"dref_{nm}"], 1e-18) for r in xs], ls_, marker="o", label=f"{nm} w={om}")
ax[0].set_xlabel("L"); ax[0].set_title("||S(L)-S_ref|| (rtol=1e-12)"); ax[0].legend(fontsize=6)
for om in (1.4, 1.6):
    for nm, ls_ in (("dir", "--"), ("ric", "-")):
        xs = [r for r in conv if r["omega"] == om and r["L"] == 12]
        ax[1].loglog([r["rtol"] for r in xs], [max(r[f"unit_{nm}"], 1e-18) for r in xs], ls_, marker="o", label=f"{nm} w={om}")
ax[1].set_xlabel("rtol"); ax[1].set_title("residuo unitarieta' vs tolleranza (L=12)"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_convergence.png"), dpi=130); plt.close()

# ---------------------------------------------------- 6. rapporto
unit_ric16 = mx("unit_ric16"); unit_dir8 = mx("unit_dir8"); dS = mx("dS_dir8_ric8"); dtr = mx("dS_ric8_ric16"); rec16 = mx("recip_ric16")
cl = lambda c: "VERIFICATO NUMERICAMENTE" if c else "NON DIMOSTRATO"
rep = []
rep.append("# RGS-PSU Fase 2 -- Rapporto finale (generato dallo script, numeri reali)\n")
rep.append(f"Modello: 3 canali, mu=(1,1.5,2), alpha={ALPHA}, beta={BETA}, accoppiamenti centrati in u={SHIFT}; DOP853 rtol={RT}, atol={AT}; {len(g)} frequenze in [1.05,1.95].\n")
rep.append("## Validazioni del codice")
rep.append(f"- disaccoppiato (PT l=1, riflessionless): max|r| = {VALID_DEC:.2e}  -> {cl(VALID_DEC < 1e-8)}")
rep.append(f"- 2 canali storico w=1.6: residuo flusso pesato (calcolato) = {res_mine:.2e}; con le cifre a 6 decimali del testo = {res_doc:.2e} (il valore 1.42e-10 citato nel testo e' incompatibile con dati a 6 cifre)\n")
rep.append("## Audit matching diretto vs Riccati (scansione completa)")
rep.append(f"- unitarieta' max: Riccati L=16 {unit_ric16:.2e}; matching L=8 {unit_dir8:.2e}")
rep.append(f"- reciprocita' max (Riccati L=16): {rec16:.2e}")
rep.append(f"- max ||S_dir(L=8)-S_ric(L=8)|| = {dS:.2e};  max ||S_ric(8)-S_ric(16)|| (troncamento) = {dtr:.2e}")
rep.append(f"- max cond(M) matching L=8 = {mx('cond_dir8'):.2e}; max|Y| Riccati L=16 = {mx('Ymax_ric16'):.2e}; eccezioni = {len(fails)}")
rep.append(f"- convergenza Riccati (L>=12, rtol<=1e-10): max ||S-S_ref|| = {conv_ric:.2e}")
rep.append(f"- unitarieta' S (Riccati, L=16, tutta la scansione): {cl(unit_ric16 < 1e-7)} (soglia 1e-7)")
rep.append(f"- reciprocita' S: {cl(rec16 < 1e-7)} (soglia 1e-7)")
rep.append(f"- accordo matching/Riccati a L=8: {cl(dS < 1e-5)} (soglia 1e-5); il matching diretto degrada a L grande per canali chiusi (vedi cond e convergence_table.csv)\n")
rep.append("## Soglia w=1.5")
rep.append(f"- cuspide di Wigner (P_tot canale 1): c1={cusp[20]['c1']:.4g} ({sig:.1f} sigma), stabilita' fra L=12/16/20 {c1_stab:.1e}, continuita' {cont:.2f} sigma -> " + ("CANDIDATO/VERIFICATO secondo criterio" if CUSP_OK else "NON DIMOSTRATO"))
rep.append(f"- soppressione trasmissione canale 1 rispetto al disaccoppiato: max(1-T)={supp.max():.3f}; il disaccoppiato e' T=1 esatto, quindi e' soppressione da accoppiamento (riflessione + conversione), non tunnelling in barriera.\n")
rep.append("## Poli di S(w)")
for p in ptab:
    st = "CANDIDATO (convergente in L, F~0, fase det S)" if (p in res_ok and abs(p["dphase_over_pi"]) > 0.3) else "NON DIMOSTRATO"
    rep.append(f"- w = {p['re']:.8f} {p['im']:+.3e}i  Gamma={p['width']:.3e}  spread_L={p['spread_L8_10_12']:.1e}  dphase/pi={p['dphase_over_pi']:.3f} -> {st}")
rep.append("\n## Non dimostrato / limiti")
rep.append("- BIC: nessun polo reale nel continuo e' stato cercato oltre la ricerca generale; un polo con Im=0 non e' stato trovato/classificato come BIC.")
rep.append("- La scansione e' limitata a w<=1.95 (terza soglia a w=2 non esplorata); sheet dei poli fissata dal segno Re w vs mu_j (poli vicini a soglia, |Re w-mu_j|<0.004-0.006, esclusi).")
rep.append("- Il confronto con il caso disaccoppiato usa che PT l=1 e' riflessionless; il modello e' scelto da me (accoppiamenti sech^2 centrati in s=0.8) e non coincide per costruzione col potenziale RGS-PSU originale a 3 canali, che non esiste nel materiale recuperato.")
open(os.path.join(OUT, "report.md"), "w").write("\n".join(rep) + "\n")
log("\n".join(rep))
log(f"\nTempo totale {time.time()-t_start:.0f}s")
LOGF.close()
