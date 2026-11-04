#!/usr/bin/env python3
"""Monte Carlo SU(2) Wilson 4D — versione auditata."""
import argparse, csv, json, time
from pathlib import Path
import numpy as np
try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

def qmul(a, b):
    a0, av = a[..., :1], a[..., 1:]
    b0, bv = b[..., :1], b[..., 1:]
    scalar = a0*b0 - np.sum(av*bv, axis=-1, keepdims=True)
    vector = a0*bv + b0*av + np.cross(av, bv, axis=-1)
    return np.concatenate((scalar, vector), axis=-1)

def qdag(a):
    out = a.copy(); out[..., 1:] *= -1.0; return out

def qnorm(a): return np.sqrt(np.sum(a*a, axis=-1))

def qidentity(shape):
    out = np.zeros(tuple(shape)+(4,)); out[..., 0] = 1.0; return out

def shift_by(field, offsets):
    return np.roll(field, shift=tuple(-int(v) for v in offsets), axis=(0,1,2,3))

def shift_one(field, mu, amount=1):
    o=[0,0,0,0]; o[mu]=amount; return shift_by(field, o)

def make_parity_masks(L):
    c = np.indices((L,L,L,L)); p = np.sum(c, axis=0)%2
    return p==0, p==1

def random_hot_start(L, rng):
    U = rng.normal(size=(4,L,L,L,L,4))
    return U/np.linalg.norm(U, axis=-1, keepdims=True)

def plaquette_field(U, mu, nu):
    p = qmul(U[mu], shift_one(U[nu], mu, +1))
    p = qmul(p, qdag(shift_one(U[mu], nu, +1)))
    return qmul(p, qdag(U[nu]))

def average_plaquette(U):
    return float(np.mean([np.mean(plaquette_field(U,mu,nu)[...,0])
        for mu in range(4) for nu in range(mu+1,4)]))

def spatial_scalar_operator(U):
    L=U.shape[1]; O=np.zeros(L); c=0
    for mu in range(3):
        for nu in range(mu+1,3):
            O += np.mean(plaquette_field(U,mu,nu)[...,0], axis=(0,1,2)); c+=1
    return O/c

def staple(U, mu):
    K = np.zeros_like(U[mu])
    for nu in range(4):
        if nu==mu: continue
        forward = qmul(shift_one(U[nu], mu, +1), qdag(shift_one(U[mu], nu, +1)))
        forward = qmul(forward, qdag(U[nu]))
        om=[0,0,0,0]; om[nu]=-1
        omn=[0,0,0,0]; omn[mu]=+1; omn[nu]=-1
        backward = qmul(qdag(shift_by(U[nu], omn)), qdag(shift_by(U[mu], om)))
        backward = qmul(backward, shift_by(U[nu], om))
        K += forward + backward
    return K

def sample_a0(kappa, rng):
    kappa=np.asarray(kappa,float); n=len(kappa)
    res=np.empty(n); pend=np.ones(n,bool)
    while np.any(pend):
        idx=np.flatnonzero(pend)
        a=2*rng.random(len(idx))-1
        w=np.sqrt(np.maximum(0,1-a*a))*np.exp(np.minimum(0,kappa[idx]*(a-1)))
        acc=rng.random(len(idx))<w
        g=idx[acc]; res[g]=a[acc]; pend[g]=False
    return res

def random_su2_boltzmann(kappa, rng):
    kappa=np.asarray(kappa,float); a0=sample_a0(kappa,rng)
    d=rng.normal(size=(len(kappa),3)); d/=np.linalg.norm(d,axis=1,keepdims=True)
    r=np.sqrt(np.maximum(0,1-a0*a0))
    return np.column_stack((a0, r[:,None]*d))

def heatbath_link(U, mu, mask, beta, rng):
    K=staple(U,mu); Ks=K[mask]; n=qnorm(Ks)
    v=n>1e-14
    if not np.any(v): return
    Ks=Ks[v]; nv=n[v]; q=random_su2_boltzmann(beta*nv, rng)
    new=qmul(q, qdag(Ks)/nv[:,None])
    old=U[mu][mask]; old[v]=new; U[mu][mask]=old

def overrelax_link(U, mu, mask):
    K=staple(U,mu); Ks=K[mask]; n=qnorm(Ks)
    v=n>1e-14
    if not np.any(v): return
    Us=U[mu][mask]
    q=qmul(Us[v], Ks[v]/n[v,None])
    qr=q.copy(); qr[...,1:]*=-1.0
    Us[v]=qmul(qr, qdag(Ks[v])/n[v,None])
    U[mu][mask]=Us

def sweep(U, beta, masks, rng, nor=2):
    for mask in masks:
        for mu in range(4): heatbath_link(U,mu,mask,beta,rng)
    for _ in range(nor):
        for mask in masks:
            for mu in range(4): overrelax_link(U,mu,mask)

def wilson_rectangle(U,R,T,spatial_mu=0,time_mu=3):
    L=U.shape[1]; P=qidentity((L,L,L,L))
    for k in range(R):
        o=[0,0,0,0]; o[spatial_mu]=k; P=qmul(P,shift_by(U[spatial_mu],o))
    for k in range(T):
        o=[0,0,0,0]; o[spatial_mu]=R; o[time_mu]=k; P=qmul(P,shift_by(U[time_mu],o))
    for k in range(R-1,-1,-1):
        o=[0,0,0,0]; o[spatial_mu]=k; o[time_mu]=T; P=qmul(P,qdag(shift_by(U[spatial_mu],o)))
    for k in range(T-1,-1,-1):
        o=[0,0,0,0]; o[time_mu]=k; P=qmul(P,qdag(shift_by(U[time_mu],o)))
    return float(np.mean(P[...,0]))

def measure_loops(U, rmax=3, tmax=3):
    L=U.shape[1]; res={}
    for R in range(1,min(rmax,L//2)+1):
        for T in range(1,min(tmax,L//2)+1):
            res[f"W{R}_{T}"]=float(np.mean([wilson_rectangle(U,R,T,spatial_mu=m) for m in range(3)]))
    return res

def creutz_ratios(lm, rmax=3, tmax=3):
    out={}
    for R in range(2,rmax+1):
        for T in range(2,tmax+1):
            ks=[f"W{R}_{T}",f"W{R-1}_{T-1}",f"W{R-1}_{T}",f"W{R}_{T-1}"]
            if not all(k in lm for k in ks): continue
            a,b,c,d=[lm[k] for k in ks]
            out[f"chi{R}_{T}"]=np.nan if min(a,b,c,d)<=0 else float(-np.log(a*b/(c*d)))
    return out

def scalar_correlator(op):
    x=np.asarray(op,float); c=x-x.mean(); Lt=x.shape[1]
    return np.array([np.mean(c*np.roll(c,-t,axis=1)) for t in range(Lt)])

def cosh_model(t,A,m,Lt): return A*(np.exp(-m*t)+np.exp(-m*(Lt-t)))

def fit_glueball_mass(corr,Lt,tmin=1,tmax=None):
    if not HAVE_SCIPY: return {"mass_lattice":None,"fit_status":"SciPy non installato"}
    if tmax is None: tmax=Lt//2
    ts=np.arange(max(1,tmin),min(tmax,Lt//2)); ys=np.asarray(corr[ts])
    ok=np.isfinite(ys)&(ys>0); ts,ys=ts[ok],ys[ok]
    if len(ts)<3: return {"mass_lattice":None,"fit_status":"Punti positivi insufficienti","fit_range":[int(tmin),int(tmax)]}
    try:
        popt,pcov=curve_fit(lambda t,A,m:cosh_model(t,A,m,Lt),ts,ys,p0=(float(ys[0]),0.5),
                            bounds=([0,1e-5],[np.inf,10]),maxfev=20000)
        err=float(np.sqrt(pcov[1,1])) if pcov.shape==(2,2) and pcov[1,1]>=0 else None
        return {"mass_lattice":float(popt[1]),"mass_fit_error_formal":err,
                "fit_status":"ok; verificare finestra","fit_range":[int(ts[0]),int(ts[-1])],"n":int(len(ts))}
    except Exception as e:
        return {"mass_lattice":None,"fit_status":f"fit fallito: {e}"}

def jackknife_mean(values):
    x=np.asarray(values,float); x=x[np.isfinite(x)]
    if len(x)<2: return {"mean":float(np.mean(x)) if len(x) else None,"error":None,"n":int(len(x))}
    n=len(x); m=float(np.mean(x)); loo=(np.sum(x)-x)/(n-1)
    return {"mean":m,"error":float(np.sqrt((n-1)/n*np.sum((loo-m)**2))),"n":int(n)}

def autocorr_time(values, max_lag=100):
    x=np.asarray(values,float); n=len(x)
    if n<8: return None
    x=x-x.mean(); var=np.dot(x,x)/n
    if var<=0: return 0.5
    tau=0.5
    for lag in range(1,min(max_lag,n//4)):
        rho=np.dot(x[:-lag],x[lag:])/((n-lag)*var)
        if rho<=0: break
        tau+=rho
    return float(tau)

def run(args):
    L=args.L; rng=np.random.default_rng(args.seed); masks=make_parity_masks(L)
    U = random_hot_start(L,rng) if args.hot else qidentity((L,L,L,L))
    print(f"SU(2) 4D beta={args.beta} L={L} seed={args.seed} hot={args.hot}", flush=True)
    t0=time.time()
    for i in range(args.therm): sweep(U,args.beta,masks,rng,args.nor)
    plaq=[]; loops=[]; ops=[]
    rmax=min(args.rmax,L//2); tmax=min(args.tmax,L//2)
    for i in range(args.nmeas):
        for _ in range(args.interval): sweep(U,args.beta,masks,rng,args.nor)
        plaq.append(average_plaquette(U))
        loops.append(measure_loops(U,rmax,tmax))
        ops.append(spatial_scalar_operator(U))
        if (i+1)%max(1,args.nmeas//10)==0:
            print(f"meas {i+1}/{args.nmeas} plaq={plaq[-1]:.6f} t={time.time()-t0:.0f}s", flush=True)
    plaq=np.asarray(plaq); ops=np.asarray(ops)
    keys=sorted(loops[0]); lm=np.array([[r[k] for k in keys] for r in loops])
    lmean={k:float(np.nanmean(lm[:,j])) for j,k in enumerate(keys)}
    lstats={k:jackknife_mean(lm[:,j]) for j,k in enumerate(keys)}
    creutz=creutz_ratios(lmean,rmax,tmax)
    corr=scalar_correlator(ops)
    mf=fit_glueball_mass(corr,Lt=L,tmin=args.fit_tmin,tmax=args.fit_tmax)
    sig=[v for k,v in creutz.items() if k.startswith("chi") and np.isfinite(v)]
    sigma=float(np.mean(sig)) if sig else None
    mass=mf.get("mass_lattice")
    ratio=float(mass/np.sqrt(sigma)) if mass is not None and sigma and sigma>0 else None
    out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,beta=np.array(args.beta),L=np.array(L),plaquette=plaq,
        loop_keys=np.array(keys),loop_samples=lm,operator_samples=ops,
        scalar_correlator=corr,creutz_json=np.array(json.dumps(creutz)))
    rep={"beta":args.beta,"L":L,"seed":args.seed,"hot":args.hot,
         "plaquette":jackknife_mean(plaq),"wilson_loops":lstats,"creutz_ratios":creutz,
         "sigma_lattice_candidate":sigma,"glueball_0pp":mf,
         "m0pp_over_sqrt_sigma_candidate":ratio,
         "tau_int_plaquette_estimate":autocorr_time(plaq),
         "runtime_seconds":time.time()-t0}
    out.with_suffix(".json").write_text(json.dumps(rep,indent=2))
    print("DONE", json.dumps({k:rep[k] for k in ["plaquette","sigma_lattice_candidate","glueball_0pp","tau_int_plaquette_estimate","runtime_seconds"]}), flush=True)

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="command",required=True)
    r=sub.add_parser("run")
    r.add_argument("--beta",type=float,required=True); r.add_argument("--L",type=int,default=8)
    r.add_argument("--therm",type=int,default=200); r.add_argument("--nmeas",type=int,default=200)
    r.add_argument("--interval",type=int,default=1); r.add_argument("--nor",type=int,default=2)
    r.add_argument("--seed",type=int,default=12345); r.add_argument("--out",type=str,default="mc.npz")
    r.add_argument("--rmax",type=int,default=3); r.add_argument("--tmax",type=int,default=3)
    r.add_argument("--fit-tmin",type=int,default=1); r.add_argument("--fit-tmax",type=int,default=None)
    r.add_argument("--hot",action="store_true")
    args=p.parse_args(); run(args)

if __name__=="__main__": main()
