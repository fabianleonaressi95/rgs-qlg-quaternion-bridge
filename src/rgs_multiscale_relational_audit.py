"""Finite RGS–PSU–Q6 multiscale/relational-time audit.

This is a reproducible finite-model bridge, not a derivation of quantum gravity.
Dependencies: Python >=3.10 and NumPy.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
from typing import Any
import numpy as np

SEED = 20261104
ATOL = 1.0e-11
RTOL = 1.0e-10

def qn_adjacency(n: int) -> np.ndarray:
    if n < 1: raise ValueError("n must be positive")
    a = np.zeros((1 << n, 1 << n), dtype=float)
    for x in range(1 << n):
        for bit in range(n):
            y = x ^ (1 << bit)
            a[x, y] = 1.0
    return a

def _sha256_array(a: np.ndarray) -> str:
    x = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(str(x.shape).encode()); h.update(x.dtype.str.encode()); h.update(x.tobytes())
    return h.hexdigest()

def q6_audit() -> dict[str, Any]:
    a = qn_adjacency(6); lap = 6*np.eye(64)-a
    got = np.linalg.eigvalsh(lap)
    expected = np.repeat(2*np.arange(7, dtype=float), [math.comb(6,j) for j in range(7)])
    residual = float(np.max(np.abs(got-expected)))
    return {"vertices":64, "undirected_edges":int(np.count_nonzero(np.triu(a,1))),
            "degree_min":int(a.sum(1).min()), "degree_max":int(a.sum(1).max()),
            "betti_1_graph":129, "laplacian_spectrum_max_residual":residual,
            "adjacency_sha256":_sha256_array(a),
            "pass":bool(residual<ATOL and np.all(a.sum(1)==6))}

def relational_clock_audit() -> dict[str, Any]:
    ac, ass = qn_adjacency(2), qn_adjacency(4)
    ec, vc = np.linalg.eigh(ac); es, vs = np.linalg.eigh(ass)
    pairs=[]; used=set()
    for i,e in enumerate(ec):
        candidates=[j for j,x in enumerate(es) if j not in used and abs(x+e)<ATOL]
        if not candidates: raise AssertionError(f"no Q4 partner for energy {e}")
        j=candidates[0]; used.add(j); pairs.append((i,j))
    coeff=np.zeros((4,16),complex)
    for i,j in pairs: coeff[i,j]=1/math.sqrt(len(pairs))
    residual=float(np.linalg.norm(np.diag(ec)@coeff+coeff@np.diag(es),ord="fro"))
    def conditional(theta):
        return vs @ (coeff.T @ np.exp(1j*ec*theta))
    p0,p1=conditional(0),conditional(.37)
    overlap=abs(np.vdot(p0,p1))/(np.linalg.norm(p0)*np.linalg.norm(p1))
    distinguishability=float(math.sqrt(max(0,1-min(1,float(overlap))**2)))
    return {"clock_factor":"Q2","system_factor":"Q4","clock_dimension":4,"system_dimension":16,
            "clock_distinct_energies":sorted(set(np.round(ec,12).tolist())),
            "global_state_terms":len(pairs),"global_constraint_residual":residual,
            "conditional_state_distinguishability_theta_0_vs_0p37":distinguishability,
            "clock_period":math.pi,
            "interpretation_limit":"finite relational-clock toy; theta is not calibrated proper time; only three distinct frequencies",
            "pass":bool(residual<ATOL and distinguishability>1e-6)}

def _rates_q6(driven: bool, affinity: float=1.2) -> np.ndarray:
    w=qn_adjacency(6)
    if driven:
        fwd,rev=math.exp(affinity/2),math.exp(-affinity/2)
        for i,j in [(0,1),(1,3),(3,2),(2,0)]:
            w[i,j],w[j,i]=fwd,rev
    return w

def markov_entropy_audit(seed: int=SEED, trajectories: int=2500, horizon: float=3.0) -> dict[str,Any]:
    rng=np.random.default_rng(seed); pi=np.full(64,1/64); out={}
    for label,driven in [("equilibrium",False),("driven_cycle",True)]:
        w=_rates_q6(driven); flux=pi[:,None]*w; current=flux-flux.T
        lr=np.zeros_like(w)
        for i in range(64):
            for j in range(64):
                if w[i,j]>0: lr[i,j]=math.log(w[i,j]/w[j,i])
        ep=.5*float(np.sum(current*lr))
        stationarity=float(np.max(np.abs(current.sum(0))))
        db=float(np.max(np.abs(current)))
        sigmas=[]; jumps=0
        for _ in range(trajectories):
            state=int(rng.integers(64)); t=0.; sigma=0.
            while True:
                rates=w[state].copy(); rates[state]=0.; total=float(rates.sum())
                dt=float(rng.exponential(1/total))
                if t+dt>horizon: break
                nxt=int(rng.choice(64,p=rates/total))
                sigma+=math.log(w[state,nxt]/w[nxt,state]); state=nxt; t+=dt; jumps+=1
            sigmas.append(sigma)
        out[label]={"stationarity_residual_uniform_pi":stationarity,
            "detailed_balance_max_flux_residual":db,
            "entropy_production_rate_nats_per_model_time":ep,
            "sample_mean_path_entropy_nats":float(np.mean(sigmas)),
            "sample_path_entropy_min":float(np.min(sigmas)),
            "sample_path_entropy_max":float(np.max(sigmas)),
            "trajectories":trajectories,"horizon_model_units":horizon,"total_jumps":jumps,
            "pass":bool(stationarity<ATOL and (ep<ATOL if not driven else ep>0))}
    out.update({"seed":seed,"model_note":"Driven square affinity is prescribed, not derived from PSU geometry.",
                "pass":bool(out["equilibrium"]["pass"] and out["driven_cycle"]["pass"])})
    return out

def acoustic_metric_audit(cs:float=1.,radial_speed:float=0.,circulation:float=1.,radius:float=2.) -> dict[str,Any]:
    if cs<=0 or radius<=0 or radial_speed<0: raise ValueError("invalid acoustic parameters")
    vr,vt=-radial_speed/radius,circulation/radius
    g=np.array([[-cs**2+vr**2+vt**2,-vr,-vt],[-vr,1.,0.],[-vt,0.,1.]])
    eig=np.linalg.eigvalsh(g); det=float(np.linalg.det(g)); resid=abs(det+cs**2)
    return {"metric_2plus1":g.tolist(),"metric_determinant":det,
        "determinant_identity_residual_abs_det_plus_cs2":resid,"metric_eigenvalues":eig.tolist(),
        "lorentzian_signature_one_negative_two_positive":bool(np.sum(eig<0)==1 and np.sum(eig>0)==2),
        "acoustic_horizon_radius":radial_speed/cs if radial_speed>0 else None,
        "ergoradius":math.sqrt(radial_speed**2+circulation**2)/cs,
        "parameters_dimensionless_illustration":{"sound_speed":cs,"draining_speed_parameter":radial_speed,
            "circulation_parameter":circulation,"evaluation_radius":radius},
        "physical_scope":"analogue acoustic metric, not Einstein gravity; pure circulation has no radial horizon",
        "pass":bool(resid<ATOL and np.sum(eig<0)==1 and np.sum(eig>0)==2)}

def discrete_scale_toy_audit(delta:float=math.log(2)) -> dict[str,Any]:
    if delta<=0: raise ValueError("delta must be positive")
    x=np.geomspace(1e-5,1e2,2000); scale=math.exp(delta)
    f=np.sin(2*math.pi*np.log(x)/delta)
    shifted=np.sin(2*math.pi*np.log(scale*x)/delta)
    resid=float(np.max(np.abs(f-shifted)))
    return {"log_echoing_period":delta,"scale_factor":scale,"max_echo_residual":resid,
        "status":"TOY_IDENTITY_ONLY_NOT_A_REPRODUCTION_OF_ECKER_ET_AL","pass":resid<1e-10}

def gravity_scale_audit() -> dict[str,Any]:
    G,hbar,c=6.67430e-11,1.054571817e-34,299792458.
    mp=math.sqrt(hbar*c/G)
    masses={"electron":9.1093837139e-31,"proton":1.67262192595e-27,"planck_mass":mp}
    rows={}
    for name,m in masses.items():
        alpha=G*m*m/(hbar*c); r=hbar/(m*c); compact=2*G*m/(r*c*c)
        rows[name]={"mass_kg":m,"gravitational_coupling_alpha_g":alpha,
            "reduced_compton_length_m":r,"compactness_at_reduced_compton_length":compact,
            "compactness_identity_residual_abs_C_minus_2alpha":abs(compact-2*alpha)}
    return {"planck_mass_kg":mp,"species":rows,
        "interpretation":"Particle self-gravity is tiny in this measure; Planck coupling is order one. This does not explain temporal or thermodynamic fluctuations.",
        "pass":bool(all(x["compactness_identity_residual_abs_C_minus_2alpha"]<1e-12 for x in rows.values())
            and rows["electron"]["gravitational_coupling_alpha_g"]<1e-40
            and abs(rows["planck_mass"]["gravitational_coupling_alpha_g"]-1)<1e-12)}

def run_audit(seed:int=SEED) -> dict[str,Any]:
    q6=q6_audit(); clock=relational_clock_audit(); markov=markov_entropy_audit(seed)
    acoustic=acoustic_metric_audit(); dsi=discrete_scale_toy_audit(); gravity=gravity_scale_audit()
    sections={"q6_exact_graph":q6,"relational_clock_q2_q4":clock,
        "markov_fluctuation_controls":markov,"superfluid_analogue_metric":acoustic,
        "discrete_scale_toy":dsi,"gravity_scale_diagnostic":gravity}
    return {"audit":"RGS_PSU_Q6_MULTISCALE_RELATIONAL_TIME_V1",
        "classification":"finite-model audit; no claim of a unified physical theory",
        "seed":seed,"tolerances":{"absolute":ATOL,"relative":RTOL},
        "input_hashes":{"q6_adjacency_sha256":q6["adjacency_sha256"]},
        "sections":sections,"pass":all(v["pass"] for v in sections.values()),
        "open_physics_questions":["derive a physical clock observable rather than choose a clock factor",
            "derive Markov rates from an RGS-PSU action or data",
            "derive a continuum vierbein/metric from the discrete state space",
            "reproduce the large-D Einstein-Klein-Gordon solution and finite-D limits",
            "calibrate analogue-fluid metric without identifying it with real gravity",
            "derive scale-dependent gravitational/thermodynamic coupling rather than fit it"]}

def main():
    p=argparse.ArgumentParser(); p.add_argument("--output",default="rgs_multiscale_relational_audit.json")
    p.add_argument("--seed",type=int,default=SEED); args=p.parse_args()
    report=run_audit(args.seed)
    Path(args.output).write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"pass":report["pass"],"output":args.output,
        "section_pass":{k:v["pass"] for k,v in report["sections"].items()}},indent=2))
    if not report["pass"]: raise SystemExit(1)

if __name__=="__main__": main()
