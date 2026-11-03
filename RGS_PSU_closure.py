#!/usr/bin/env python3
"""RGS-PSU multi-background diagnostic suite.

This is deliberately a DIAGNOSTIC profile evaluator, not a converged BVP solver.
The three historical initial-data sets are kept separate. Geometric quantities
are computed for ds^2=-N(u)^2 dt^2+du^2+R(u)^2 dOmega^2. Teleparallel
identities are NOT marked physically verified until T and the torsion vector
are derived consistently from a specified tetrad and inertial spin connection.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp

OUT = Path(__file__).resolve().parent

CONFIGS = [
    {"id": "config_1_percorso_A", "label": "Percorso A non-triviale",
     "R0": 2.31, "phi0": 0.80, "N0": 1.0, "source": "reported Percorso A baseline"},
    {"id": "config_2_formulazione_4D", "label": "FORMULAZIONE 4D background",
     "R0": 0.8292, "phi0": 0.3450, "N0": 1.4529, "source": "reported FORMULAZIONE 4D data"},
    {"id": "config_3_siiiiii_unified", "label": "siiiiii unified framework",
     "R0": 2.0, "phi0": 0.3450, "N0": 0.9821, "source": "reported siiiiii background data"},
]


def deriv(y, u, order=1):
    """Numerical derivative for quantities without an analytic derivative."""
    out = np.asarray(y, dtype=float)
    for _ in range(order):
        out = np.gradient(out, u, edge_order=2)
    return out


def profile_exact(u, b, phi0, N0, scalar_width):
    """Return the analytic diagnostic ansatz and exact derivatives through order 2."""
    R = np.sqrt(u*u + b*b)
    Rp = u / R
    Rpp = b*b / R**3
    w2 = scalar_width**2
    phi = phi0*np.exp(-u*u/w2)
    phip = -2*u/w2*phi
    phipp = (4*u*u/w2**2 - 2/w2)*phi
    N = np.full_like(u, N0, dtype=float)
    return {"R": R, "Rp": Rp, "Rpp": Rpp, "phi": phi,
            "phip": phip, "phipp": phipp, "N": N,
            "Np": np.zeros_like(u), "Npp": np.zeros_like(u)}


class GeometricIdentities:
    """Metric geometry for proper-radial coordinate u."""
    @staticmethod
    def metric_scalars(N, R, u):
        Rp, Rpp = deriv(R, u), deriv(R, u, 2)
        Np, Npp = deriv(N, u), deriv(N, u, 2)
        # Ricci scalar for signature (-,+,+,+), standard Riemann convention.
        Rlc = -2*Npp/N - 4*Rpp/R - 4*Np*Rp/(N*R) - 2*(Rp/R)**2 + 2/R**2
        # Mixed Einstein tensor components.
        Gtt = 2*Rpp/R + (Rp/R)**2 - 1/R**2
        Guu = 2*Np*Rp/(N*R) + (Rp/R)**2 - 1/R**2
        Gth = Npp/N + Np*Rp/(N*R) + Rpp/R
        return {"R_LC": Rlc, "G_t_t": Gtt, "G_u_u": Guu,
                "G_theta_theta": Gth, "R_prime": Rp, "R_second": Rpp,
                "N_prime": Np, "N_second": Npp}

    @staticmethod
    def report_unverified_teleparallel_identity():
        return {
            "status": "NOT_VERIFIED",
            "max_abs_residual": None,
            "rms_residual": None,
            "reason": ("A metric alone does not fix the teleparallel torsion scalar: "
                       "specify tetrad, inertial spin connection, torsion-vector convention, "
                       "and sign conventions, then derive all terms from that same choice."),
            "warning": "Do not manufacture PASS by defining one term from the identity."
        }


class PhysicalQuantities:
    @staticmethod
    def misner_sharp_mass(R, u, G=1.0):
        Rp = deriv(R, u)
        return R/(2*G) * (1.0 - Rp**2)

    @staticmethod
    def scalar_potential(phi, V0=0.1, lam=1.0, phi_v=0.5):
        V = V0 + (lam/4.0)*(phi**2 - phi_v**2)**2
        Vp = lam*phi*(phi**2 - phi_v**2)
        Vpp = lam*(3*phi**2 - phi_v**2)
        return V, Vp, Vpp

    @staticmethod
    def em_density(R, Q=1.0, ell0=1.0):
        # This is the supplied phenomenological density; no NLED action is inferred.
        return Q**2/(2.0*np.sqrt(R**4 + ell0**4))


def l2(y, u):
    return float(np.sqrt(np.trapezoid(np.asarray(y)**2, u)))


def evaluate_config(cfg, u, params):
    # Analytic diagnostic ansatz only; no claim that it solves field equations.
    prof = profile_exact(u, cfg["R0"], cfg["phi0"], cfg["N0"], params["scalar_width"])
    R, Rp_exact, Rpp_exact = prof["R"], prof["Rp"], prof["Rpp"]
    N, phi = prof["N"], prof["phi"]
    geo = GeometricIdentities.metric_scalars(N, R, u)
    # Exact derivatives for this analytic ansatz avoid finite-difference noise.
    geo["R_prime"], geo["R_second"] = Rp_exact, Rpp_exact
    geo["R_LC"] = -4*Rpp_exact/R - 2*(Rp_exact/R)**2 + 2/R**2
    geo["G_t_t"] = 2*Rpp_exact/R + (Rp_exact/R)**2 - 1/R**2
    geo["G_u_u"] = (Rp_exact/R)**2 - 1/R**2
    geo["G_theta_theta"] = Rpp_exact/R

    # Scalar residual of the displayed Percorso-A form, with N constant and
    # the source's proposed torsion expressions explicitly treated as unverified.
    # This is a diagnostic residual, not a variationally certified field equation.
    T_claim = -2*Rpp_exact/R - 2*(Rp_exact/R)**2 - 4/R**2
    divT_claim = -2*Rpp_exact/R - 2*(Rp_exact/R)**2 - 1/R**2
    phip = prof["phip"]
    # Since N is constant, □phi = phi'' + 2(R'/R)phi'.
    boxphi = prof["phipp"] + 2*(Rp_exact/R)*phip
    _, Vp, _ = PhysicalQuantities.scalar_potential(phi, params["V0"], params["lam"], params["phi_v"])
    Fp = 2*params["xi"]*phi
    # ∇(F'(phi) T^u) cannot be certified without a consistently derived T^u;
    # this diagnostic proxy is labeled accordingly.
    scalar_proxy = boxphi - Vp - 0.5*Fp*T_claim + params["xi"]*phi*divT_claim

    mms = PhysicalQuantities.misner_sharp_mass(R, u, params["G"])
    rho = PhysicalQuantities.em_density(R, params["Q"], params["ell0"])
    # Full two-sided proper-volume integral over the finite diagnostic domain.
    Eem = float(np.trapezoid(4*np.pi*rho*R**2, u))
    center = int(np.argmin(np.abs(u)))
    # For this ansatz m_MS -> 0 as R -> infinity. This is not an ADM validation,
    # and the supplied rho_EM ~ Q^2/(2R^2) makes its integrated energy diverge
    # linearly with domain size. Therefore no ADM mass is reported as physical.
    return {
        "id": cfg["id"], "label": cfg["label"], "source": cfg["source"],
        "status": "PROFILE_DIAGNOSTIC_NOT_BVP_SOLUTION",
        "method": "analytic_profile_with_residuals",
        "parameters_at_throat": {"R0_target": cfg["R0"], "phi0_target": cfg["phi0"], "N0": cfg["N0"]},
        "measured_at_u_near_zero": {"u": float(u[center]), "R0": float(R[center]),
            "phi0": float(phi[center]), "N0": float(N[center]),
            "R_prime0": float(Rp_exact[center]), "phi_prime0": float(phip[center])},
        "geometry": {
            "R_LC_at_throat": float(geo["R_LC"][center]),
            "G_t_t_at_throat": float(geo["G_t_t"][center]),
            "G_u_u_at_throat": float(geo["G_u_u"][center]),
            "G_theta_theta_at_throat": float(geo["G_theta_theta"][center]),
            "teleparallel_identity": GeometricIdentities.report_unverified_teleparallel_identity(),
        },
        "diagnostic_residuals": {
            "scalar_equation_proxy_L2_unverified_torsion": l2(scalar_proxy, u),
            "Einstein_vacuum_G_t_t_L2_not_matter_subtracted": l2(geo["G_t_t"], u),
            "note": "These are profile diagnostics; full RGS-EMS variational equations are not implemented."
        },
        "physical_quantities": {
            "m_MS_throat": float(mms[center]),
            "m_MS_asymptotic_at_domain_edge": float(mms[-1]),
            "E_EM_finite_domain": Eem,
            "E_EM_status": "FINITE_DOMAIN_ONLY; asymptotic divergence for supplied rho_EM ~ Q^2/(2R^2)",
            "M_ADM_candidate": None,
            "M_ADM_status": "NOT_REPORTED: no validated asymptotic solution/convergence study"
        },
        "profile": {"u": u.tolist(), "R": R.tolist(), "phi": phi.tolist(), "N": N.tolist(),
                    "m_MS": mms.tolist(), "rho_EM": rho.tolist()}
    }


def domain_convergence(params, L_values=(5.0, 10.0, 20.0, 40.0), n=4001):
    """Domain-size audit: demonstrate whether finite-domain EM energy stabilizes."""
    rows = []
    for cfg in CONFIGS:
        for L in L_values:
            u = np.linspace(-L, L, n)
            prof = profile_exact(u, cfg["R0"], cfg["phi0"], cfg["N0"], params["scalar_width"])
            rho = PhysicalQuantities.em_density(prof["R"], params["Q"], params["ell0"])
            Eem = float(np.trapezoid(4*np.pi*rho*prof["R"]**2, u))
            m_edge = float(PhysicalQuantities.misner_sharp_mass(prof["R"], u, params["G"])[-1])
            rows.append({"config_id": cfg["id"], "domain_half_width": L,
                         "grid_points": n, "E_EM_finite_domain": Eem,
                         "m_MS_at_edge": m_edge})
    return rows



def solve_scalar_fixed_background_bvp(cfg, params, L=20.0, n=401, tol=1e-6):
    """
    Solve only the minimally coupled static scalar equation on a FIXED geometry.

    Equation: phi'' + 2 R'/R phi' - V'(phi) = 0, with N constant,
    phi'(0)=0 (reflection symmetry) and phi(L)=0. This is deliberately not
    the full RGS-EMS/tetrad BVP: torsion coupling and metric backreaction are absent.
    The nonconvex potential can admit multiple branches, so domain/grid stability
    is reported and no physical solution status is inferred from solver success.
    """
    b = float(cfg["R0"])
    u = np.linspace(0.0, float(L), int(n))
    def ode(x, y):
        phi, phip = y
        Vp = params["lam"] * phi * (phi**2 - params["phi_v"]**2)
        return np.vstack((phip, Vp - 2.0*x/(x*x+b*b)*phip))
    def bc(ya, yb):
        return np.array([ya[1], yb[0]])
    width = params["scalar_width"]
    amp = float(cfg["phi0"])
    guess_phi = amp*np.exp(-(u/width)**2)
    guess_phip = -2.0*u/(width**2)*guess_phi
    sol = solve_bvp(ode, bc, u, np.vstack((guess_phi, guess_phip)),
                    tol=tol, max_nodes=30000, verbose=0)
    # Independently evaluate the ODE residual on the returned adaptive mesh.
    x = sol.x
    phi, phip = sol.y
    phipp = params["lam"]*phi*(phi**2-params["phi_v"]**2) - 2*x/(x*x+b*b)*phip
    residual = phipp + 2*x/(x*x+b*b)*phip - params["lam"]*phi*(phi**2-params["phi_v"]**2)
    return {
        "config_id": cfg["id"], "status": "SOLVER_CONVERGED_REDUCED_SCALAR_BVP" if sol.success else "SOLVER_FAILED",
        "equation": "phi'' + 2 R'/R phi' - V'(phi)=0 on fixed R=sqrt(u^2+R0^2), N=constant",
        "omissions": ["teleparallel torsion coupling", "metric backreaction", "electromagnetic stress tensor", "full gravitational equations"],
        "domain": {"u_min": 0.0, "u_max": float(L), "initial_grid_points": int(n), "adaptive_mesh_points": int(sol.x.size)},
        "boundary_conditions": {"phi_prime_at_throat": float(phip[0]), "phi_at_outer_boundary": float(phi[-1])},
        "solver": {"success": bool(sol.success), "message": str(sol.message), "tolerance": float(tol),
                   "max_rms_interval_residual": float(np.max(sol.rms_residuals)) if sol.rms_residuals.size else None,
                   "max_abs_reconstructed_ode_residual": float(np.max(np.abs(residual)))},
        "solution_diagnostics": {"phi_throat": float(phi[0]), "phi_max_abs": float(np.max(np.abs(phi))),
                                 "phi_at_L_over_2": float(np.interp(L/2, x, phi))},
        "profile": {"u": x.tolist(), "phi": phi.tolist(), "phi_prime": phip.tolist()}
    }

def scalar_bvp_stability_audit(cfg, params):
    """Vary domain and initial grid; branch changes count as instability, not convergence."""
    rows=[]
    for L, n in [(10.0,201),(10.0,401),(20.0,201),(20.0,401),(40.0,201),(40.0,401)]:
        r=solve_scalar_fixed_background_bvp(cfg, params, L=L, n=n, tol=1e-6)
        rows.append({"L":L,"n_initial":n,"success":r["solver"]["success"],
                     "adaptive_points":r["domain"]["adaptive_mesh_points"],
                     "phi_throat":r["solution_diagnostics"]["phi_throat"],
                     "max_rms_interval_residual":r["solver"]["max_rms_interval_residual"]})
    vals=np.array([r["phi_throat"] for r in rows])
    nontrivial=(np.max(np.abs(vals))>1e-4)
    stable=(np.ptp(vals)<1e-3)
    return {"status":"BRANCH_STABLE_ON_TESTED_DOMAINS" if stable else "BRANCH_NOT_STABLE_UNDER_DOMAIN_GRID_VARIATION",
            "interpretation":"A small ODE residual is not sufficient: compare the selected nonlinear branch as L and the initial mesh vary.",
            "rows":rows,"phi_throat_range":float(np.ptp(vals)),"nontrivial_branch_seen":bool(nontrivial)}

def main():
    params = {"G": 1.0, "Q": 1.0, "ell0": 1.0, "xi": 0.5,
              "V0": 0.1, "lam": 1.0, "phi_v": 0.5, "scalar_width": 1.0}
    u = np.linspace(-10.0, 10.0, 2001)
    results = [evaluate_config(c, u, params) for c in CONFIGS]
    convergence = domain_convergence(params)
    scalar_bvp = [solve_scalar_fixed_background_bvp(c, params, L=20.0, n=401) for c in CONFIGS]
    scalar_stability = [scalar_bvp_stability_audit(c, params) for c in CONFIGS]
    output = {
        "suite": "RGS-PSU diagnostics plus reduced fixed-background scalar BVP",
        "date": "2026-10-28",
        "overall_status": "REDUCED_SCALAR_BVP_SOLVED; FULL_RGS_EMS_BVP_NOT_CLOSED",
        "parameters": params,
        "methodology": {
            "profile": "R(u)=sqrt(u^2+R0^2), phi(u)=phi0 exp(-u^2/width^2), N(u)=N0",
            "solves_field_equations": False,
            "teleparallel_identity": "NOT_VERIFIED until tetrad and inertial spin connection are fixed consistently",
            "ADM": "not reported as a physical result",
            "electromagnetic_sector": "phenomenological density only; NLED action/source equations absent"
        },
        "configurations": results,
        "reduced_scalar_bvp": scalar_bvp,
        "reduced_scalar_bvp_stability_audit": scalar_stability,
        "domain_convergence_audit": {
            "status": "FINITE_DOMAIN_DIAGNOSTIC_NOT_ASYMPTOTIC_CONVERGENCE",
            "rows": convergence,
            "interpretation": ("E_EM is expected to grow approximately linearly with domain size "
                               "for rho_EM ~ Q^2/(2R^2); this is a model/asymptotic issue, "
                               "not a converged finite total energy.")
        },
        "analytic_checks": {
            "R_ansatz_identity": "R'^2 + R*R'' = 1 for R=sqrt(u^2+b^2)",
            "m_MS_ansatz_exact": "m_MS=R0^2/(2R), hence m_MS(throat)=R0/2 and tends to 0 at large |u|",
            "scalar_profile_derivatives": "phi' and phi'' are evaluated analytically",
            "teleparallel_identity": "not asserted: requires a consistently derived tetrad and inertial spin connection"
        },
        "open_requirements": [
            "Derive the full metric/tetrad equations by varying one explicitly fixed action.",
            "Specify tetrad and inertial spin connection and derive T and T^mu from them.",
            "Derive NLED field equation and stress tensor from one Lagrangian; verify Coulomb asymptotics.",
            "Solve the coupled BVP with throat and asymptotic boundary conditions and demonstrate mesh/domain convergence.",
            "Only then extract ADM mass and label any background a solution."
        ]
    }
    out_json = OUT / "RGS_PSU_triple_configs_results.json"
    out_json.write_text(json.dumps(output, indent=2))
    print("="*78)
    print("RGS-PSU CLOSURE — DIAGNOSTICS + REDUCED SCALAR BVP")
    print("="*78)
    for i, r in enumerate(results, 1):
        p, g, q = r["parameters_at_throat"], r["geometry"], r["physical_quantities"]
        print(f"[{i}] {r['label']}: R0={p['R0_target']:.6f}, phi0={p['phi0_target']:.6f}, N0={p['N0']:.6f}")
        print(f"    Teleparallel identity: NOT_VERIFIED (not forced algebraically)")
        print(f"    Scalar proxy L2={r['diagnostic_residuals']['scalar_equation_proxy_L2_unverified_torsion']:.6e}")
        print(f"    m_MS(throat)={q['m_MS_throat']:.6f}; E_EM(domain)={q['E_EM_finite_domain']:.6f}")
        print(f"    ADM: {q['M_ADM_status']}")
    print("-"*78)
    print("REDUCED SCALAR BVP ON FIXED GEOMETRY (NOT FULL RGS-EMS):")
    for r, stab in zip(scalar_bvp, scalar_stability):
        print(f"    {r['config_id']}: success={r['solver']['success']}, phi(0)={r['solution_diagnostics']['phi_throat']:.8e}, "
              f"ODE residual={r['solver']['max_abs_reconstructed_ode_residual']:.3e}, stability={stab['status']}")
    print("-"*78)
    print("DOMAIN-SIZE CHECK (Config 1, E_EM finite-domain):")
    for row in convergence:
        if row["config_id"] == CONFIGS[0]["id"]:
            print(f"    L={row['domain_half_width']:5.1f}: E_EM={row['E_EM_finite_domain']:.6f}, m_MS(edge)={row['m_MS_at_edge']:.6f}")
    print("-"*78)
    print(f"JSON written: {out_json}")
    print("Status: diagnostic only; no converged field-equation solution claimed.")

if __name__ == "__main__":
    main()
