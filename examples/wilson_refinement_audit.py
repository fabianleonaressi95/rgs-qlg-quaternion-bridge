"""Wilson refinement, finite overlap/GW, and Q6 Dirac-Kahler audit.

Run from repository root:
    python examples/wilson_refinement_audit.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from src.dirac_kahler import assemble_q6_dirac_kahler
from src.overlap_dirac import assemble_overlap_from_wilson\nfrom src.mosco_torus import torus_mosco_sampling_audit
from src.wilson_refinement import (
    euclidean_wilson_singular_value, hermitian_scalar_sum_branches,
    refinement_audit, scalar_sum_accidental_zero_angle, wilson_brillouin_corners,
)


def main():
    refinements_fixed_r=refinement_audit(
        sizes=(16,32,64,128,256),dimension=3,physical_length=1.0,
        mode=(1,0,0),r0=1.0,r_exponent=0.0
    )
    refinements_running_r=refinement_audit(
        sizes=(16,32,64,128,256),dimension=3,physical_length=1.0,
        mode=(1,0,0),r0=1.0,r_exponent=0.5
    )
    kahler,kahler_audit=assemble_q6_dirac_kahler()\n    mosco_torus=torus_mosco_sampling_audit(\n        sizes=(16,32,64,128,256),dimension=3,physical_length=1.0,\n        mode=(1,0,0),r0=1.0,r_exponent=0.0\n    )

    # A transparent finite-dimensional overlap/GW algebra check. This is a
    # matrix identity test, not an RGS-Q6 overlap spectrum.
    gamma5=np.diag([1.0,-1.0]).astype(complex)
    naive=np.array([[0.0,1.0],[1.0,0.0]],dtype=complex)
    wilson=0.3*np.eye(2,dtype=complex)
    overlap,overlap_audit=assemble_overlap_from_wilson(
        naive,wilson,gamma5,lattice_spacing=1.0,rho=1.0
    )

    corners=wilson_brillouin_corners(3,1.0/256.0,1.0)
    theta_star=scalar_sum_accidental_zero_angle(1.0)
    scalar_branches=hermitian_scalar_sum_branches([theta_star],1.0,1.0)
    euclidean_sigma_star=euclidean_wilson_singular_value([theta_star],1.0,1.0)
    report={
        "audit":"RGS_QLG_PHASE9_REFINEMENT_AND_FERMION_FORMULATIONS",
        "wilson_refinement_fixed_r":refinements_fixed_r,
        "wilson_refinement_running_r_alpha_0_5":refinements_running_r,
        "three_dimensional_brillouin_corners_at_a_1_over_256":{
            "total_corners":len(corners),
            "physical_corner_count":sum(c["pi_directions"]==0 for c in corners),
            "doubler_corner_count":sum(c["pi_directions"]>0 for c in corners),
            "masses_by_pi_directions":{
                str(k):2.0*k*256.0 for k in (1,2,3)
            }
        },
        "one_dimensional_convention_audit":{
            "r":1.0,"accidental_zero_angle_of_hermitian_D0_plus_W":theta_star,
            "hermitian_scalar_sum_branches_at_angle":list(scalar_branches),
            "euclidean_wilson_singular_value_at_same_angle":euclidean_sigma_star,
            "interpretation":"D0+W has an accidental zero; i*D0+W does not"
        },
        "flat_torus_mosco_prototype":mosco_torus,\n        "q6_dirac_kahler":kahler_audit,
        "finite_overlap_ginsparg_wilson":overlap_audit,
        "epistemic_status":{
            "fixed_physical_momentum_consistency":"NUMERICALLY_AUDITED_FROM_EXACT_FOURIER_SYMBOL",
            "doubler_mass_scaling":"ANALYTICALLY_DERIVED_AND_NUMERICALLY_SAMPLED",
            "q6_dirac_kahler":"FINITE_UNWEIGHTED_COMBINATORIAL_OPERATOR",
            "overlap":"FINITE_MATRIX_GW_IDENTITY_TEST_ON_TOY_KERNEL",
            "rgs_q6_overlap_operator":"NUMERICALLY_AUDITED_FROM_ARCHIVE; SEE rgs_q6_overlap_audit.json",
            "mosco_convergence":"NOT_PROVED",
            "continuum_locality_of_overlap":"NOT_PROVED"
        }
    }
    out=Path("results/json/wilson_refinement_audit.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    if not (refinements_fixed_r["simultaneous_scaling_window"]
            and refinements_running_r["simultaneous_scaling_window"]
            and kahler_audit["d_squared_residual"]<1e-12
            and kahler_audit["square_to_hodge_laplacian_residual"]<1e-10
            and overlap_audit["ginsparg_wilson_residual"]<1e-10
            and abs(scalar_branches[0])<1e-10
            and euclidean_sigma_star>1e-8\n            and mosco_torus["fixed_mode_form_converges"]\n            and mosco_torus["doubler_mass_diverges"]):
        raise SystemExit("Phase 9 audit failed")


if __name__=="__main__":
    main()
