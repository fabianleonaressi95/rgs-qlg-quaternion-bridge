"""Reproducible Q6 hypercube / Green--Kubo / supergolden audit."""
import json
import numpy as np
from src.q6_hypercube_gk_supergolden import (
    green_kubo_ou,q6_closed_two_cochain_projector,q6_hypercube_cosine_law,
    supergolden_constant,supergolden_geometric_sum,supergolden_refinement_scales,
)
def main():
    curvature=q6_hypercube_cosine_law(np.linspace(-.015,.015,240),w=2.)
    p=q6_closed_two_cochain_projector()
    # Toy OU values only: not a calibrated physical Q6 transport model.
    gk=green_kubo_ou(np.diag([.5,1.,2.]),np.diag([1.,.5,.25]),
                     np.array([1.,.5,.25]),temperature=1.,volume=1.)
    psi=supergolden_constant()
    out={"hypercube_law":curvature,
         "closed_projector":{"rank":int(np.linalg.matrix_rank(p)),
             "idempotence_residual":float(np.linalg.norm(p@p-p))},
         "green_kubo_toy_ou":{"eta":gk.eta,"correlation_integral":gk.correlation_integral,
             "rates":gk.eigenvalues.tolist(),"modal_amplitudes":gk.modal_amplitudes.tolist(),
             "convergent":gk.convergent,"scope":"toy OU example, not physical viscosity"},
         "supergolden":{"psi":psi,"polynomial_residual":psi**3-psi**2-1.,
             "infinite_geometric_sum_a0_1":supergolden_geometric_sum(1.),
             "first_six_scales":supergolden_refinement_scales(1.,6).tolist()}}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__": main()
