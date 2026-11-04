# Cut-and-project temporal-axis status: candidate versus no-go audit

The golden cut-and-project decomposition is Euclidean and kinematic:
\[
\mathbb R^6=E_\parallel\oplus E_\perp.
\]
It does not by itself supply a clock, a time derivative, a causal order, or a Lorentzian signature.

The recovered QAT temporal-closure audit reports:
- the Q6 Wilson plaquette Hessian is a spatial operator;
- the fourth A14 constraint normal is the fixed-radius/fixed-scale radial normal, not time;
- the quaternionic/SU(2) fibre supplies internal complex structures and transport, but does not select a preferred evolution parameter or Lorentzian sign;
- a wave equation such as \(\partial_s^2 A+c_*^2H_{U(1)}A=0\) remains conditional while \(s\) and \(c_*\) are external.

The torsion-coupled SVD routine in `src/rgs_psu_cut_project.py` is a *candidate selection experiment*, not a counterexample to that no-go result. It asks whether an independently measured 64-site torsion trace statistically selects a stable direction in \(E_\perp\). It does not establish that this direction is physical time. It must be tested against row-permutation nulls and bootstrap angular stability.

The routine therefore refuses to invent its missing input. In particular, the archived 192x4 edge coframe is not a 64xm torsion trace and is not silently substituted. The prior exploratory output directory `RGS_PRISM_temporal_axis_audit` is not present as a committed, hash-verified source in the current repository or retrieved Library files.

To close the temporal branch, first recover and hash the independent torsion trace or derive it from the declared RGS connection/tetrad. Then compare the SVD axis with a causal transfer operator derived from quaternionic transport. A Lorentzian principal symbol with signature \((-+++ )\), positivity/unitarity, and a continuum 3D isotropic limit are separate mandatory tests.

**Status:** spatial cut-and-project = implemented finite construction; torsion-selected internal axis = conditional candidate; emergent time and Lorentzian spacetime = open.
