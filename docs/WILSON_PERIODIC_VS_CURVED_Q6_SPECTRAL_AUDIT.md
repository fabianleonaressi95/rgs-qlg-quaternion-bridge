# Wilson spectral audit: periodic reference versus curved Q6

Date: 2026-11-06  
Status: periodic Fourier-symbol corner calculation implemented; finite curved-Q6 Wilson algebra available; curved continuum and resolvent convergence remain OPEN.

## 1. Scope

This audit separates (i) a translation-invariant periodic lattice, where a Brillouin-zone symbol exists, from (ii) the finite reconstructed curved Q6/RGS graph, where vertex-dependent coframes, edge lengths and spin transporters generally destroy translation invariance. A finite spectrum on the second object is meaningful, but it is not a Brillouin-zone spectrum and cannot by itself establish a continuum species count.

## 2. Periodic flat reference

Use the Euclidean Wilson symbol in dimension d:
```
D_E,a(p) = i (1/a) sum_mu gamma_mu sin(p_mu a)
           + (r(a)/a) sum_mu (1-cos(p_mu a)) I.
```
For Hermitian Euclidean Clifford matrices,
```
sigma_a(p)^2 = 1/a^2 sum_mu sin^2(p_mu a)
             + r(a)^2/a^2 [sum_mu (1-cos(p_mu a))]^2.
```
At a Brillouin corner with exactly k components equal to pi/a, the naive kinetic term vanishes and the Wilson mass is m_W(k)=2 k r(a)/a. There are binomial(d,k) such corners. In d=3, the eight naive zeros split into one physical corner k=0 and seven nonphysical corners k=1,2,3. For r(a)>0, all seven have nonzero Wilson mass.

The existing flat-torus prototype records r(a)=r0(a/L)^alpha. For this symbol, fixed physical momentum consistency requires alpha>-1, while divergent corner masses require alpha<1. The joint window is -1<alpha<1. This is specific to the stated periodic symbol and scaling, not a theorem for curved Q6.

New implementation:
- src/wilson_periodic_symbol.py
- tests/test_wilson_periodic_symbol.py

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_wilson_periodic_symbol.py -q
```

The audit enumerates every corner, shell multiplicity, Wilson mass and formula residual. It is an exact finite algebraic check up to floating-point evaluation of trigonometric values.

## 3. Curved finite Q6/RGS graph: existing implementation

The repository already contains a separate graph construction in src/discrete_dirac.py:
- assemble_wilson_laplacian(...) builds the positive edgewise covariant graph Wilson form for unitary spin links and r>=0;
- assemble_spatial_dirac(...) builds the finite Hermitian spatial Dirac candidate;
- assemble_euclidean_wilson_dirac(...) builds D_E=iD_0+W and audits gamma_5-Hermiticity;
- tests check positivity, Hermiticity, local frame covariance, scale covariance, a periodic ring's pi-mode lifting, and the distinction between the Hermitian graph sum and the standard Euclidean Wilson operator.

For a fixed curved graph, these checks establish finite-dimensional algebraic properties only. Positivity follows edgewise from unitary transport and r>=0. Gauge covariance follows from block conjugation of link matrices. Gamma_5-Hermiticity is an algebraic identity for the implemented Euclidean convention. Eigenvalues can be computed, but there is no canonical momentum label or Brillouin corner unless an additional periodic structure is constructed.

The polar lift of coframe changes to SU(2) links is currently a candidate transport prescription. It is not yet proved to approximate parallel transport of a continuum spin connection on a curved refinement family.

## 4. Relation to Mosco convergence

The flat periodic Fourier model in docs/MOSCO_PROTOTYPE_REGULAR_TORUS.md is a reference for a common-Hilbert-space form limit. The scalar P1 theorem in docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md concerns a different form and must not be substituted for the Wilson-Dirac form.

The curved scalar programme in docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md remains OPEN. A curved Wilson/Dirac extension additionally requires:
1. a specified curved refinement family K_h, h->0, with uniform geometric control;
2. compatible spinor bundles and discrete spin transports derived from a declared connection;
3. common Hilbert-space embeddings for spinors and a precise operator/form domain;
4. uniform estimates for the Wilson kernel and control of low-energy modes;
5. liminf and recovery estimates for the relevant closed forms, or an operator theorem sufficient for strong resolvent convergence;
6. separate analysis of chirality/index and possible boundary or local-system terms.

Scalar Dirichlet-form Mosco convergence does not automatically imply convergence of the spinorial Wilson operator. A finite spectral gap on one Q6 graph does not prove compactness or resolvent convergence as h->0.

## 5. Umeda's suggested point alpha=(pi/2,-pi/2)

The momentum point quoted in the correspondence must be kept distinct from the scaling exponent alpha used in the periodic-torus note. They are different symbols with different meanings. Before using Umeda's point, the precise operator convention and momentum-coordinate map from the cited paper must be checked. The present Q6 graph code does not identify that point with a Q6 eigenmode or curved-graph spectral sector. No claim about its spectral estimates is made here.

## 6. Evidence ledger

| Statement | Status |
|---|---|
| Periodic Wilson symbol and D_E^dagger D_E formula | ANALYTICALLY DERIVED under stated Clifford convention |
| Brillouin-corner masses 2kr/a and multiplicities binomial(d,k) | ANALYTICALLY DERIVED |
| Corner enumeration for selected d,a,r | NUMERICALLY VERIFIED by the new audit |
| Finite graph Wilson positivity for unitary links, r>=0 | PROVED edgewise; tested |
| Finite graph Wilson gauge covariance | PROVED algebraically; tested |
| Finite graph Euclidean gamma_5-Hermiticity | PROVED algebraically for D_E=iD_0+W; tested |
| Eight continuum species or species count on curved Q6 | OPEN / NOT CLAIMED |
| Curved-Q6 Brillouin symbol | NOT DEFINED |
| Curved spin-transport consistency under refinement | OPEN |
| Curved Wilson Mosco/strong-resolvent convergence | OPEN |
| Wilson-kernel gap uniform in refinement | OPEN |
| Ginsparg-Wilson relation and index convergence | OPEN |

## 7. Next proof gate

Do not infer a curved doubler count from the 64-vertex graph. Define a periodic reference family and a separate curved refinement family, then compare low-energy operators through explicit identification maps.

At each curved refinement level, record the smallest and largest eigenvalues of the positive Wilson kernel, multiplicity below a declared physical cutoff, gamma_5-Hermiticity residual, link unitarity and gauge-covariance residuals, dependence on mesh scale and geometric distortion, and stability of low-energy spectral projectors. These finite-level tests are diagnostic only; uniform estimates and a common-space convergence proof are still required.
