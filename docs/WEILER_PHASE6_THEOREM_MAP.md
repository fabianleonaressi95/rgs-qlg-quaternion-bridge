# Phase 6 — Weiler-compatible geometry and discrete Dirac operator

## Mathematical contract

The listed Weiler results are translated into finite-graph obligations. Their Euclidean/manifold hypotheses are not silently transferred to Q6.

| Weiler result | Constraint | Phase 6 status |
|---|---|---|
| 3.2.1, 4.3.1 | Equivariant linear maps are relative-position kernel transforms; graph matrices are not automatically convolutions. | Finite-support field transform implemented; metric graph automorphism classification open |
| 3.2.2, 4.3.2 | Bias belongs to the invariant subspace of spatial and fibre actions. | No bias added; invariant-bias projector open |
| 3.2.3, 4.3.3 | Pointwise nonlinearities commute with fibre representations and are position-independent. | No nonlinear layer added |
| 3.2.4–3.2.8, 4.3.4–4.3.9 | Pooling/subsampling require invariant windows, appropriate permutation/unitary reps and correct measure factors. | Not used in Dirac construction |
| 4.5.1 | Aff(G) group convolution equals regular-field steerable convolution under its hypotheses. | Not invoked: Q6 is not assumed to be a continuous affine group |
| 5.3.1 | Steerable kernels are built from representation intertwiners / Wigner–Eckart data. | Selected 2I Hom spaces already tested; complete Q6 kernel basis open |
| 12.2.6 | Compact support and geometric assumptions ensure well-defined kernel-field transforms. | Finite sums are algebraically well-defined; this is not a continuum integral theorem |
| 12.2.8–12.2.9 | Coordinate formula includes local gauges and transporters. | Compact-support transform explicitly includes source-to-target transport |
| 13.1.3–13.1.4 | Isometries act compatibly on local gauges and transporter pullbacks. | Need actual metric-preserving automorphism subgroup of reconstructed X,E |
| 13.2.4–13.2.7 | Equivariance iff the kernel field is invariant under the relevant isometries. | Exact audit open until that subgroup and kernel field are specified |
| 13.3.1–13.3.3 | Quotient representatives and homogeneous-space convolution require quotient/homogeneity hypotheses. | Not invoked before transitivity is proved |
| 15.1.2–15.1.6 | Affine chart transitions induce compatible G-gauge atlases. | Local coframe convention fixed; global atlas/cocycle open |
| 15.2.1–15.2.2 | Coordinate-independent Euclidean GM convolution reduces to coordinate form and remains affine equivariant. | Not claimed for this finite non-Euclidean graph |
| 17.2.1–17.2.2 | Spherical tangent kernels correspond to steerable kernels under geodesic hypotheses. | No spherical geodesic model supplied for the 64 vertices |
| 17.3.1 | Gnomonic projection is a warped exponential map under spherical hypotheses. | Not used |
| B.3.7 | Orbit-stabilizer controls finite group orbits/stabilizers. | Required to classify actual metric-preserving subgroup |
| B.5.10 | Schur's lemma constrains intertwiners of irreducibles. | Selected 2I Hom spaces tested; not a complete catalogue |
| B.5.14, B.5.16 | Compact unitary representations are completely reducible. | Applies to finite unitary 2I reps, not arbitrary GL(3) matrices |
| B.5.22 | Peter–Weyl decomposes functions on compact groups. | Finite-group analogue exists; full Q6 kernel harmonic analysis open |
| J.0.1 | Regular feature fields correspond to scalar functions on the G-structure. | Conceptual target; total-space implementation open |

## Vierbein and curved Clifford field

Input E has shape (64,3,3), as defined in the RGS reconstruction. The static ansatz is
e^A_mu(i)=diag(1,E_i), with eta=diag(-1,+1,+1,+1).
The inverse coframe defines Gamma^mu(i)=e_A^mu(i) gamma^A. The audit
checks all 16 anticommutators at every vertex. The time leg is assumed static;
no lapse, shift, or time-dependent vierbein is derived.

## Discrete spatial Dirac candidate

For canonical edge i<j, define M_ji=E_j E_i^{-1}, take its right-polar
rotational factor R_ji in SO(3), and lift it to S_ji in SU(2). This is an
explicit discrete frame-transport ansatz, not yet the Levi-Civita spin connection.

For d_ij=x_j-x_i, let n_i=E_i d_ij/||E_i d_ij|| and similarly at j.
With c_i(n)=n_a gamma_i^a, set

H_ji = i/(2||d_ij||) [ c_j(n_j) S_ji + S_ji c_i(n_i) ],
H_ij = H_ji^dagger.

For vertex volumes v_i=det(E_i)>0, define D=V^(-1/2) H V^(-1/2).
Thus D=D^dagger by construction. This proves finite-dimensional
self-adjointness for the stated weighted inner product, not continuum
consistency, absence of fermion doubling, or Mosco convergence.

## Evidence labels and open obligations

- PROVED algebraically: the curved Clifford relation follows from the inverse tetrad when E_i is invertible and orientation-preserving.
- NUMERICALLY VERIFIED by tests: pointwise Clifford identities, SU(2) unitarity of polar lifts, Hermiticity and coordinate-dilation scaling of the finite operator.
- CANDIDATE: polar-factor spin transport and nearest-neighbour spatial Dirac operator.
- OPEN: convergence of links to Levi-Civita spin connection; smooth-field consistency; spectral convergence; fermion doubling; 4D time-dependent Lorentzian operator; Mosco convergence.
- OPEN: identify the subgroup of Q6 automorphisms preserving the actual reconstructed X_i,E_i. Only this subgroup can support an exact Weiler-style isometry-equivariance audit.
- OPEN: complete 2I/Clebsch–Gordan kernel basis and invariant bias/pooling spaces.

No claim of physical unification or of a proven continuum Dirac equation is made.


## Phase 8 — Gauge-covariant Wilson term and fermion-doubling audit

The spatial Dirac candidate now has an optional Wilson graph Laplacian. For each
canonical edge i<j, with spin transporter S_ji mapping the source fibre i to
the target fibre j, and positive edge length l_ij, define

`w_ij = r/(2 l_ij)`, with r >= 0, and the unnormalized blocks

`L_ii += w_ij I_4`, `L_jj += w_ij I_4`,
`L_ji -= w_ij S_ji`, `L_ij -= w_ij S_ji^dagger`.

The volume-weighted positive graph Laplacian is `W = V^(-1/2) L V^(-1/2)`.
The legacy option `assemble_spatial_dirac(..., wilson_r=r)` returns the
Hermitian scalar perturbation `D_H = D_naive + W`; it is **not** the standard
Euclidean Wilson-Dirac matrix in this repository's Hermitian convention.

The correct Euclidean Wilson operator is

`D_E = i D_naive + W`,

implemented by `assemble_euclidean_wilson_dirac(...)`. It satisfies
`D_E^dagger = Gamma5 D_E Gamma5` when the naive operator anticommutes with
Gamma5 and W commutes with Gamma5. The overlap construction uses this
gamma5-Hermitian Euclidean matrix. The default `wilson_r=0` preserves the
Phase 6 naive operator.

### Algebraic properties

- **Hermiticity:** opposite edge blocks are adjoints and V is positive diagonal,
  hence W is Hermitian.
- **Gauge covariance:** for local fibre transformations G_i, the link law is
  `S_ji -> G_j S_ji G_i^dagger`; therefore
  `W -> G W G^dagger`, where G is block diagonal over vertices.
- **Positivity:** before volume normalization, the quadratic form is a sum of
  `w_ij ||psi_j - S_ji psi_i||^2`; hence W is positive semidefinite.
- **Chiral breaking:** W is scalar in spinor space and commutes with
  Gamma5. Thus `D_E=iD_naive+W` is gamma5-Hermitian but does not have the
  naive exact anticommutation. The legacy Hermitian sum `D_H=D_naive+W`
  also breaks the anticommutator, but this fact alone does not make it the
  standard Euclidean Wilson operator.
- **Scaling:** if all edge lengths are dilated by alpha while E and r are fixed,
  W scales as alpha^(-1), matching the scaling of the implemented first-order
  graph operator. A continuum Wilson coefficient normally carries an explicit
  lattice spacing; its refinement scaling must be specified separately.

### Periodic one-dimensional check

On a uniform periodic ring with spacing a and trivial links,

`W(p) = (r/a) [1 - cos(pa)]`, and the naive Hermitian symbol magnitude is
`E0(p)=(2/a)|sin(pa)|`.

The Euclidean Wilson matrix `D_E=iD0+W` has singular-value magnitude
`sqrt(E0(p)^2+W(p)^2)`: it vanishes at p=0 and equals 2r/a at p=pi/a.
This is the standard Wilson lifting mechanism in the current convention.

**Important correction:** the legacy Hermitian sum `D_H=D0+W` has branches
`W(p) +/- E0(p)`. In one dimension it has an additional accidental zero
at `pa=2 arctan(2/r)`, besides p=0. Therefore the earlier test of W(p) at
p=pi/a alone did not prove that the combined Hermitian sum is free of
additional zeros. No no-doubling claim is based on that sum.

### Evidence and limits

Tests cover Hermiticity, positive semidefiniteness, local SU(2) gauge covariance,
coordinate-dilation scaling, exact p=0 and p=pi/a Wilson eigenvalues on a
16-site ring, and nonzero chiral anticommutator after adding W.

This establishes the stated finite-graph algebra and the one-dimensional
Fourier mechanism. It does **not** by itself prove absence of doublers on the
curved Q6 graph or in a continuum limit, nor Mosco/spectral convergence.
The Wilson parameter's physical scaling under refinement and the relation
between this Hermitian finite operator and a covariant continuum Dirac
operator remain open obligations.


### RGS/Q6 archive audit (64 vertices, 192 edges)

A separate numerical evaluation using the supplied `RGS_PSU_all_matrices(1).npz`
archive and the Phase 6 vierbein reconstruction constants gives, for r=1:

- Wilson matrix dimension 256x256; Hermiticity residual
  `1.39e-15`; minimum eigenvalue `5.221e-4`; maximum eigenvalue
  `5.13545`.
- Maximum spin-link unitarity residual `7.02e-16`; local gauge-covariance
  residual under deterministic random sitewise SU(2) transformations
  `6.86e-15`.
- Naive Dirac matrix Hermiticity residual `3.16e-15` and chiral
  anticommutator residual `0`.
- Wilson-corrected matrix Hermiticity residual `3.46e-15`; chiral
  anticommutator norm `89.2471`.
- The smallest absolute eigenvalues of the corrected Q6 matrix begin at
  `0.0127313` (eightfold cluster), then `0.0293490`. The uncorrected
  operator's smallest absolute eigenvalue is about `1.25821`.

These are **numerical observations on this one finite reconstructed graph**.
They do not label any Q6 mode as a physical fermion or a doubler: Q6 has no
Brillouin torus or canonical momentum, and the chosen polar links have
nontrivial graph holonomy. In particular, the small eigenvalue cluster after
adding W is not evidence of a removed doubler. The archive-derived values are
recorded in `results/json/rgs_q6_wilson_audit.json`; the 1D periodic ring
mechanism is independently tested in CI.
