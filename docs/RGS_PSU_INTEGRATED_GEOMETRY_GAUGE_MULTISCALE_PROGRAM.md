# Integrated RGS–PSU–Q6–PRISM geometry, gauge and convergence program

**Status:** implementation and finite audits added; continuum and physical theorems remain open.  
**Design rule:** all finite complexes, gauge groups, local-frame transitions and continuum limits remain separately typed and explicitly mapped.

## 1. The geometric architecture

The integrated program has three layers:

1. **Finite combinatorial/geometric data:** the 12-vertex Möbius subcomplex, the 64-vertex Q6 2-skeleton, the 69,120-vertex PRISM product graph, the 64B coordinate seed and the 188D constrained tangent sector.
2. **Local geometry and fields:** local coframes, frame transitions, metric-compatible connection, torsion, curvature, Clifford action, Dirac operator and gauge holonomy.
3. **Refinement and continuum analysis:** multifractal partition scaling, empirical equidistribution, measured Gromov–Hausdorff control, Mosco convergence, Yang–Mills and hydrodynamic limits.

The finite pieces are connected only by explicit maps, incidence data, coframes or representations. Equal cardinalities or similar filenames are not identifications.

## 2. Möbius subcomplex versus Q6

The selected Möbius subcomplex is a genuine six-square subcomplex of the Q6 2-skeleton:

| Complex | Vertices | Edges | Square faces |
|---|---:|---:|---:|
| Möbius subcomplex | 12 | 18 | 6 |
| Full Q6 2-skeleton | 64 | 192 | 240 |

The finite audit checks:
- B1 B2 = 0;
- Euler characteristic chi = 12 - 18 + 6 = 0;
- Betti numbers (b0,b1,b2) = (1,1,0);
- one boundary component;
- an orientation-propagation conflict and nontrivial dual-cycle orientation holonomy.

The Möbius six faces are selected from, and remain distinct from, the other 234 Q6 faces. No 12/18/6 count is substituted for the full 64/192/240 Q6 complex.

### Local frames and reflections

The Möbius dual graph carries explicit Z2 transition signs. Local-frame changes act by

[
t_{ij}'=h_jt_{ij}h_i^{-1},qquad h_i,h_j,t_{ij}in{pm1}.
]

Scalar, sign and regular representations are implemented and their composition laws tested. The product around the selected dual cycle is -1.

This is an orientation-line transition system. It does **not** by itself define a physical reflection field, torsion, a Pin structure or a fermion bundle. A Pin+ or Pin- lift requires the relevant Stiefel–Whitney obstruction test (respectively w2 or w2+w1^2, with convention stated) and explicit reflection lifts satisfying the chosen Clifford relations. The current 2I spinor action is spatial Spin(3), not a reflection lift.

## 3. Seeded SU(2) gauge-covariance audit

The Möbius lattice has independent internal SU(2) link variables. This internal gauge bundle is kept separate from the orientation-line transitions.

For a link U_gh, the convention is

[
U_{gh}mapsto G_hU_{gh}G_g^dagger.
]

For a path P:v0 to vn,

[
Pmapsto G_{v_n}PG_{v_0}^dagger,
]

and for a closed loop based at v0,

[
Hmapsto G_{v_0}HG_{v_0}^dagger,qquad
operatorname{Tr}Hmapstooperatorname{Tr}H.
]

The integrated runner uses fixed seeds 7, 17, 31, 20261103, records a SHA-256 digest of each generated link/gauge input, and records maximum open-path, closed-loop conjugation, Wilson-trace and reverse-link residuals. Tests cover all six face loops and their open prefixes. These are finite randomized covariance checks, not a proof of gauge-theory continuum convergence.

## 4. Q6 incidence bundle and Q6–PRISM map

src/q6_complex.py defines the source-controlled canonical Q6 order:
- vertices 0..63, six-bit big-endian;
- edges lexicographically sorted with u<v;
- square faces ordered by axis pair a<b, then base vertex;
- B1:C1 to C0, B2:C2 to C1, with B1 B2 = 0.

The exporter writes results/artifacts/Q6_canonical_incidence_v1.npz and a JSON manifest with shapes and SHA-256. This explicit schema is reproducible, but is **not claimed byte-compatible** with the separately recovered historical NPZ until every array and ordering has been compared. Existing Q6 .npy data are preserved.

The first explicit Q6–PRISM map is the inverse-Gray-code injection

[
vin Q_6mapsto(operatorname{gray}^{-1}(v),0,0,0,0)
in C_{64}square C_6^3square C_5.
]

It maps the Gray Hamiltonian cycle of Q6 onto the C64 factor. It is injective and preserves 64 edges, exactly one third of Q6's 192 edges. It is **not** a full Q6 graph embedding or graph isomorphism. The other Q6 edges require a richer fibre map if they are to be preserved.

## 5. Binary icosahedral group, representations and intertwiners

The explicit group is the 120 unit-quaternion 600-cell vertex set 2I subset SU(2), scalar-first Hamilton convention. The audit checks:
- 120 distinct unit elements;
- all 14,400 products close in the set;
- inverse closure;
- 60 distinct induced SO(3) rotations after identifying q ~ -q;
- orthogonality and determinant of induced rotations.

The selected representations (trivial, fundamental 2D, vector 3D) are not a complete character table. Their homomorphism/unitarity audits and all pairwise Hom-space dimensions are computed. Schur-lemma expectations apply to these selected irreducibles only. The Dirac-spinor restriction to spatial Spin(3) has two equivalent fundamental doublets; its complex commutant is Mat_2(C), dimension four.

## 6. Clifford conventions, grades and steerable kernels

The Lorentzian convention is

[
eta=operatorname{diag}(-1,+1,+1,+1),qquad
{gamma^mu,gamma^
u}=2eta^{mu
u}I_4.
]

The implementation builds all 16 ordered blades, projects matrix elements into grades 0 through 4, and checks grade preservation under the adjoint action B -> S(q) B S(q)^-1 for all 120 spatial 2I lifts. Intertwiner bases and residuals are tested, including the four-dimensional Dirac commutant.

This does not implement Lorentz boosts or Pin reflections. A separate reflection action must be constructed and tested before coupling to the Möbius orientation line.

## 7. Differential operators, coframe connection and the Cartan sequence

For the finite chain complex, define

[
d_0=B_1^T:C^0	o C^1,qquad d_1=B_2^T:C^1	o C^2.
]

The audit checks d1 d0 = 0 and builds the unit-weight combinatorial Hodge Laplacians

[
Delta_0=d_0^Td_0,quad
Delta_1=d_0d_0^T+d_1^Td_1,quad
Delta_2=d_1d_1^T.
]

Weighted Hodge stars, discrete wedge/cup products and a coframe-valued exterior derivative are not yet defined by this module.

The next geometric construction must obey

[
elongrightarrowomegalongrightarrow T^a
longrightarrow R^{ab}longrightarrow D_{m Clifford}
longrightarrowmathrm{TEGR/Dirac}.
]

In a local frame, the continuum transformation laws are

[
emapstoLambda e,quad
omegamapstoLambdaomegaLambda^{-1}
+Lambda dLambda^{-1},
]

[
T^a=de^a+omega^a{}_bwedge e^b,qquad
R^a{}_b=domega^a{}_b+omega^a{}_cwedgeomega^c{}_b.
]

The discrete implementation still needs a metric-compatible edge connection, a compatible discrete wedge/cup product, and explicit Cartan torsion/curvature residuals. A polar Spin(4) alignment alone is not this construction. Teleparallel gravity additionally requires a flat inertial spin connection and nonzero torsion; TEGR equivalence to Einstein gravity depends on the convention-consistent boundary identity. No torsion or TEGR claim is inferred from Möbius non-orientability.

## 8. 188D Hessian and T2u provenance

The integrated inventory scans every data/rgs_188d/*.npy, records file SHA-256, canonical array hash, shape and dtype, compares all same-shape 188x188 operators, and checks candidate inverse pairs.

A 192x192 array is only a Hessian candidate: its formula, action, constraints and source provenance must be recovered before certifying K_tangent = Q^T H Q. The current repository evidence validates the A14 constraint Jacobian/tangent basis but explicitly does not validate the original 192x192 energy Hessian or the claimed 29/3/156 spectrum.

The T2u bundle audit records hashes of all_T2u_arrays.npz, K_tangent_rebuilt.npy, K_T2u_rebuilt.npy, B_T2U_TAN_rebuilt.npy and J0_rebuilt.npy; when the archived eigenvector and first-triplet index keys are present, it recomputes the 188x3 triplet and reports the relative residual. The historical report's sha256-not-computed placeholders are replaced by runtime hashes in the generated audit JSON. Synthetic J_A,J_B remain quarantined as unrecoverable.

## 9. Multifractality, equidistribution, GH and Mosco

src/rgs_multiscale.py provides finite diagnostics:

- Partition sums Z(q,epsilon) = sum_i mu_i(epsilon)^q, fitted tau(q), and generalized dimensions Dq = tau(q)/(q-1); the q=1 information dimension uses the entropy limit.
- Haar moment comparisons for finite point sets on S^(d-1), including S3 and the 2I 120-point set.
- Distortion of finite metric correspondences, with half-distortion bound only when the correspondence projects onto both finite spaces.
- Finite liminf/recovery sample diagnostics for quadratic forms.

A finite 2I orbit is atomic and cannot equal continuous Haar measure. Its moment symmetries are useful checks, not Weyl equidistribution. Equidistribution requires a sequence of increasingly refined point sets and weak convergence against a determining class of functions.

Likewise, a finite metric distortion is not Gromov–Hausdorff convergence of the RGS family. That needs common-scale metrics, refinement maps/correspondences, covering estimates and compactness. Mosco convergence separately requires a common Hilbert-space framework, a liminf inequality for every weakly convergent sequence, and recovery sequences for every vector in the limit form domain. The flat-torus Fourier prototype is not transferred to curved Q6/RGS without proof.

The multifractal link is operational through scale-dependent measures and spectra, not by declaring the 64 esagrams or a finite 2I orbit automatically equidistributed. The intended bridge is: Q6/PRISM refinement -> nested geometric partitions -> empirical measures and partition sums -> measured metric convergence -> closed quadratic forms and Mosco convergence -> limiting differential operators.

## 10. Yang–Mills, Navier–Stokes and teleparallel dynamics

The finite Q6 Wilson-plaquette action is a lattice gauge prototype. To obtain continuum Yang–Mills one must specify a sequence of gauge bundles/links, scale the plaquette action, prove gauge-covariant consistency and control the measure/action limit.

Navier–Stokes is a distinct hydrodynamic limit. It requires a conserved momentum variable, a scaling regime, constitutive closure, viscosity from a controlled Green–Kubo or kinetic limit, incompressibility/pressure constraints and convergence of the discrete dynamics. A common graph, a shear parameter or a shared spectral gap does not establish Yang–Mills/Navier–Stokes equivalence.

The proposed unification is therefore a shared geometric carrier with separately derived sectors:
- coframe/connection sector for Riemann–Cartan and teleparallel geometry;
- link/plaquette sector for Yang–Mills;
- matter/spinor sector for Clifford–Dirac;
- hydrodynamic sector only after a conserved stress/momentum limit and Green–Kubo derivation;
- multifractal/measure sector governing the refinement and spectral scaling.

The sectors may couple through one action only after the variables, gauge symmetries, dimensions, units and variation equations are written explicitly.

## 11. Reproducible run

    PYTHONPATH=. python scripts/rgs_psu_integrated_audit.py
    python -m pytest tests/test_integrated_mobius_psu_audit.py tests/test_mobius_complex.py tests/test_mobius_gauge_audit.py tests/test_mobius_representations.py tests/test_binary_icosahedral.py tests/test_representations_2i.py tests/test_clifford_qlg.py tests/test_q6_complex.py -v

The JSON report includes seeds, tolerances, source and data hashes, residuals, operator inventories and explicit continuum non-claims. CI uploads the JSON and canonical Q6 NPZ/manifest as a workflow artifact.

The integrated JSON writer converts complex representation entries to explicit `{real, imag}` pairs; no phase or imaginary component is discarded.

### 188D operator identity classification

The operator inventory now records SVD numerical rank/nullity and symmetry residual for each 188x188 array. For the candidate pairs (K0,K0_inv), (K0 (1),Kinv), and (K_tangent_rebuilt,Kinv_rebuilt), it reports both identity products and the four Moore–Penrose conditions:
[
ABA=A,quad BAB=B,quad (AB)^T=AB,quad (BA)^T=BA.
]
This distinguishes an ordinary inverse from a pseudoinverse (where AB and BA are orthogonal projectors) and from an unrelated or incorrectly paired matrix. A product residual of sqrt(3) against the identity can be consistent with a rank-185 projector if the three-dimensional kernel is independently confirmed; it is not by itself proof of that interpretation. The integrated JSON records the rank, projector idempotence, trace and all four residuals before assigning a classification.
