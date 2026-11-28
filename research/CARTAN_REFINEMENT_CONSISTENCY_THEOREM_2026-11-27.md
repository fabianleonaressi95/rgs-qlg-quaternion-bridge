# Cartan refinement: local consistency theorem and continuum-limit gates (2026-11-27)

## Status

**Status: CONDITIONAL local consistency proposition for smooth sampled fields; the Taylor/Magnus proof is outlined below, while fully expanded coordinate remainder estimates remain to be written before promotion to PROVED.**
**Convergence of arbitrary discrete fields or of discrete field equations: open.**
**Physical identification of the Q6/PSU complex with a 4D Lorentzian refinement family: open.**

This note advances the affine-Poincare construction without promoting the finite tetrahedral boundary word to a continuum theorem. It specifies the normalization and the exact scope of the result.

## 1. Geometric setting

Let (M) be an oriented four-dimensional manifold, let (einOmega^1(M,mathbb R^{1,3})) be a smooth nondegenerate coframe, and let (omegainOmega^1(M,mathfrak{so}(1,3))) be an independent smooth Lorentz connection. Use
[
T=de+omegawedge e,qquad
R=domega+omegawedgeomega.
]
No torsion-free condition is imposed. No teleparallel condition is imposed. The Riemann-Cartan, teleparallel ((R=0)), and torsion-free ((T=0)) sectors remain distinct.

Let ({mathcal K_h}) be a family of oriented simplicial triangulations of a relatively compact region of (M), with maximal simplex diameter (h	o0). Assume uniform shape regularity and a uniform bound on aspect ratios. The simplices are embedded in coordinate charts; all estimates below are local and uniform on compact subsets.

For an oriented edge (i	o j), let (L_{ij}in SO^+(1,3)) be the connection parallel transport with the convention
[
L'_{ij}=g_iL_{ij}g_j^{-1}.
]
Let (E_{ij}inmathbb R^{1,3}) be the coframe integral along that edge, parallel-transported to the initial vertex frame. Choose the orientation convention so that
[
L_{ji}=L_{ij}^{-1},qquad E_{ji}=-L_{ij}^{-1}E_{ij}.
]
Then (E'_{ij}=g_iE_{ij}). Package an edge as (G_{ij}=(L_{ij},E_{ij})in SO^+(1,3)ltimesmathbb R^{1,3}), with
[
(H,a)(K,b)=(HK,a+Hb).
]

## 2. Face observables and normalization

For an oriented triangular face (f=[i,j,k]), define
[
F_f=G_{ij}G_{jk}G_{ki}=(H_f,	au_f).
]
Here (H_f) is the Lorentz holonomy and (	au_f) is the transported coframe-closure vector. Both are based at the initial vertex and transform covariantly:
[
H_f' =g_iH_fg_i^{-1},qquad 	au_f'=g_i	au_f.
]

Let (f_h) be a small oriented triangle with tangent edge vectors (h u,h v) at its base point. Its oriented area bivector is
[
Sigma_{f_h}=	frac12 h^2,uwedge v.
]
Use the component evaluations (R(u,v)) and (T(u,v)), with the same orientation convention as the face boundary.

For sufficiently small (h), choose the principal matrix logarithm near the identity. Define the face densities
[
mathscr R_{f_h}:=2h^{-2}log H_{f_h},qquad
mathscr T_{f_h}:=2h^{-2}	au_{f_h}.
]
The factor 2 compensates for the parameter-simplex area (1/2); equivalently one may normalize by the actual oriented area and retain the corresponding bivector components.

## 3. Local consistency proposition (conditional; proof outline)

**Proposition (smooth-field face consistency).** Under the smoothness and shape-regularity hypotheses above, and with edge data obtained by exact sampling of the same smooth pair ((e,omega)), one has in the base frame
[
log H_{f_h}=	frac12h^2R(u,v)+O(h^3),
qquad
	au_{f_h}=	frac12h^2T(u,v)+O(h^3).
]
Consequently,
[
mathscr R_{f_h}=R(u,v)+O(h),qquad
mathscr T_{f_h}=T(u,v)+O(h).
]
The constants in the remainders are uniform on compact subsets when the required derivatives of (e) and (omega) are uniformly bounded and the triangulations are uniformly shape regular.

**Proof sketch.** In a local gauge centered at the face base point, expand the path-ordered edge transports through second order. The antisymmetric second-order term in the ordered product is the curvature (domega+omegawedgeomega); the logarithm removes the non-invariant higher products and has leading term equal to the integrated curvature. For the affine translation component, expand the three transported edge coframe integrals. The first-order edge terms cancel around the closed boundary; the leading area term is the covariant exterior derivative (de+omegawedge e). The parameter domain of the oriented triangle has area (1/2), giving the displayed coefficient. Taylor remainders are (O(h^3)); shape regularity makes the bounds uniform. This is a local consistency statement for data sampled from a pre-existing smooth Cartan geometry, not a convergence theorem for arbitrary lattice configurations.

## 4. Relation to the exact finite Bianchi word

For every oriented tetrahedron ([0,1,2,3]), set
[
A=F_{012},quad B=F_{023},quad C=F_{031},quad
D=G_{01}F_{123}G_{01}^{-1}.
]
The finite edge-cancellation identity is
[
ABC D^{-1}=mathbf 1.
]
Writing (A=(H_A,t_A)), etc., its translational component is
[
t_A+H_At_B+H_AH_Bt_C-t_D=0.
]
This is exact at every mesh size and uses no commutativity assumption.

For a refinement family sampled from a smooth ((e,omega)), insert the face expansions of Section 3 and expand the transported factors consistently to the first non-vanishing order on each tetrahedron. The resulting leading continuum identities are the Cartan Bianchi identities
[
D_omega R=0,qquad D_omega T=Rwedge e.
]
The implication is a **consistency limit for smooth sampled fields**. It is not a derivation of the continuum equations from Q6 alone, and it does not establish convergence of discrete solutions. In particular, the translational identity is curvature-transported; replacing it by an untransported sum of four face vectors is incorrect.

## 5. What this result does and does not establish

Established exactly or conditionally:
1. the affine edge variables have the stated local Lorentz gauge covariance (finite algebra);
2. the face Lorentz logarithm and affine translation have the displayed curvature/torsion leading terms under smooth sampling, subject to completing the detailed uniform remainder estimates;
3. the exact tetrahedral boundary word is a finite identity and is compatible with the continuum Cartan Bianchi identities under smooth sampling/refinement.

The local consistency statement is not promoted to PROVED until the coordinate-level remainder bounds and convention checks are written out in full.

Not established:
1. existence of a shape-regular 4D Lorentzian refinement family built from the present Q6 incidence data;
2. uniform nondegeneracy and causal character of a reconstructed coframe;
3. convergence or compactness for arbitrary discrete ((L,E)), or for minimizers/critical points of a discrete action;
4. a discrete action whose independent variations yield the intended Einstein-Cartan or teleparallel field equations;
5. a compatible spin lift and convergence of the discrete matter operator to a (mathrm{Spin}^+(1,3)) Dirac operator;
6. selection of Lorentzian time from the six-dimensional Euclidean golden cut-and-project split.

## 6. Required next theorem: stability and solution convergence

A genuine discrete-to-continuum theorem must add, rather than assume, the following:
- **mesh geometry:** a specified 4D simplicial or cellular complex family with uniform shape regularity and controlled Lorentzian metric signature;
- **field bounds:** gauge-invariant bounds on holonomies, coframes, and discrete covariant differences sufficient for compactness;
- **discrete action:** a gauge-invariant action (S_h[L,E,	ext{matter}]), with independent (L)- and (E)-variations and explicit boundary terms;
- **consistency:** convergence of the discrete Euler-Lagrange operators on smooth test fields;
- **stability:** coercivity or an appropriate constrained estimate after gauge fixing, with a separate treatment of Lorentzian hyperbolicity;
- **compactness and identification:** subsequential convergence and identification of every limit as a weak solution of the target Cartan equations.

Until these gates are proved, the correct status is **local consistency proved under smooth sampling; global continuum convergence open**.

## 7. Programme connections

This refinement statement is an interface result, not a replacement for the ten-family bridge:
- Q6 supplies a finite state/graph carrier, but its graph alone is not a 4D triangulation.
- GOLDEN supplies the Euclidean (3+3) cut-and-project frame; it does not by itself select a Lorentzian time direction.
- IVM/CELL must supply an explicit incidence construction before face and 4-cell observables can be identified with physical cells.
- PSU–Anti-Twin sheet exchange remains an internal two-sheet operation, not Cartan torsion.
- SPACETIME requires the nondegenerate Lorentzian coframe and independent connection.
- REP requires a compatible (mathrm{Spin}^+(1,3)) lift; (mathrm{Spin}(6)) is not identified with it.
- ANALYSIS supplies FEEC/DEC-type consistency, stability, and compactness tools.
- DYNAMICS requires the discrete action and its independent variations.
- ARITH and MORPH remain coupled through the already recorded theorem graph, but topical graph edges are not proof dependencies.

## References and methodological baseline

The construction is framed against discrete exterior calculus, finite element exterior calculus, Regge-type discretization, and lattice gauge holonomy. Those literatures provide methods and comparison points; they do not prove the RGS/Q6-specific refinement gates above. The existing programme bibliography remains in `research/TEN_FAMILY_DISCRETE_DIFFERENTIAL_GEOMETRY_PHYSICS_BRIDGE_2026-11-27.md`.


## 8. Manufactured-solution convergence audit

A reproducible numerical audit is implemented in `scripts/audit_cartan_refinement_consistency.py`, with focused tests in `tests/test_cartan_refinement_consistency.py`.

The manufactured smooth field is defined on a 2D coordinate plane inside the 4D Lorentzian model:
[
omega_x=yB,qquad omega_y=xA,qquad
e_x=(1,0,0,0)^T,qquad e_y=(0,1,0,0)^T+x,c,
]
where (A=0.7K_x), (B=0.4K_y) are noncommuting Lorentz boost generators and (c=(0.3,0.2,-0.1,0.4)^T). At the origin,
[
R_{xy}=A-B,qquad T_{xy}=c.
]
Every oriented edge affine transport is integrated from (dG/dt=Gmathcal A_t) by RK4 with 64 steps; the triangular face is ([0,(h,0),(0,h)]). The matrix logarithm is evaluated by its convergent power series near identity.

The observed face errors are:

| (h) | curvature face error | torsion face error | observed order R | observed order T |
|---:|---:|---:|---:|---:|
| 0.2000 | 5.279724e-05 | 7.7047789e-04 | — | — |
| 0.1000 | 3.299832e-06 | 9.8119266e-05 | 4.000 | 2.973 |
| 0.0500 | 2.062395e-07 | 1.2381329e-05 | 4.000 | 2.986 |
| 0.0250 | 1.288997e-08 | 1.5550480e-06 | 4.000 | 2.993 |
| 0.0125 | 8.056230e-10 | 1.9484565e-07 | 4.000 | 2.997 |

The curvature order four is a cancellation specific to this manufactured field; the generic local estimate targeted by the proposition is (O(h^3)) for the unnormalized face observable. The torsion error exhibits the expected (O(h^3)). After division by the face area, the corresponding density errors scale as (O(h^2)) for this curvature example and (O(h)) for torsion.

Reproduce:
```bash
PYTHONPATH=. python scripts/audit_cartan_refinement_consistency.py
PYTHONPATH=. python -m pytest tests/test_cartan_refinement_consistency.py -q
```

This experiment supports the local smooth-sampling consistency claim for one noncommuting manufactured Cartan field. It does not prove uniform estimates for all smooth fields, stability/compactness, convergence of discrete solutions, or a Q6-derived Lorentzian refinement.
