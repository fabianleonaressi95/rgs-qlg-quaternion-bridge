# RGS–PSU golden cut-and-project: theory and finite application

## 1. Exact algebraic split

Let
[
phi=rac{1+sqrt5}{2},qquad
phi^star=rac{1-sqrt5}{2}=-phi^{-1}.
]
Use the recovered six-axis icosahedral matrix
[
A(p)=
egin{pmatrix}
0&0&1&-1&p&p\
1&-1&p&p&0&0\
p&p&0&0&1&-1
end{pmatrix}.
]
It obeys
[
A(p)A(p)^T=2(2+p)I_3,qquad
A(phi)A(phi^star)^T=0.
]
Hence
[
P_parallel=rac{A(phi)}{sqrt{2(2+phi)}},qquad
P_perp=rac{A(phi^star)}{sqrt{2(2+phi^star)}}.
]
The stacked matrix (O=(P_parallel^T,P_perp^T)^T) is orthogonal:
[
P_parallel P_parallel^T=P_perp P_perp^T=I_3,quad
P_parallel P_perp^T=0,quad OO^T=I_6.
]
For (ninmathbb Z^6),
[
pi_parallel(n)=P_parallel n,qquad
n^star=pi_perp(n)=P_perp n.
]
This is the precise 3+3 cut-and-project decomposition; it replaces the old
simplified script, which incorrectly cut on the norm of the physical
projection rather than on the internal coordinate.

## 2. Acceptance window

Let (C=[-1/2,1/2]^6). The canonical window for this projection convention is
[
W=P_perp C
=left{sum_{j=1}^6 t_j(P_perp e_j): |t_j|le	frac12ight}.
]
It is a 3D zonotope. The finite convex-hull audit gives 32 vertices and 30
coplanar rhombic facets: the rhombic triacontahedron. Its boundary is a
finite union of planar polygons and therefore has 3D Lebesgue measure zero.

The earlier spherical window of radius 1.5 is not the canonical window and
is retained only as a historical exploratory choice. A scaled family
(sW) is available for sensitivity tests; (s=1) is the canonical value.

## 3. Regular model set and finite patch

The infinite cut-and-project set is
[
Lambda(W)={P_parallel n:ninmathbb Z^6, P_perp nin W}.
]
A finite physical patch is obtained by adding (|P_parallel n|le R).
This is not a second internal cut. It is the observation-region truncation.

The integer search is finite and certified. Since (O) is orthogonal,
[
|n|^2=|P_parallel n|^2+|P_perp n|^2.
]
Writing (ho_W=max_{win W}|w|), every accepted point in the radius-(R)
patch satisfies
[
|n_j|le|n|lesqrt{R^2+ho_W^2}.
]
The implementation enumerates the resulting integer box and then applies
both exact geometric cuts. The regular-model-set hypotheses are: compact
window with nonempty interior, boundary of measure zero, and the irrational
3+3 embedding. Here the absence of integer vectors in either projection
kernel follows directly from the golden irrationality. For example,
\(A(\\phi)n=0\), with \(n\\in\\mathbb Z^6\), separates into integer and
\(\\phi\)-coefficients:
\[
n_2-n_3=0,\quad n_4+n_5=0,\quad
n_0-n_1=0,\quad n_2+n_3=0,\quad
n_4-n_5=0,\quad n_0+n_1=0,
\]
which forces \(n=0\). The same argument applies to \(A(\\phi^\\star)\),
since \(\\phi^\\star\) is also irrational and satisfies
\((\\phi^\\star)^2=\\phi^\\star+1\). Thus both projections are injective
on the lattice. Moreover, the internal projection of \(\\mathbb Z^6\) is
dense: a nontrivial annihilating character would give a nonzero integer
vector in \(E_\\perp\), contradicting the same kernel argument. The
numerical finite audit does not replace a proof of all dynamical
diffraction consequences.

## 4. Applying the cut to the 64 Q6 states

The canonical Q6 state (gin{0,1}^6) is centered as
(x_g=g-(1/2,ldots,1/2)). Its two images are
[
X_g=P_parallel x_g,qquad Y_g=P_perp x_g.
]
All 64 (Y_g) lie in (W=P_perp C) by construction. The 64 physical
images are a finite projected hypercube seed, not the whole infinite model
set. The module also compares its unlabelled pair-distance spectrum with
the hash-reconstructed RGS 64-point cloud; no equality is assumed.

## 5. Internal time direction from torsion coupling

The previously specified RGS criterion is retained. Given the 64 internal
coordinates (Yinmathbb R^{64	imes3}) and an independently derived
torsion observable (Tinmathbb R^{64	imes m}), center and normalize each
column to unit Euclidean norm, obtaining (Y_z,T_z). Form
[
C_{	au T}=rac1{64}Y_z^T T_z=USigma V^T.
]
The candidate internal temporal direction is (u_	au=U_{:,1}), and the
coupled torsion mode is (v_	au=V_{1,:}). The scalar coordinate is
[
	au_g=(Y_g-ar Y)cdot u_	au.
]
The leading singular value is the maximal normalized cross-coupling over
unit directions in the two spaces. A row-permutation null and a row-bootstrap
angular stability analysis are mandatory before retaining the candidate.

This direction is internal to (E_perp); it is not an externally appended
time coordinate. Nevertheless it is a data-selected candidate, not a proof
that physical time emerges from torsion. The required (64	imes m) torsion
trace is deliberately an explicit input: the archived 192x4 edge coframe is
not the same observable and is not silently substituted.

## 6. The 4D bridge

After a stable (u_	au) is selected, split
[
E_perp=operatorname{span}(u_	au)oplus E_{perp,mathrm{space}},
qquad E_{perp,mathrm{space}}=u_	au^perp.
]
A candidate 4D chart is ((	au,X)inmathbb R	imes E_parallel).
A Lorentzian metric must still be specified and tested independently; the
Euclidean cut-and-project orthogonality does not itself produce signature
((-+++)), a vierbein, TEGR torsion, or causal dynamics.

## 7. Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_psu_cut_project.py -q
PYTHONPATH=. python -m src.rgs_psu_cut_project --output results/rgs_psu_cut_project --physical-radius 2.5
```

When a verified torsion trace is available:
```bash
PYTHONPATH=. python -m src.rgs_psu_cut_project --output results/rgs_psu_cut_project --rgs-coordinates PATH_TO_X_OH.npy --torsion-trace PATH_TO_TORSION_TRACE.npy
```

## 8. Evidence ledger

- **ANALYTIC IDENTITY:** golden conjugacy and the stated matrix Gram/cross
  identities.
- **NUMERICALLY VERIFIED:** orthogonal 3+3 frame, projected-cube window
  (32 vertices/30 facets), Q6 projections, finite-radius model-set patch.
- **TO BE TESTED ON ARCHIVE:** congruence/similarity of the cut-projected Q6
  seed and the reconstructed RGS 64-point coordinates.
- **CANDIDATE / DATA-DEPENDENT:** torsion-coupled internal time axis; requires
  the independent torsion trace, permutation significance and bootstrap
  stability.
- **OPEN:** physical interpretation, 4D Lorentzian metric, tetrad/TEGR,
  continuum limits, diffraction/transport predictions and experimental test.
