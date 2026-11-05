# Integral splitting theorem for the Q6 annihilator torus

## 1. Conventions and hypotheses

Let (X) be the connected Q6 graph, with (V=64) vertices and (E=192) oriented edges. Let
[
Dinmathbb Z^{E	imes V}
]
be its vertex-edge incidence matrix, with the convention that a vertex phase
(	hetainmathbb T^V) produces the edge gradient (D	heta). Thus
(operatorname{rank}D=V-1=63), and the integer cycle lattice is
[
Z_1=ker(D^T:mathbb Z^E	omathbb Z^V).
]
Choose an integral cycle basis (Cinmathbb Z^{129	imes192}), whose rows
form a (mathbb Z)-basis of (Z_1). In particular,
[
operatorname{rank}C=129,qquad CD=0.
]

Let (P_5) be the exact rational projector onto the selected five-dimensional
coexact eigenspace, and let
[
M:mathbb Z^{129}longrightarrowmathbb Q^5,qquad
clongmapsto cCP_5
]
denote the five-mode map in any rational basis of that eigenspace. Put
[
K=ker_{mathbb Z}M.
]
The audited Q6 construction has (operatorname{rank}M=5) and
(operatorname{rank}K=124). Let (B_Kinmathbb Z^{124	imes129})
be a saturated row basis of (K), and define
[
L=B_KCinmathbb Z^{124	imes192}.
]
The exact Smith audit for this lattice has 124 unit invariant factors.
Equivalently, (Lmathbb Z^{192}) is a primitive row sublattice of rank
124. This primitivity is a hypothesis in the structural theorem below and
must be checked against the exact matrices used in a given revision.

Define the annihilator subgroup and the vertex-gauge subgroup by
[
H={xinmathbb T^E:Lx=0pmod{2pi}},qquad
G=operatorname{im}(D:mathbb T^V	omathbb T^E).
]

## 2. Theorem

**Theorem (integral gauge quotient and five-torus splitting).** Under the
hypotheses above:

1. (H) is a connected 68-dimensional torus.
2. (G) is a connected 63-dimensional subtorus of (H).
3. The quotient (H/G) is a connected five-dimensional torus.
4. Its character lattice is canonically identified with the five-mode
   integer quotient:
   [
   widehat{H/G}cong Z_1/(K C)congmathbb Z^{129}/K
   cong M(mathbb Z^{129})subsetmathbb Q^5.
   ]
   Here (KC) denotes the edge-flux image of the kernel lattice.
5. The exact sequence of compact groups
   [
   0longrightarrow Glongrightarrow Hlongrightarrow
   H/Glongrightarrow0
   ]
   splits noncanonically as a sequence of tori. Equivalently, there is a
   unimodular integral basis of the cocharacter lattice of (H) whose first
   63 vectors span the gauge sublattice and whose remaining five vectors
   represent a complementary quotient lattice.

The conclusion is an abstract integral splitting. It does not identify the
five quotient coordinates with a preferred orthonormal basis of the
coexact eigenspace; that requires fixing and auditing the rational coordinate
map (M).

## 3. Proof

**Step 1: (H) is a connected 68-torus.** The homomorphism
(L:mathbb T^{192}	omathbb T^{124}) has kernel (H). Its character
lattice is
[
widehat H=mathbb Z^{192}/operatorname{row}_{mathbb Z}(L).
]
Because the row lattice of (L) is primitive of rank 124, this quotient
is free abelian of rank (192-124=68). Pontryagin duality therefore gives
(Hcongmathbb T^{68}), with no finite component group.

**Step 2: the gauge subgroup lies in (H).** By construction,
[
LD=B_KCD=0.
]
Consequently (L(D	heta)=0) modulo (2pi) for every
(	hetainmathbb T^{64}), so (Gsubset H). Since the graph is connected,
the kernel of (D:mathbb T^{64}	omathbb T^{192}) consists precisely of
the constant vertex phases. Hence (Gcongmathbb T^{63}) and is connected.

**Step 3: the quotient is a five-torus.** A quotient of a compact connected
abelian Lie group by a connected closed subgroup is again a compact
connected abelian Lie group. Its dimension is
[
dim(H/G)=68-63=5.
]
Every compact connected abelian Lie group of dimension five is isomorphic
to (mathbb T^5). This proves the group-level dimension statement without
confusing graph cycle rank with the first Betti number of the full cubical
complex.

**Step 4: identify the character lattice.** A character of (H) is
represented by an integer edge vector (ninmathbb Z^{192}), modulo the
row lattice (L). It descends to (H/G) exactly when it is trivial on
every gauge gradient:
[
n^TD	heta=0pmod{2pi}quadorall	heta
quadLongleftrightarrowquad D^Tn=0.
]
Thus
[
widehat{H/G}
=ker(D^T:mathbb Z^{192}	omathbb Z^{64})/
  operatorname{row}_{mathbb Z}(L)
=Z_1/(KC).
]
Since the rows of (C) are a (mathbb Z)-basis of (Z_1), the map
(cmapsto cC) identifies (Z_1) with (mathbb Z^{129}), and sends
(K) exactly to (KC). Therefore
[
Z_1/(KC)congmathbb Z^{129}/K.
]
By the first isomorphism theorem for the integer mode map,
(mathbb Z^{129}/Kcong M(mathbb Z^{129})). This proves that the
five-torus quotient has precisely the character lattice of the five-mode
integer quotient, not merely the same dimension.

**Step 5: integral splitting.** The character-lattice sequence dual to
(0	o G	o H	o H/G	o0) is an exact sequence of free abelian groups
[
0longrightarrowwidehat{H/G}
longrightarrowwidehat H
longrightarrowwidehat G
longrightarrow0.
]
The quotient lattice (widehat Gcongmathbb Z^{63}) is free, so this
sequence splits as a sequence of abelian groups. Dualizing gives a
noncanonical splitting of the torus sequence. Equivalently, the
cocharacter lattice of (G) is a primitive rank-63 sublattice of the
rank-68 cocharacter lattice of (H), and a basis of it can be extended
to a unimodular basis of the latter. (square)

## 4. Explicit unimodular construction

The proof gives a finite exact algorithm for a basis, not only an existence
claim.

1. Compute a Smith decomposition (U_L L V_L=S_L), with
   (U_L,V_L) unimodular. Since the nonzero invariant factors are all one,
   the last 68 columns of (V_L) form an integral basis
   (Q_Hinmathbb Z^{192	imes68}) of (ker_{mathbb Z}L).
2. Solve the exact integer system
   [
   Q_H X_G=D,qquad X_Ginmathbb Z^{68	imes64}.
   ]
   Existence follows from (Gsubset H); the rank is 63.
3. Compute a Smith decomposition
   [
   U_G X_G V_G=S_G.
   ]
   Since the gauge subtorus is connected, its cocharacter lattice is
   primitive in that of (H); hence all 63 nonzero invariant factors of
   (S_G) must equal one.
4. The columns of (U_G^{-1}) form a unimodular basis in (H)-coordinates.
   Its first 63 columns span the gauge lattice; its last five columns give
   an integral complement. Mapping these columns through (Q_H) gives the
   requested 192-dimensional edge-coordinate basis.

The implementation must verify, using exact integer arithmetic:
[
LQ_H=0,quad operatorname{rank}Q_H=68,quad
Q_HX_G=D,quad
operatorname{SNF}(X_G)=(1,ldots,1,0),
]
and that the resulting 68-column basis is unimodular in the chosen
cocharacter coordinates. Floating-point nullspaces are not acceptable for
these checks.

## 5. Haar measure and Wilson interpretation

The quotient map (p:H	o H/Gcongmathbb T^5) pushes normalized Haar
measure on (H) to normalized Haar measure on (mathbb T^5). The
disintegration is Haar measure on the gauge fibres, with no Jacobian
factor because this is a homomorphism of compact groups with normalized
Haar measures.

This statement concerns Haar measure only. The Wilson integrand used in
the annihilator representation is generally not gauge invariant before
the appropriate character/constraint construction. Therefore one may not
replace an integral over (H) by an integral over (H/G) for that
integrand without first specifying its gauge average or conditional fibre
integral. The exact character-lattice identification above does, however,
give the correct five-dimensional *character group* on which the
fibre-marginal Wilson weight is indexed.

It still does not prove that the fibre-marginal Wilson law is absolutely
continuous on a cotangent chart, equals a Liouville density, or produces
the D5 Hua factor (N(z)^{-5}). Those are separate measure-matching
questions.

## 6. Notation correction required in the global note

With (Dinmathbb Z^{192	imes64}), the gradient map is
(D:mathbb T^{64}	omathbb T^{192}), the cycle condition is
(D^Tn=0), and the matrix identity is (CD=0). Thus the gauge subgroup
must be written (operatorname{im}D), not (operatorname{im}D^T);
the compatible identity is (LD=B_KCD=0). Any expression such as
(CD^T) is dimensionally invalid with these conventions.

## 7. Evidence status

- **PROVED under the stated exact lattice hypotheses:** connectedness and
  dimension of (H); inclusion and dimension of (G); abstract
  (H/Gcongmathbb T^5); canonical character-lattice identification
  with (mathbb Z^{129}/K); existence of an integral splitting.
- **EXACT REPOSITORY AUDIT already reported:** rank-124 saturated kernel,
  rank-124 image lattice with unit Smith invariants, rank-five mode map.
  These facts must be rerun against the same commit as any explicit basis.
- **NOT YET COMPUTED in this note:** the explicit 68-by-68 unimodular
  complement matrix and its machine-readable certificate.
- **OPEN:** compatibility with a preferred orthonormal D5 frame, the
  conditional Wilson density on the quotient, the cotangent lift, and
  the Wilson-to-Hua Radon–Nikodym identity.

## References

- S. A. Morris, *Pontryagin Duality and the Structure of Locally Compact
  Abelian Groups*, Cambridge University Press (1977).
- T. Tao, “Using the Smith normal form to manipulate lattice subgroups
  and closed torus subgroups” (2022).
- Standard Smith/Hermite normal form theory over (mathbb Z), as
  implemented by FLINT: https://flintlib.org/doc/algorithms/linear-algebra/hermite-smith.html.
