# GEC-6: trigram block decomposition and antipodal action

**Status:** exact representation-theoretic decomposition of the existing GEC-6 maps. The source is the Q6 graph cycle lattice, not the filled cubical 6-cube.

## 1. I-Ching block convention

Keep the repository's bit convention: coordinates 0,1,2 are the lower trigram and 3,4,5 are the upper trigram. Write
\[
M=\mathbb Z^6=L\oplus U,\quad
L=\langle e_0,e_1,e_2\rangle,\quad
U=\langle e_3,e_4,e_5\rangle.
\]
The 15 pair directions split canonically as
\[
\boxed{\Lambda^2M=
\Lambda^2L\oplus(L\otimes U)\oplus\Lambda^2U}
\]
with ranks \(3+9+3=15\). These are respectively lower-internal pairs, lower/upper mixed pairs, and upper-internal pairs.

This is a decomposition of the exterior-square target. It is not a decomposition into three independent copies of the final rank-five quotient.

## 2. Symmetry group preserving the trigram blocks

The block-preserving subgroup is \(H_0=S_3^{(L)}\times S_3^{(U)}\). Allowing exchange of lower and upper blocks gives
\[
H=(S_3\times S_3)\rtimes C_2\cong S_3\wr C_2.
\]
The 15-dimensional exterior-square module is H-invariant. Over \(\mathbb Q\), let \(\mathbf 1_L, W_L\) denote the trivial and two-dimensional standard representations of the lower S3, and similarly for U. Then
\[
\Lambda^2L\cong W_L\oplus\mathrm{sgn}_L,\qquad
\Lambda^2U\cong W_U\oplus\mathrm{sgn}_U,
\]
and
\[
L\otimes U\cong
(\mathbf1_L\otimes\mathbf1_U)
\oplus(W_L\otimes\mathbf1_U)
\oplus(\mathbf1_L\otimes W_U)
\oplus(W_L\otimes W_U),
\]
of dimensions \(1+2+2+4=9\). The block exchange interchanges the two internal-pair sectors. On the mixed sector, under the identification \(L\otimes U\to U\otimes L\), it acts as minus the tensor-factor swap because \(u\wedge l=-l\wedge u\).

## 3. Restriction of the A5 quotient

The contraction from GEC-6 is
\[
c:\Lambda^2(\mathbb Z^6)\to A_5,\qquad
c(e_i\wedge e_j)=e_j-e_i\quad(i<j).
\]
As an H0-module, the rational A5 representation restricts as
\[
A_5\otimes\mathbb Q
\cong W_L\oplus W_U\oplus\mathbf1_{-},
\]
where \(\mathbf1_-\) is the one-dimensional contrast representation, represented by vectors constant on each block with opposite block sums. The block exchange swaps \(W_L,W_U\) and acts by -1 on \(\mathbf1_-\).

The three pieces of the exterior square behave differently under c:

* \(c(\Lambda^2L)\) is the rank-two lower-trigram difference lattice \(A_2^{(L)}\); the map has rank 2 and a rank-one kernel.
* \(c(\Lambda^2U)\) is the rank-two upper-trigram difference lattice \(A_2^{(U)}\); the map has rank 2 and a rank-one kernel.
* \(c(L\otimes U)\) is already all of \(A_5\), of rank 5. In matrix form, a coefficient matrix \(B=(b_{ia})\in\mathbb Z^{3\times3}\) maps to lower coordinates minus row sums and upper coordinates plus column sums. Its kernel consists exactly of matrices with every row sum and every column sum zero, a rank-four lattice. Over \(\mathbb Q\), this kernel is \(W_L\otimes W_U\).

Consequently, the full contraction is surjective, but its kernel is **not** the direct sum of the three separate kernels just listed: contributions from the two internal sectors can be cancelled by mixed-sector contributions. The total kernel has rank 10, as already established in GEC-6.

The mixed sector alone reaches the entire rank-five quotient. This makes the nine lower/upper pairings structurally important: they carry both trigram-difference modules and the block-contrast mode.

## 4. Antipodal complement

The antipodal map on Q6 vertices is \\(a(v)=v\\oplus63\\), complementing all six bits. The orientation convention matters: graph edges are globally oriented from bit 0 to bit 1. Under a bit translation \\(v\\mapsto v\\oplus t\\), the induced orientation on a coordinate-pair square \\(Q_{ij}\\) changes by
\\[
\\boxed{\\epsilon_{ij}(t)=(-1)^{t_i+t_j}.}
\\]
Thus a general element of \\(C_2^6\\) does **not** act trivially on \\(\\Lambda^2(\\mathbb Z^6)\\). In fact, for a translation flipping only coordinate 0, the relation
\\[
r=e_0\\wedge e_1+e_1\\wedge e_2-e_0\\wedge e_2\\in\\ker c
\\]
is sent to \\(-e_0\\wedge e_1+e_1\\wedge e_2+e_0\\wedge e_2\\), whose contraction is \\(2(e_2-e_1)\\ne0\\). Therefore this general translation action does not descend to the quotient \\(A_5\\).

For the antipodal translation \\(t=(1,1,1,1,1,1)\\), however, \\(t_i+t_j=2\\) for every pair, so \\(\\epsilon_{ij}(t)=+1\\). The identity and the global antipode are precisely the two bit translations acting trivially on every pair generator. Consequently the Q6 antipodal involution acts as the identity on the GEC-6 quotient. This does **not** mean it fixes every individual Q6 vertex, edge, or graph cycle.

## 5. What this does and does not establish

| Claim | Status |
|---|---|
| Pair directions split into 3 lower, 9 mixed, 3 upper | PROVED |
| H0 and block-swap preserve the decomposition | PROVED |
| Rational H0-module decomposition above | PROVED |
| Internal sectors map onto their respective A2 difference lattices | PROVED |
| Mixed sector alone maps surjectively onto A5 | PROVED |
| Antipodal complement acts trivially on the GEC-6 quotient | PROVED |
| A preferred I-Ching interpretation of the 9 mixed pairs | OPEN |
| A physical coupling law or dynamics from this decomposition | OPEN |

The block decomposition respects the lower/upper trigram distinction but does not yet encode the six-line ordering's traditional yin/yang semantics, King Wen sequence, or a physical interaction. Those require additional structure and must not be inferred from the representation decomposition alone.

## 6. Reproduction

Run:
~~~bash
PYTHONPATH=. python -m pytest tests/test_gec6_trigram_decomposition.py -q
PYTHONPATH=. python examples/gec6_trigram_decomposition_audit.py
~~~


## 7. Icosahedral selection of a matching family

The existing six-axis icosahedral representation distinguishes a five-element orbit of perfect matchings from a ten-element orbit. The diagonal lower-to-upper matching belongs to the five-element orbit, but each 3+3 partition admits two special cross-matchings, interchanged by its six-element A5 stabilizer. The full signed A5 action does not descend to the GEC-6 contraction quotient: only six of its 60 elements preserve `ker(c)`, exactly the selected block-partition stabilizer. See [the full orbit and compatibility audit](GEC6_ICOSAHEDRAL_MATCHING_ORBITS.md).


The distinguished five-element matching orbit has five distinct stabilizers of order 12, each with element-order profile 1 identity, 3 elements of order 2 and 8 of order 3; these are the five (A_4) subgroups of (A_5). A generic matching has stabilizer order 6. This gives the five special matchings a group-theoretic tetrahedral interpretation, but still does not select one matching from the family. See [the full orbit and compatibility audit](GEC6_ICOSAHEDRAL_MATCHING_ORBITS.md).
