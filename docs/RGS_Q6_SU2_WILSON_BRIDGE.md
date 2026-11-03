# RGS–Q6–SU(2)–Wilson bridge

## 1. The missing map and the obstruction

The recovered point set has a faithful coordinate action of the signed-permutation group
\[
O_h=\{S\in O(3):S\text{ has one signed nonzero entry in each row and column}\},
\qquad |O_h|=48.
\]
The canonical Q6 square complex has 64 vertices, 192 edges and 240 square 2-cells.
For each induced vertex permutation \(\pi_S\), test separately whether it preserves
the edge set and the square-cell set.

The direct finite audit gives:
\[
|H|=16,\qquad H=\{S\in O_h:\pi_S(E)=E,\ \pi_S(F_2)=F_2\}.
\]
Thus the full \(O_h\) action on the recovered coordinates does **not** descend to
automorphisms of the independently labelled Q6 2-complex. This is a genuine
compatibility obstruction, not a numerical failure of the coordinate action.
The proper subgroup has order
\[
|H^+|=|H\cap SO(3)|=8.
\]
The other eight elements of \(H\) are improper transformations. They act on the
base 2-complex, but do not have lifts in \(SU(2)\); a spinorial extension would
require a Pin-type group and its representation.

## 2. The connection configuration space

Choose one canonical orientation for each undirected Q6 edge. A lattice connection is
\[
\mathcal A_{Q_6}=SU(2)^{E},\qquad |E|=192.
\]
There are 192 independent link matrices \(U_{ij}\in SU(2)\), with reverse links
defined by \(U_{ji}=U_{ij}^{-1}=U_{ij}^{\dagger}\). The vertex gauge group is
\[
\mathcal G_{Q_6}=SU(2)^{V},\qquad |V|=64,
\]
and, in the repository's column-vector convention,
\[
U_{ij}\mapsto U'_{ij}=G_jU_{ij}G_i^\dagger.
\]
This is an exact definition of the finite lattice gauge configuration space. The
values of the 192 links are connection data; they are not determined by the 64
coordinates alone. The current audit uses a seeded generic connection to test the
maps and invariants, not to claim a uniquely reconstructed physical connection.

## 3. The 240 plaquette holonomies

For the oriented face \(p=(v_0,v_1,v_2,v_3)\), define
\[
U_p=U_{v_3v_0}U_{v_2v_3}U_{v_1v_2}U_{v_0v_1}.
\]
The product order follows the repository path-transport convention. Under a gauge
transformation,
\[
U_p\mapsto G_{v_0}U_pG_{v_0}^{\dagger}.
\]
Therefore the fundamental trace is gauge invariant. There are exactly 240 such
plaquettes, one for each Q6 square 2-cell.

## 4. Wilson action

For the fundamental representation of \(SU(2)\), define
\[
S_W[U]=\beta\sum_{p\in F_2}
\left(1-\frac12\operatorname{Re}\operatorname{Tr}U_p\right),
\qquad \beta\ge0.
\]
This action is invariant under the local vertex gauge group. It is also invariant
under every automorphism of the Q6 2-complex: the base symmetry permutes plaquettes,
and orientation reversal replaces a holonomy by its inverse, which has the same
real trace in \(SU(2)\).

For the eight proper elements of \(H^+\), a lift \(\widetilde S\in SU(2)\) is
constructed and the induced spinorial action is tested by
\[
U_{\pi_S(i)\pi_S(j)}
=\widetilde S\,U_{ij}\,\widetilde S^{-1},
\]
with inversion when the image edge orientation is reversed. The Wilson action
is invariant under this conjugation. This is a tested diagonal spin/base action,
not a derivation that the internal gauge group must equal the spatial spin group.

## 5. What is established and what remains open

**Exact finite structure**
- \(Q_6\) has 64 vertices, 192 edges and 240 square faces.
- The connection space is \(SU(2)^{192}\), with gauge group \(SU(2)^{64}\).
- The Wilson action and plaquette definitions are explicit.

**Numerically verified by the reproducible audit**
- 16 of the 48 induced coordinate symmetries preserve the Q6 2-complex.
- 8 of these 16 are proper and admit \(SU(2)\) lifts.
- 192 sampled links and 240 plaquette holonomies satisfy the \(SU(2)\) constraints.
- Local gauge covariance, Wilson gauge invariance, base-complex symmetry invariance
  and the proper-subgroup lifted action pass at the stated tolerances.

**Not established**
- A unique physical connection derived from the RGS coordinates/coframe.
- A Pin(3) lift for the improper subgroup.
- A full \(O_h\)-equivariant SU(2) theory on the fixed Q6 labelling.
- A continuum Yang–Mills limit, a physical value of \(\beta\), or a relation to
  Standard Model gauge dynamics.

## 6. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_su2_wilson.py -q
PYTHONPATH=. python -m src.rgs_q6_su2_wilson --seed 20261103 --beta 1.0
\`\`\`

The JSON report records the symmetry subgroup orders, dimensions, gauge covariance
residuals, plaquette residuals and Wilson-action invariance residuals.
