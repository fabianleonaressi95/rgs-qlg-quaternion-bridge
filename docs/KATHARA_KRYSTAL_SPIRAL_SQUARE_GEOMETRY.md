# Kathara / Krystal Spiral: exact square-similarity reconstruction

**Status:** exact algebraic reconstruction of the stated planar stage rule, conditional on the source's geometric description. It is not a claim that the associated cosmological or energetic interpretations are established physics.

## 1. Source data used

The source description states:

1. the Kathara diagram has 12 numbered centers and a central AzurA point;
2. the central column from center 1 to center 12 and the width from center 6 to center 7 are in ratio \(2:1\);
3. successive grids rotate clockwise by \(45^\circ\) and expand by \(\sqrt2\);
4. the new center 9 overlays the previous grid's center 12;
5. after two stages, the new center 6 overlays the center 12 from two stages earlier.

Source: A&A Deane / Emerald Awakening, “Spirals of Creation”, especially the Kathara description and stage rule:
https://emerald24.org/spirals-of-creation/
The source page identifies the Kathara diagram and shows the 12-center geometry and stage overlays in its two diagrams. This file formalizes those stated geometric constraints; it does not endorse the source's non-standard physical claims.

## 2. The square similarity

Let \(R_{-\pi/4}\) be clockwise rotation by \(45^\circ\). A square diagonal has length \(\sqrt2\) times its side, so the combined rotation and scale map is
\[
T=\sqrt2 R_{-\pi/4}
=\begin{pmatrix}1&1\\-1&1\end{pmatrix}.
\]
Its columns are the two orthogonal square diagonals \((1,-1)\) and \((1,1)\). Therefore
\[
T^TT=2I,\qquad \det T=2.
\]
Two stages give
\[
T^2=\begin{pmatrix}0&2\\-2&0\end{pmatrix}
=2R_{-\pi/2}.
\]
Thus two stages are exactly a clockwise quarter-turn and a doubling of length, not merely an approximate spiral fit.

In complex notation, \(T\) is multiplication by \(1-i\):
\[
z\longmapsto(1-i)z,\qquad
(1-i)^2=-2i.
\]

## 3. Recovering the center-6 / center-9 / center-12 relation

Put AzurA at the origin and let \(w>0\) be the distance between centers 6 and 7. The midpoint condition gives
\[
p_6=(-w/2,0),\qquad p_7=(w/2,0).
\]
The \(2:1\) proportion and the statement that AzurA is the midpoint of the central 1-to-12 column give
\[
p_1=(0,-w),\qquad p_{12}=(0,w).
\]

Translate the two source overlay statements into equations. “Stage-1 center 9 overlays stage-0 center 12” means
\[
Tp_9=p_{12}.
\]
“Stage-2 center 6 overlays stage-0 center 12” means
\[
T^2p_6=p_{12}.
\]
Since \(T\) is invertible, the two equations imply
\[
\boxed{p_9=Tp_6=(-w/2,w/2).}
\]
A direct check gives
\[
Tp_9=(0,w)=p_{12},\qquad T^2p_6=(0,w)=p_{12}.
\]

The three radii are therefore
\[
\|p_6\|=w/2,\qquad
\|p_9\|=w/\sqrt2,\qquad
\|p_{12}\|=w,
\]
so each stage multiplies the radius by exactly \(\sqrt2\), while rotating clockwise by \(45^\circ\).

This is a genuine three-center geometric relation extracted from the stated overlay constraints. It does not yet determine the coordinates of all remaining Kathara centers: their exact positions require the complete line/center diagram and its incidence convention.

The three vectors also obey the exact second-order relation
\[
p_{12}-2p_9+2p_6=0,
\]
or \(p_{n+2}=2p_{n+1}-2p_n\) along this three-point orbit. Its characteristic polynomial is \(\lambda^2-2\lambda+2\), not the supergolden cubic. The three-center spiral orbit alone therefore cannot be the desired tribonacci refinement.


## 4. Gaussian lattice consequence

For \((x,y)\in\mathbb Z^2\),
\[
T(x,y)=(x+y,-x+y).
\]
The two output coordinates have the same parity. Conversely, if \(u\equiv v\pmod2\), then
\[
x=(u-v)/2,\qquad y=(u+v)/2
\]
are integers. Hence
\[
T(\mathbb Z^2)=\{(u,v)\in\mathbb Z^2:u\equiv v\pmod2\},
\]
an index-two checkerboard sublattice. Also \(T^2(\mathbb Z^2)=2\mathbb Z^2\), of index four.

This is the exact arithmetic content of the square-diagonal similarity. It is not a three-state tribonacci transfer.

## 5. Relation to the trigram/Krist-All programme

The map \(T\) acts on the planar geometric coordinates. The Walsh matrix \(E\) acts on functions on the eight trigram labels. These are different representation spaces. A bridge must specify a geometric labeling map and a compatible action on oriented cells before a compression such as \(E^T\mathcal R E\) is meaningful.

The source-derived square map has characteristic polynomial
\[
\lambda^2-2\lambda+2,
\]
not \(\lambda^3-\lambda^2-1\). Therefore the \(45^\circ/\sqrt2\) stage rule alone does not derive the supergolden constant. The remaining task is to reconstruct the complete oriented Kathara cell/line incidence and then derive its induced action on the trigram channels. No tribonacci matrix is inserted by hand.

## 6. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_kathara_square_similarity.py -q
PYTHONPATH=. python examples/kathara_square_similarity_audit.py
\`\`\`
