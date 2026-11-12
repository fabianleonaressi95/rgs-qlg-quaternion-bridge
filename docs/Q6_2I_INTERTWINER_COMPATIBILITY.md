# Q6–2I intertwiner compatibility

## Result

The icosahedral signed-axis action on the six Q6 coordinate directions is an
\(A_5\)-representation. Pull it back along
\[
2I\longrightarrow 2I/\{\pm1\}\cong A_5.
\]
The central element \(-1\in2I\) therefore acts trivially on the six-dimensional
Q6 direction space and on its induced exterior module
\(\Lambda^*(\mathbb C^6)\), of dimension 64.

By contrast, the defining fundamental representation has
\[
\rho_{\mathbf2}(-1)=-I_2.
\]
For an intertwiner \(T:\mathbf2\to\Lambda^*(\mathbb C^6)\), equivariance at
the central element would require
\[
T(-I_2)=I_{64}T,
\]
so \(T=-T\), hence \(T=0\). Thus the untwisted Hom space is exactly zero.
This is a representation-theoretic obstruction, not a numerical near miss.

## Minimal central-character twist

Define the tensor-product representation
\[
\rho_{\mathrm{tw}}(q)=
\rho_{\mathbf2}(q)\otimes\rho_{\Lambda^*}(q)
\quad\text{on}\quad
\mathbb C^2\otimes\Lambda^*(\mathbb C^6).
\]
Its dimension is 128 and its central character is now nontrivial. Since the
degree-zero vacuum \(1\in\Lambda^0(\mathbb C^6)\) is fixed by every induced
exterior action, the map
\[
J:\mathbb C^2\longrightarrow
\mathbb C^2\otimes\Lambda^*(\mathbb C^6),
\qquad J(v)=v\otimes1
\]
is a nonzero \(2I\)-equivariant injection. This proves that the twist removes
the central-character obstruction and supplies an explicit intertwiner.

It does **not** prove that this tensor factor is physically or geometrically
required. It is the minimal algebraic repair, with the degree-zero sector
providing a canonical embedding.

## Implementation and tests

Implementation: src/q6_2i_intertwiners.py
Tests: tests/test_q6_2i_intertwiners.py

The audit matches each quaternion's SO(3) rotation to the 60 signed-axis
actions, induces the exterior action on all 64 basis states, checks the
central actions, and checks the explicit twisted intertwiner against all 120
elements of \(2I\).

Reproduction command:
PYTHONPATH=. python -m pytest tests/test_q6_2i_intertwiners.py tests/test_representations_2i.py tests/test_rgs_icosahedral_q6_intertwiner.py -q

## Separation from Pin+

The preceding Pin+(4) construction lifts the specified Möbius orientation
transition graph. The present result concerns the pullback of the proper
icosahedral \(A_5\) action to \(2I\). No homomorphism identifying these two
actions, no common bundle transition atlas, and no Pin-to-\(2I\) coupling has
been constructed. The central-character calculation tells us exactly why the
untwisted Q6 exterior module cannot carry the fundamental doublet as a
subrepresentation.
