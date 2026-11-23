# Q6–2I intertwiner compatibility

## Exact representation-theoretic result

The signed-axis action of the proper icosahedral group \(A_5\) on the six Q6
directions induces an action on
\[
W=\Lambda^*(\mathbb C^6),\qquad \dim_{\mathbb C}W=64.
\]
Pull it back along the abstract quotient
\[
2I\longrightarrow 2I/\{\pm1\}\cong A_5.
\]
The central element \(-1\in2I\) then acts trivially on \(W\). In contrast, the
fundamental doublet \(S=\mathbb C^2\) has \(\rho_S(-1)=-I_2\).

For an intertwiner \(T:S\to W\), equivariance at \(-1\) gives
\[
T(-I_2)=I_{64}T,
\]
hence \(T=-T\), so \(T=0\). Therefore
\[
\operatorname{Hom}_{2I}(S,W)=0.
\]
This is an exact central-character obstruction, not a numerical near miss.

## Explicit tensor twist

On \(S\otimes W\), use the tensor-product action
\[
\rho_{\mathrm{tw}}(g)=\rho_S(g)\otimes\rho_W(g).
\]
The degree-zero vacuum \(1\in\Lambda^0(\mathbb C^6)\) is fixed by every
induced exterior action. Thus
\[
J:S\to S\otimes W,\qquad J(v)=v\otimes1
\]
is an explicit injective intertwiner. This proves the twist removes the
central-character obstruction algebraically.

It does **not** prove that the added fundamental factor is geometrically or
physically required.

## Finite implementation audit

Implementation: src/q6_2i_intertwiners.py
Tests: tests/test_q6_2i_intertwiners.py

The code induces the signed-permutation action on all 64 exterior basis states,
checks all 60 A5 signed-axis actions and all 120 fundamental 2I matrices, and
checks the tensor-product embedding using the invariant vacuum. It deliberately
does not match individual quaternion matrices to individual signed-axis
matrices: an explicit isomorphism between those two A5 matrix presentations
has not been constructed.

Reproduction:
PYTHONPATH=. python -m pytest tests/test_q6_2i_intertwiners.py tests/test_representations_2i.py tests/test_rgs_icosahedral_q6_intertwiner.py -q

## Separation from Pin+

The Pin+(4) construction lifts the specified Möbius orientation-transition
graph. This result concerns the pullback of the proper icosahedral A5 action
to 2I. No common bundle atlas, Pin-to-2I homomorphism, torsion law or QLG
interaction is asserted.
