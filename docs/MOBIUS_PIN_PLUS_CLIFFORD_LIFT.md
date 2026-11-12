# Finite Möbius Pin+(4) lift and Clifford transport

## Implemented object

The six-face Möbius dual graph carries the existing orientation transition signs
\(t_{ij}\in\{\pm1\}\), with product \(-1\) around
\(0\to1\to2\to5\to4\to3\to0\).

Use the Euclidean Clifford algebra \(\mathrm{Cl}(4,0)\) with Hermitian
generators \(\gamma_a\gamma_b+\gamma_b\gamma_a=2\delta_{ab}\).
For the chosen reflected axis \(a\), lift the nontrivial sign by
\[
\widehat{-1}=\gamma_a,\qquad \widehat{+1}=1.
\]
Because \(\gamma_a^2=+1\), this is a homomorphism
\(\mathbb Z_2\to\mathrm{Pin}^+(4)\). Its twisted adjoint is
\[
v\longmapsto-\gamma_a v\gamma_a^{-1},
\]
which is exactly the orthogonal reflection that reverses axis \(a\).

## Clifford-algebra transport

Let \(\Gamma=\gamma_0\gamma_1\gamma_2\gamma_3\) be the volume element.
In dimension four it implements the grade automorphism:
\(\alpha(X)=\Gamma X\Gamma^{-1}\), and \(\alpha(v)=-v\) for vectors.
For a negative transition define
\[
\Phi_{-1}(X)=\gamma_a\alpha(X)\gamma_a^{-1};
\qquad \Phi_{+1}(X)=X.
\]
Then on degree-one Clifford generators,
\[
\Phi_{-1}(v)=-\gamma_a v\gamma_a^{-1}=R_a v.
\]
Since \(\Phi_{-1}\) is a composition of algebra automorphisms, it preserves
Clifford products. Also \(\Phi_{-1}^2=\mathrm{id}\), so the induced Clifford
transport is a representation of the transition group.

## Cocycle checks

The implementation checks:
- \(\widehat{-1}^{\,2}=1\), the Pin+ relation;
- twisted-adjoint projection to the selected \(O(4)\) reflection;
- grade automorphism and Clifford-vector covariance;
- lifted local-frame covariance
  \(\hat h_j\hat t_{ij}\hat h_i^{-1}=\widehat{h_jt_{ij}h_i}\);
- the ordered dual-cycle lift equals the lift of its base holonomy \(-1\);
- the induced Clifford maps obey the \(\mathbb Z_2\) composition law.

The six-edge cycle has five identity transitions and one reflection transition,
so its lifted holonomy is \(\gamma_a\), not the identity. This is the expected
lift of the nontrivial orientation holonomy, not a cocycle failure.

## Scope

This is a Pin+ lift of the **specified finite orientation-transition system on
the six-face dual graph**, together with its Clifford-algebra transport. It is
not yet a proof of a Pin+ structure on a separately specified tangent bundle
over every cell of the 2-complex: that requires a full bundle atlas and all
overlap/cell-extension data. It does not identify the orientation cocycle
with the independent SU(2) gauge connection, and it does not imply torsion.

## Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_pin_mobius_lift.py tests/test_mobius_representations.py tests/test_discrete_spin_geometry.py -q
\`\`\`
