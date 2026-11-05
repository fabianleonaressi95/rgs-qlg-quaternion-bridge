# Q6–A4–STF(3): linear representation intertwiner

## Result being implemented

The Q6 cubical Hodge construction already identifies its top coexact
\(\lambda=12\) five-mode space as the irreducible real five-dimensional
\(A_5\)-module with character \((5,1,-1,0,0)\). Restricting this module
to an icosahedral subgroup \(A_4\) gives

\[
V_5\vert_{A_4}\simeq \mathbf 1'\oplus\mathbf 1''\oplus\mathbf 3
\simeq E_{\mathbb R}\oplus T_{\mathbb R}.
\]

The tetrahedral action on symmetric traceless tensors in three dimensions
has the same decomposition:

\[
\operatorname{Sym}_0(\mathbb R^3)\simeq E_{\mathbb R}\oplus T_{\mathbb R}.
\]

The target action is \(X\mapsto R(g)XR(g)^T\). In the Frobenius-orthonormal
basis used by the implementation, the first two basis tensors are diagonal
traceless matrices and the last three are off-diagonal symmetric matrices.

## Intertwining equation

For the restricted Q6 representation \(\rho_Q\) and the STF representation
\(\rho_{\rm STF}\), solve

\[
\rho_{\rm STF}(g)T=T\rho_Q(g),\qquad g\in A_4.
\]

The real intertwiner space has dimension three: the real two-dimensional
irreducible component has a two-dimensional real commutant, while the
irreducible three-dimensional component has a one-dimensional commutant.
A generic element of this space is invertible. The implementation finds a
deterministic invertible element and takes its orthogonal polar factor.

The A4 subgroup is found algorithmically inside the 60 proper icosahedral
rotations as the normalizer of a Klein four subgroup. This avoids relying
on undocumented group-element indices.

## Audit

Run:

    pytest -q tests/test_rgs_q6_a4_stf_intertwiner.py

The audit checks subgroup order 12, source and target orthogonality,
intertwiner nullity 3, rank 5, the restricted character, and equivariance
on every element of the selected subgroup.

## Scope

This is a finite-dimensional representation-theoretic identification. It
does not prove that the five Q6 modes are physically tensor components of
a continuum metric, does not produce a global Q6-to-manifold embedding,
and does not establish a quaternionic structure or Yang–Mills dynamics.
