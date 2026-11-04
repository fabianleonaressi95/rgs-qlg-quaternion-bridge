# SU(2) links and the Q6 graph Dirac-Kähler operator

## Link convention

For each canonical edge (e=(a,b)), (a<b), store (U_{ab}\in SU(2))
as transport from the fibre at (b) to the fibre at (a). Reverse traversal
uses (U_{ba}=U_{ab}^\dagger). A vertex gauge transformation acts as
(\psi_v\mapsto G_v\psi_v) and
(U_{ab}\mapsto G_aU_{ab}G_b^\dagger).

The covariant coboundary is
[
(D_U\psi)_{ab}=U_{ab}\psi_b-\psi_a.
]
It transforms at the edge tail, so
(D_{U'}=P_1D_UP_0^\dagger). The Hermitian block operator
[
\mathcal D_{DK}=\begin{pmatrix}0&D_U^\dagger\\D_U&0\end{pmatrix}
]
therefore obeys
(\mathcal D'_{DK}=P\mathcal D_{DK}P^\dagger).

## Curvature observable

For every oriented square, multiply the four forward/reverse links in
boundary order to obtain (W_f\in SU(2)). Under gauge transformations
(W_f\mapsto G_{v_0}W_fG_{v_0}^\dagger), so
\(\mathrm{tr}(W_f)\) is gauge invariant. This is the finite plaquette
holonomy, not yet a continuum field-strength tensor.

## Scope

The operator is a graph Dirac-Kähler incidence operator on vertex and edge
cochains with an internal SU(2) fibre. It is **not** a spinorial Clifford
Dirac operator, and SU(2) here is an internal gauge group, not the Lorentz
or Weitzenböck connection. No fermion doubling cure, physical QCD action,
continuum limit or TEGR coupling is claimed.
