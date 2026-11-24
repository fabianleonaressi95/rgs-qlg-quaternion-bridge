# Q6 metric-dependent plaquette response from the recovered RGS embedding

## Geometric input: no invented Q6-to-R3 map

The implementation reuses `load_rgs_seed_data()` and
`build_rgs_q6_induced_metric()`: the archived 64-point map
(X:Q_6^{(0)}\to\mathbb R^3), canonical 240 faces and measured embedded
face Gram matrices. It does not interpret the 4-component edge coframe as
a vertex tetrad, nor does it assert a nondegenerate 6D metric. The recovered
pullback metric is rank at most three.

## Declared finite face-Hodge action

For a canonical face (p=(v,v_a,v_{ab},v_b)), set
(a_p=X_{v_a}-X_v), (b_p=X_{v_b}-X_v) and
[
G_p=\begin{pmatrix}a_p\cdot a_p&a_p\cdot b_p\\
a_p\cdot b_p&b_p\cdot b_p\end{pmatrix},\qquad
A_p=\sqrt{\det G_p}.
]
The reference area (A_p^0) is fixed at the undeformed recovered geometry.
We define the conditional metric-dependent plaquette action
[
S_g[X,U]=\beta\sum_{p=1}^{240}\frac{A_p^0}{A_p(X)}
\left(1-\tfrac12\operatorname{ReTr}U_p\right).
]
This is the face-area part of a discrete 2-form Hodge weight with a fixed
reference dual measure. The repository seed does not specify circumcentric
dual lengths/volumes, so the expression is explicitly a *conditional
finite face-Hodge model*, not claimed to be the unique DEC Hodge star or
the continuum Yang--Mills action.

For ambient shear (X(\gamma)=(I+\gamma H)X), with
(H=(E_{xy}+E_{yx})/2), define
[
\dot G_p=
\begin{pmatrix}
2a_p\cdot Ha_p&a_p\cdot Hb_p+Ha_p\cdot b_p\\
a_p\cdot Hb_p+Ha_p\cdot b_p&2b_p\cdot Hb_p
\end{pmatrix}.
]
Jacobi's formula gives the exact identity
[
\frac{\dot A_p}{A_p}
=\frac12\operatorname{tr}(G_p^{-1}\dot G_p).
]
Holding the link holonomies and reference areas fixed,
[
\boxed{\Pi_{xy}^{g}
=-\frac1V\dot S_g
=\frac{\beta}{V}\sum_p
\frac{A_p^0}{A_p}
\left(1-\tfrac12\operatorname{ReTr}U_p\right)
\frac12\operatorname{tr}(G_p^{-1}\dot G_p).}
]
The code computes the analytic derivative and audits it against a
central finite difference.

## What this closes and what it does not

- **Exact finite identity:** face-area derivative and action derivative
  for the declared shear and fixed holonomies.
- **Numerical verification:** tests compare the analytic result with
  independent central differences on the hash-recovered RGS geometry.
- **Still open:** the dual-face measure, a complete 3D/4D discrete Hodge
  star compatible with the graph, the relation of the link variables to
  physical edge transport, equilibrium dynamics/Gibbs measure, FDT
  normalization, thermodynamic limit, Green--Kubo convergence and
  Stokes matching.

The new action is not silently merged into the existing plaquette-only
Wilson action. It is a separately named conditional metric extension
that can be promoted only after the dual measure and continuum scaling
are derived.
