# G6: Mosco convergence of the flat periodic P1 Dirichlet forms

**Status:** analytic Mosco proof for the regular flat Freudenthal P1 family.  
**Not claimed:** Mosco convergence for a curved RGS-derived family.

## 1. Common-space forms

Fix the flat torus T_L^d, d=2 or 3, and the periodic Freudenthal
triangulations K_h with h=L/N and N tending to infinity. Let V_h be
the continuous periodic P1 finite-element space and let
J_h: R^{V_h} -> L2(T_L^d) be the G4 nodal reconstruction. Define
[
mathcal Q_h(v)=
egin{cases}
displaystyleint_{mathbb T_L^d}|
abla v_h|^2,dx,
&v=J_hu_h=v_hin V_h,\
+infty,&v
otin V_h.
end{cases}
]
The integral is exactly u_h^* K_h u_h; the Hilbert norm is exactly
u_h^* M_h u_h. The candidate limit is
[
mathcal Q(v)=
egin{cases}
displaystyleint_{mathbb T_L^d}|
abla v|^2,dx,
&vin H^1(mathbb T_L^d),\
+infty,&v
otin H^1(mathbb T_L^d).
end{cases}
]
This is the P1 scalar Dirichlet form. It is distinct from the
Fourier-embedded Wilson form in MOSCO_PROTOTYPE_REGULAR_TORUS.md,
which has a different operator symbol and normalization.

## 2. M1: weak liminf

Let v_h converge weakly to v in L2 and suppose
liminf Q_h(v_h) is finite. Choose a subsequence realizing that liminf.
For all its sufficiently relevant indices, v_h belongs to V_h.
Weak L2 convergence bounds the L2 norms, and the energy bound gives
[
sup_hig(|v_h|_{L^2}^2+
|
abla v_h|_{L^2}^2ig)<infty.
]
Thus the reconstructed functions are bounded in H1(T_L^d).
After taking a subsequence, they converge weakly in H1 to some w.
The compact embedding H1(T_L^d) into L2(T_L^d) gives strong L2
convergence along a further subsequence. Its limit must equal the
given weak L2 limit v, so w=v and v belongs to H1.

By weak lower semicontinuity of the L2 norm of the gradient,
[
mathcal Q(v)=|
abla v|_{L^2}^2
leliminf_{h	o0}|
abla v_h|_{L^2}^2
=liminf_{h	o0}mathcal Q_h(v_h).
]
If the liminf is infinite, the inequality is immediate. This proves M1.

## 3. M2: strong recovery

First let v be smooth and periodic. Let I_h v be its nodal P1
interpolant. The family is uniformly shape-regular and quasi-uniform:
each element is a scaled copy of one of finitely many reference
Freudenthal simplices. The standard P1 interpolation estimate gives
[
|I_hv-v|_{L^2}	o0,qquad
|
abla I_hv-
abla v|_{L^2}	o0.
]
Consequently,
[
I_hv	o vquad	ext{strongly in }L^2,qquad
mathcal Q_h(I_hv)	omathcal Q(v).
]
For a general v in H1, choose smooth periodic v^(m) converging to v
in H1. Apply the smooth recovery construction for each m and choose
a diagonal sequence m=m(h) tending to infinity sufficiently slowly.
This gives strong L2 convergence and the required upper bound.
For v outside H1, Q(v)=+infinity and the recovery upper bound is
vacuous; any strongly convergent sequence is admissible.

This proves M2.

## 4. Conclusion and assumptions

Under the stated periodic, shape-regular, quasi-uniform Freudenthal
family, exact Euclidean element geometry and common P1 L2 embedding,
M1 and M2 hold. Therefore the P1 Dirichlet forms Mosco-converge to
the continuum torus Dirichlet form.

The proof uses:
- the G4 common-Hilbert-space embedding and exact consistent mass;
- G5 uniform shape regularity and Rellich compactness;
- exact P1 stiffness as the continuum gradient energy;
- standard P1 interpolation estimates and density of smooth periodic
  functions in H1.

The numerical recovery test is a reproducibility diagnostic, not a
replacement for the proof.

## 5. RGS transfer remains open

This theorem is for a flat periodic finite-element family. It does not
establish the curved RGS-Q6 result. The curved transfer still requires
an explicitly defined sequence of coframe-derived meshes, uniform
bi-Lipschitz element maps, metric ellipticity and measure consistency,
compatible common-space embeddings, and liminf/recovery estimates for
the resulting curved forms. A fixed Q6 graph and the flat-torus
theorem do not supply those hypotheses.

Accordingly:
- flat periodic P1 G6: PROVED under the hypotheses above;
- curved RGS G6: OPEN;
- spinorial, SU(2), Wilson/overlap and Yang-Mills Mosco limits: OPEN.
