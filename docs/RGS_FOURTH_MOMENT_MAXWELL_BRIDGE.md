# RGS fourth moment → Maxwell constitutive bridge

## 1. Input already established by the 64-site audit

For the hash-verified RGS coordinate seed (X_n\in\mathbb R^3), (n=1,\ldots,64), define
\[
M^{(4)}_{ijkl}=\sum_n X_{n,i}X_{n,j}X_{n,k}X_{n,l}.
\]
Its Frobenius-orthogonal (SO(3))-isotropic projection is
\[
M^{(4),iso}_{ijkl}=c(\delta_{ij}\delta_{kl}+\delta_{ik}\delta_{jl}
+\delta_{il}\delta_{jk}),\qquad
c=\frac{M^{(4)}_{iikk}}{3(3+2)},
\]
and (R=M^{(4)}-M^{(4),iso}) is fully symmetric, trace-free on every pair contraction, and (O_h)-invariant. The existing audit verifies these finite-seed identities and the (E_g\oplus T_{2g}) splitting.

## 2. The tensorial no-go: no direct (RFF) Maxwell term

The spatial electromagnetic field strength is antisymmetric:
\[
F_{ij}=-F_{ji}.
\]
Since (R_{ijkl}=R_{jikl}), contraction in the first pair gives
\[
R_{ijkl}F_{ij}F_{kl}=0.
\]
This is an exact index-symmetry obstruction, not a numerical cancellation. Thus the tempting direct rank-four correction (R_{ijkl}F_{ij}F_{kl}) is identically zero. The fully symmetric fourth moment cannot be identified with a Riemann-type constitutive tensor or teleparallel torsion.

## 3. Correct intermediate: an STF order parameter

Let (Q_{ij}=Q_{ji}), (Q_{ii}=0), be a separately defined internal/coframe order parameter. The fourth moment induces
\[
A_{ij}(Q)=R_{ijkl}Q_{kl}.
\]
Because (R) is symmetric and trace-free, (A(Q)\in\operatorname{Sym}_0(3)). The (O_h)-equivariance is
\[
A_{g\cdot R}(gQg^T)=gA_R(Q)g^T.
\]
The anisotropy can therefore enter the electric and magnetic response matrices:
\[
Z_E(Q)=z_E I+g_E A(Q),\qquad
Z_B(Q)=z_B I+g_B A(Q),
\]
with (z_E,z_B>0) and positivity conditions
\[
|g_E|\,\|A(Q)\|_2<z_E,\qquad
|g_B|\,\|A(Q)\|_2<z_B.
\]
The 3+1 constitutive density is
\[
\mathcal L_{EM}(Q)=\frac12 E_i(Z_E)_{ij}E_j
-\frac12 B_i(Z_B)_{ij}B_j.
\]
This is gauge invariant because it depends on (F=dA), and it is covariant under simultaneous transformations of the internal frame, (R), (Q), (E), and (B). The implementation checks positive definiteness rather than clipping negative eigenvalues.

## 4. What this does and does not imply for Wilson → Maxwell

The compact Q6 Wilson marginal and its five retained modes are already defined in the repository. To insert the RGS fourth moment, one still needs a specified map
\[
T:\mathbb R^5_{\text{Q6 retained}}\longrightarrow\operatorname{Sym}_0(3),
\qquad Q=Tq,
\]
with the correct group action and equivariance. Equal dimensions (five) do not establish such an intertwiner: the Q6 (A_5) mode representation and the RGS cubic (O_h) representation are not to be identified by dimension alone.

If such a map is derived, (Z_E(Tq)) and (Z_B(Tq)) define a field-dependent constitutive response to be combined with the Wilson effective action. The existing quadratic Gaussian Wilson marginal gives (12\beta_s I_5) and only a (q)-independent determinant; it does not generate this tensor or a Hua factor. The nonlinear compact marginal must supply the measure and dynamics for (q) (or (Q)).

There is also a useful consequence for the vacuum normalization. If (Q) is integrated with a centered, (Q\mapsto-Q) symmetric measure, then (\langle A(Q)\rangle=0): there is no first-order anisotropic correction to the vacuum Maxwell quadratic term. Integrating a Gaussian (Q) coupled linearly to (F^2) produces an (F^4) term at leading nontrivial order, not a shift of the scalar (F^2) coefficient. A nonzero background (Q_0), a non-centered measure, or Wilson-induced field dependence would have to be derived before a quadratic anisotropic response appears.

Consequently this bridge provides a well-defined **candidate constitutive sector**, but it does not yet determine the scalar (Z_F), the renormalized electric charge, or
\[
\alpha=\frac{1}{4\pi Z_F}
\]
(in rationalized conventions). No Wyler normalization or Hua cutoff is inserted.

## 5. Evidence status

- **Exact:** the direct (RFF) contraction vanishes by index symmetry; (R:Q) is STF; the positivity bounds above are sufficient and sharp for the chosen (Q); simultaneous-frame covariance follows from tensor transformation.
- **Numerically audited:** the actual 64-site seed, its fourth-moment decomposition, (O_h) invariance, and the induced constitutive matrices for the documented test (Q).
- **Candidate:** the (Q)-dependent electric/magnetic constitutive matrices as an RGS material response.
- **Open:** derive (Q) from the coframe/PSU–Anti-Twin geometry; construct an equivariant Q6 (A_5\to\operatorname{Sym}_0(3)) map if that is the intended identification; derive the (Q)-measure from the complete nonlinear Wilson fibre; prove the continuum limit; determine the scalar Maxwell normalization and only then compare with Wyler.

## 6. Reproduction and references

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_rgs_fourth_moment_maxwell_bridge.py -q
PYTHONPATH=. python -m src.rgs_fourth_moment_maxwell_bridge
\`\`\`

The rank-four seed construction is documented in [RGS_64_FOURTH_MOMENT_BRIDGE.md](RGS_64_FOURTH_MOMENT_BRIDGE.md). The Wilson lattice-gauge and transfer-matrix framework is Kogut, *Rev. Mod. Phys.* 51 (1979), 659–713, DOI: 10.1103/RevModPhys.51.659. The general constitutive-tensor formulation is treated in Hehl & Obukhov, *Foundations of Classical Electrodynamics* (Birkhäuser, 2003), DOI: 10.1007/978-94-010-0205-8. These references support the general frameworks; they do not prove the RGS-specific coupling or its physical interpretation.
