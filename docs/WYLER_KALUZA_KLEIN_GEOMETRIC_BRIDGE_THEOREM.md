# Theorem: Kaluza–Klein circle reduction and the Wyler geometric proposal

**Status:** exact kinematic Kaluza–Klein results; conditional dimensional-reduction statement; no established Kaluza–Klein derivation of Wyler's fine-structure expression.

This note formalizes the useful mathematical relation between a compact fifth dimension and a four-dimensional \(U(1)\) gauge field, then states precisely what is—and is not—connected to Wyler's bounded-domain expression. It supplements, and does not replace:

- \`docs/WYLER_D5_INVARIANT_MEASURE_AND_SCALE_THEOREM.md\`;
- \`docs/WYLER_FINE_STRUCTURE_CONSTANT_RGS_THEOREM.md\`;
- \`docs/RGS_GODEL_WYLER_ONTOLOGY_BRIDGE.md\`.

## 1. Geometric setup

Let \(M^4\) be a four-dimensional spacetime and let the fifth coordinate be an angular variable
\[
\theta\in\mathbb R/(2\pi\mathbb Z).
\]
For a constant circle radius \(R>0\), consider the Kaluza–Klein metric ansatz
\[
ds_5^2
=g_{\mu\nu}(x)\,dx^\mu dx^\nu
+R^2\bigl(d\theta+q A_\mu(x)\,dx^\mu\bigr)^2,
\tag{1}
\]
where \(q\) is a convention-dependent normalization constant. The ansatz is a geometric decomposition, not yet a complete dynamical model. A general reduction also contains a scalar radion if the circle radius is allowed to vary.

## 2. Theorem KK1 — the circle isometry induces a \(U(1)\) gauge transformation

**Statement.** The metric (1) is invariant under the local change of fibre coordinate
\[
\theta'=\theta-q\lambda(x),\qquad
A'=A+d\lambda.
\tag{2}
\]
Consequently, the one-form \(A\) transforms as an Abelian gauge potential, and its curvature
\[
F=dA
\tag{3}
\]
is invariant.

**Proof.** Direct substitution gives
\[
d\theta'+qA'
=d\theta-q\,d\lambda+qA+q\,d\lambda
=d\theta+qA.
\]
Thus the one-form appearing in the metric is unchanged. Moreover,
\(F'=d(A+d\lambda)=dA=F\), since \(d^2\lambda=0\). \(\square\)

**Interpretation.** The four-dimensional \(U(1)\) gauge redundancy is the local description of the circle-fibre coordinate freedom in this ansatz. This is a precise geometric statement; it does not by itself fix the physical electric charge or the measured electromagnetic coupling.

## 3. Theorem KK2 — compactness quantizes the fifth momentum

**Statement.** Let a field be single-valued on the circle, so that
\[
\Psi(x,\theta+2\pi)=\Psi(x,\theta).
\]
Then it has a Fourier expansion
\[
\Psi(x,\theta)=\sum_{n\in\mathbb Z}\psi_n(x)e^{in\theta},
\]
and the eigenvalues of the physical fifth-direction momentum are
\[
p_5=\frac{n\hbar}{R},\qquad n\in\mathbb Z.
\tag{4}
\]
For a free higher-dimensional field, the corresponding four-dimensional Kaluza–Klein modes acquire a mass contribution proportional to \(n^2/R^2\) (with the appropriate factors of \(\hbar\) and \(c\)).

**Proof.** Periodicity requires the angular wave number to be an integer \(n\). The physical arc-length coordinate is \(y=R\theta\), so the momentum operator \(-i\hbar\partial_y\) acting on \(e^{in\theta}\) has eigenvalue \(n\hbar/R\). Squaring this momentum gives the stated \(R^{-2}\) dependence. \(\square\)

**Boundary.** Equation (4) fixes the mode spacing once \(R\) is specified. It does not determine \(R\).

## 4. Theorem KK3 — dimensional reduction produces a gauge kinetic term, conditionally

**Statement.** Start from a specified five-dimensional Einstein–Hilbert action and use the circle-invariant metric ansatz (1), with a fixed radion and fixed conventions. Integrating over the compact circle yields a four-dimensional action containing a Maxwell-type term
\[
S_4\supset
-\frac{1}{4g_4^2}\int_{M^4}F_{\mu\nu}F^{\mu\nu}\sqrt{|g|}\,d^4x.
\tag{5}
\]
The coefficient \(g_4\) is determined by the higher-dimensional gravitational normalization, the circle radius, and the normalization \(q\) used in (1).

**Justification.** The five-dimensional Ricci scalar for a circle-invariant metric contains the four-dimensional Ricci scalar and a term quadratic in the curvature \(F=dA\), together with radion terms when the radius is dynamical. The five-dimensional volume element contributes the circle factor. Integrating the circle coordinate therefore gives (5), with a coefficient fixed only after the five-dimensional action, radius, field convention, and radion treatment have all been specified. \(\square\)

This is a conditional reduction statement, not a numerical prediction for \(g_4\). In particular, a rescaling of \(A\) changes the displayed coefficient and simultaneously changes the charge convention; the matter representation must be normalized consistently before the coefficient is identified with the physical electric coupling.

## 5. Wyler's expression and its independent scale issue

The expression commonly attributed to Wyler is
\[
\alpha_W=
\frac{9}{8\pi^4}
\left(\frac{\pi^5}{2^4\,5!}\right)^{1/4},
\qquad
\alpha_W^{-1}=137.0360824481643\ldots.
\tag{6}
\]
The finite Euclidean volume factor used in this expression is
\[
V_E(D_5;1)=\frac{\pi^5}{1920}
\]
in the unit-radius normalization of the existing repository theorem.

For a ten-real-dimensional Euclidean realization dilated by \(z\mapsto R_Dz\),
\[
V_E(D_5;R_D)=R_D^{10}V_E(D_5;1),
\qquad
\alpha_W(R_D)=R_D^{5/2}\alpha_W(1).
\tag{7}
\]
Thus the Euclidean-volume expression retains a radius dependence. The finite Euclidean coordinate volume must not be conflated with the total invariant Bergman volume, which is infinite on the complete noncompact domain. The scale obstruction and its proof are recorded in \`WYLER_D5_INVARIANT_MEASURE_AND_SCALE_THEOREM.md\`.

## 6. Theorem KW1 — what follows, and what does not follow, from the two constructions

**Statement.** The Kaluza–Klein circle ansatz establishes a \(U(1)\) gauge redundancy and, for a compact circle of specified radius, a quantized fifth-momentum spectrum. These facts alone do **not** determine Wyler's bounded domain, its Euclidean radius normalization, or the fine-structure constant.

**Proof.** KK1 and KK2 depend on the circle bundle, its metric ansatz, and the radius \(R\). Wyler's expression depends on a separate bounded-domain volume construction and, by (7), on its Euclidean scale \(R_D\). No identification \(R_D=R\), no canonical map from the Kaluza–Klein quotient to the Wyler domain, and no relation between the reduced gauge kinetic coefficient and the particular volume combination in (6) is supplied by KK1 or KK2. Therefore the stated Kaluza–Klein results do not logically imply (6) or select its normalization. \(\square\)

This is a non-implication result about the currently specified inputs. It does not prove that a deeper relation is impossible; it specifies the missing theorem.

## 7. Explicit bridge obligations

A valid Kaluza–Klein–Wyler derivation would have to prove all of the following without fitting to the measured \(\alpha\):

1. **Domain map:** construct a canonical map from the reduced geometric data to the particular bounded symmetric domain used in the Wyler expression. Identify the acting group, isotropy subgroup, and global domain.
2. **Scale relation:** derive a relation between the circle radius \(R\), the domain scale \(R_D\), and the metric normalization. It must follow from independent geometric or physical data.
3. **Measure identity:** prove the exact Jacobian relating the reduced measure to the Euclidean measure used in (6), while keeping it distinct from invariant Bergman measure.
4. **Gauge normalization:** derive the \(F^2\) coefficient from a specified higher-dimensional action and fix the normalization of the matter charge representation.
5. **Radion and spectrum:** account for the radion dynamics, the Kaluza–Klein tower, and any decoupling or truncation approximation.
6. **Four-dimensional limit:** establish the effective action and its regime of validity, including the renormalization scale at which it is compared with \(\alpha(0)\).
7. **Independent test:** predict at least one additional observable or scale dependence after all geometric normalizations have been fixed.

A numerical equality obtained by choosing \(R\), \(R_D\), \(q\), or a prefactor to match \(\alpha\) is a fit, not a derivation.

## 8. RGS–PSU–PRISM integration plan

The existing repository already contains exact finite-dimensional Q6 constructions and separate D5/Wyler scale audits. The present theorem adds no new identification between them. A disciplined next stage is:

- use the RGS vierbein and connection to define a candidate higher-dimensional metric, if such a fifth direction is actually present in the model;
- specify whether the circle is a genuine fibre, a quotient variable, or only an analogy (these are mathematically different);
- derive the \(U(1)\) action and its charge normalization from the chosen action;
- test whether a well-defined quotient or continuum limit produces the relevant bounded domain;
- prove the scale and measure identities before evaluating any candidate coupling.

The existing Q6 \(SU(2)\) link sector is not, by itself, the electromagnetic \(U(1)\) sector. Likewise, the PRISM finite spectrum does not by itself establish a compact fifth dimension or a Kaluza–Klein tower.

## 9. Evidence ledger

| Claim | Status |
|---|---|
| Circle-fibre coordinate freedom induces \(U(1)\) gauge transformations | Exact from the metric ansatz |
| Periodicity gives \(p_5=n\hbar/R\) | Exact for single-valued circle modes |
| Einstein–Hilbert reduction contains an \(F^2\) term | Standard conditional dimensional-reduction result |
| The four-dimensional coupling is numerically fixed in the current RGS model | Not established |
| Kaluza–Klein selects Wyler's bounded domain | Not established |
| The KK radius fixes Wyler's Euclidean radius | Not established |
| Wyler's expression derives the measured fine-structure constant | Not established |
| A Kaluza–Klein–Wyler bridge is mathematically impossible | Not claimed |

## References

- T. Kaluza, “Zum Unitätsproblem der Physik,” *Sitzungsberichte der Preussischen Akademie der Wissenschaften* (1921), 966–972. A modern translation/republication is available as arXiv:1803.08616.
- O. Klein, “The Atomicity of Electricity as a Quantum Theory Law,” *Nature* **118** (1926), 516.
- O. Klein, “Quantentheorie und fünfdimensionale Relativitätstheorie,” *Zeitschrift für Physik* **37** (1926), 895–906.
- A. Wyler, “L'espace symétrique du groupe des équations de Maxwell,” *Comptes Rendus de l'Académie des Sciences de Paris* **269A** (1969), 743–745.
- B. Robertson, “Wyler's Expression for the Fine-Structure Constant \(\alpha\),” *Physical Review Letters* **27** (1971), 1545–1547, https://doi.org/10.1103/PhysRevLett.27.1545.
- R. Gilmore, “Scaling of Wyler's Expression for \(\alpha\),” *Physical Review Letters* **28** (1972), 462.
- S. L. Adler, “Theories of the Fine Structure Constant,” Fermilab report 72/059-T (1972).
