# Wyler, the fine-structure constant, and the Q6/RGS programme

## 1. Physical target

The target is the low-energy (Thomson-limit) electromagnetic coupling

\[
\alpha(0)=\frac{e^2}{4\pi\varepsilon_0\hbar c},
\qquad
\alpha(0)^{-1}=137.035999177(21)
\]

using the 2022 CODATA adjustment. This is a dimensionless observable; at other momentum scales the effective QED coupling runs. A geometric derivation must therefore state which renormalization/physical scale it predicts.

## 2. Wyler's expression

The expression commonly called Wyler's constant is

\[
\alpha_W=
\frac{9}{8\pi^4}
\left(\frac{\pi^5}{2^4\,5!}\right)^{1/4},
\qquad
\alpha_W^{-1}=137.0360824481643\ldots
\]

It was motivated by ratios of volumes associated with bounded symmetric domains and geometries related to the conformal symmetry of relativistic wave equations. Its numerical proximity to the low-energy fine-structure constant is striking, but the expression is not an established derivation of QED's coupling.

With CODATA 2022, the inverse-value difference is

\[
\alpha_W^{-1}-\alpha(0)^{-1}
=8.3271164\times10^{-5},
\]

or about 0.608 ppm relative to the CODATA inverse value. This is a numerical comparison, not an uncertainty estimate for Wyler's formula.

## 3. The decisive mathematical objection

Robertson's 1971 comment in *Physical Review Letters* points out that the agreement requires a radius normalization of the bounded spaces to be set to one, while the symmetry-group relation itself does not fix that radius. No physical argument for that choice was supplied in Wyler's construction.

Therefore a valid RGS reconstruction cannot merely reproduce the same ratio of volumes. It must derive, from independently specified geometric data:

1. the relevant homogeneous/bounded domain and its invariant measure;
2. the metric and its normalization, including why no arbitrary radius remains;
3. the map from the geometric measure to the coefficient of the electromagnetic kinetic term;
4. the normalization of the gauge field and electric charge;
5. the physical scale at which the resulting coupling equals \(\alpha(0)\).

If a free scale or fitted prefactor is introduced to force the numerical match, the derivation is not predictive.

## 4. Theorem W1 — normalization obstruction

**Statement.** Suppose a proposed geometric expression for a dimensionless coupling is assembled from volumes of homogeneous spaces, but the geometric construction leaves an unconstrained metric scale (or equivalent normalization parameter) that changes the proposed volume combination. Then the symmetry data alone do not determine a unique numerical coupling.

**Proof.** A homogeneous space may admit a family of invariant metrics related by an overall scale. The symmetry action identifies the homogeneous structure but, absent an additional normalization condition, does not select one member of that family. If the proposed volume functional is not invariant under this rescaling, its numerical value varies with the undetermined scale. Hence the symmetry data cannot determine a unique number. An independent physical or geometric normalization is necessary. \(\square\)

This is the precise obstruction that an RGS derivation must overcome. It does not prove that Wyler's number is accidental; it proves that the numerical expression alone does not close the derivation.

## 5. Theorem W2 — what a genuine geometric derivation must establish

**Statement (criterion).** Let a discrete or continuum geometric model produce a quadratic electromagnetic action
\[
S_{\rm em}[A]=\frac{1}{4g^2_{\rm geom}}\int F_{\mu\nu}F^{\mu\nu}\,d\mu_g
\]
with all field and measure conventions fixed. If (i) the measure normalization is fixed independently of the target value of \(\alpha\), (ii) the gauge field is canonically normalized relative to the matter charge representation, and (iii) the continuum/renormalized limit exists at a specified scale, then the model predicts a definite coupling
\[
\alpha_{\rm geom}(\mu)=\frac{g^2_{\rm geom}(\mu)}{4\pi}
\]
in rationalized natural units.

**Proof.** Under the stated hypotheses, the coefficient of the gauge kinetic term and the normalization of the charge representation are fixed, so the dimensionless ratio defining the coupling is no longer adjustable by a field rescaling. Existence of the specified limit gives a definite value at the stated scale. \(\square\)

The theorem is a conditional criterion, not a claim that the current Q6 construction already satisfies the hypotheses.

## 6. RGS–Q6 research route

The useful bridge is not to insert Wyler's formula into Q6 by hand. Instead:

- **Geometry:** identify the actual homogeneous space arising from the RGS coframe and symmetry action. Do not assume that a six-cube automatically yields Wyler's bounded domain.
- **Measure:** construct its invariant volume/Haar normalization and prove scale independence or derive the missing scale from an independent RGS condition.
- **Gauge sector:** derive the U(1) electromagnetic kinetic term as a specified sector/reduction of the gauge connection. The existing Q6 SU(2) links and Dirac–Kähler operator do not by themselves produce electromagnetism.
- **Matter normalization:** specify the charge representation and its normalization; only then can the coefficient be interpreted as electric charge.
- **Continuum and running:** prove convergence of the discrete action/operator and match the coupling at a stated scale. The low-energy CODATA value is not interchangeable with the coupling at the electroweak scale.
- **Independent predictions:** after fixing the geometry without using \(\alpha\), calculate at least one further observable or scale dependence. This guards against a one-number fit.

## 7. Evidence ledger

**Established external facts**
- Wyler's expression and its numerical value are as stated above.
- Robertson identified the unresolved radius normalization.
- CODATA 2022 gives \(\alpha^{-1}=137.035999177(21)\).

**Exact mathematics in this note**
- The normalization-obstruction theorem W1.
- The conditional action-normalization criterion W2.

**Not established**
- That the RGS/Q6 geometry selects Wyler's bounded domain.
- That its invariant measure yields the electromagnetic kinetic coefficient.
- That \(\alpha\) has been derived from Q6, PSU, 2I, or the coframe.
- That Wyler's numerical agreement is other than a striking numerical coincidence.

## References

- A. Wyler, “L'espace symétrique du groupe des équations de Maxwell,” *Comptes Rendus de l'Académie des Sciences de Paris* **269**, A743–A745 (1969).
- B. Robertson, “Wyler's Expression for the Fine-Structure Constant \(\alpha\),” *Physical Review Letters* **27**, 1545 (1971), DOI: 10.1103/PhysRevLett.27.1545.
- S. L. Adler, “Theories of the Fine Structure Constant,” Fermilab report 72/059-T (1972).
- NIST/CODATA 2022, recommended value of the fine-structure constant, \(\alpha^{-1}=137.035999177(21)\).
