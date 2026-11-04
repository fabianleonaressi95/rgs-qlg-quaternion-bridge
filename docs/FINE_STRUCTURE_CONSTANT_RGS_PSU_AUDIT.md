# Fine-structure constant: RGS–PSU–Q6 audit and falsifiable derivation programme

**Date:** 4 November 2026  
**Reference:** CODATA 2022 (latest recommended adjustment currently published by NIST; the 2026 adjustment is not yet complete).

## 1. Physical target and convention

The fine-structure constant is the dimensionless electromagnetic coupling. In the low-energy (zero-momentum/Thomson) convention,

\[
\alpha=\frac{e^2}{4\pi\epsilon_0\hbar c},\qquad
\alpha^{-1}=137.035999177(21).
\]

The CODATA value and uncertainty are taken from the NIST 2022 adjustment. The 2026 adjustment is scheduled to use data published by 31 December 2026; it is not a replacement recommended value yet.

An RGS derivation must predict this dimensionless number from independently fixed mathematical/physical inputs. It must not insert \(e\), \(\alpha\), the measured value of \(\alpha\), or a fitted normalization into the geometric construction and then claim a prediction.

Also distinguish:
- the low-energy fine-structure constant above;
- the running QED coupling \(\alpha(q^2)\), which depends on momentum/energy scale;
- an effective lattice coupling \(g^2/(4\pi)\) or a Wilson-action parameter, whose relation to continuum QED must be derived by matching and renormalization.

## 2. Existing RGS/PSU numerical candidates recovered from prior work

### 2.1 Candidate involving the Feigenbaum constant

A previously explored expression is

\[
\boxed{\alpha_{\rm cand}^{-1}
=\frac{2\pi^2(\pi^\pi-1)F_a^2}{32}},
\qquad F_a\simeq2.5029078750958928,
\]

where \(F_a\) is the Feigenbaum alpha constant (the period-doubling universal constant).

Recalculation gives

\[
\alpha_{\rm cand}^{-1}=137.03599969886864.
\]

Against the CODATA reference:

- absolute difference: \(5.2186863\times10^{-7}\);
- relative difference: \(3.8083\times10^{-9}\), about 3.81 ppb;
- using the CODATA standard uncertainty \(2.1\times10^{-8}\), the discrepancy is about 24.9 reference standard uncertainties, before assigning any theoretical uncertainty to the candidate.

This is a striking numerical proximity, but it is **not** a derivation and it does not agree within the CODATA uncertainty. The factor \(2\pi^2=\mathrm{Vol}(S^3)\) has a clear geometric meaning. The appearance of \(32\), \(\pi^\pi-1\), and especially the Feigenbaum constant squared as the electromagnetic normalization has not been derived from an RGS/PSU action, a group representation, a spectral determinant, or a renormalization calculation. Choosing the expression after looking at the target value creates a multiple-testing/selection-bias problem.

**Status: NUMERICAL CANDIDATE ONLY; NOT A PREDICTION.**

### 2.2 Wyler expression

The historically discussed expression

\[
\alpha_W^{-1}
=2^{19/4}3^{-7/4}5^{1/4}\pi^{11/4}
\simeq137.0360824481643
\]

is not derived by the present RGS–PSU construction. Its discrepancy from the CODATA reference is about \(8.3271\times10^{-5}\) (roughly 0.607 ppm). It is retained only as a comparison baseline, not as evidence for the theory.

### 2.3 Earlier script in the main repository

The existing file 'experiments/analysis_sequence/scriptanal27.py' defines

- \(W_{geo}=0.188562\),
- \(\hbar_{geo}=2.828427\),
- \(n_{edges}=12\),
- \(n_{faces}=24\),
- \(v_{fraction}=2/3\).

It then computes

\[
\alpha^{-1}_{script}
=\frac{16\pi^2\hbar_{geo}}{W_{geo}n_{edges}}
\simeq197.391878656,
\]

and a second expression

\[
\alpha^{-1}_{Schmidt}
=\frac{4\pi^3}{(1/4+1/12)W_{geo}}
\simeq1973.225359105.
\]

Neither gives the target \(137.036\). The code's comments call these attempts a derivation, but the numerical values do not support that label. In addition, \(W_{geo}\) is a hard-coded decimal whose independent definition, units, uncertainty and derivation are not encoded in this script. The variable 'term_base' is computed but not used in the first expression. This script must be labelled exploratory and must not be treated as a validated alpha engine.

### 2.4 Golden-angle comparison

The golden angle \(360^\circ/\varphi^2\simeq137.508^\circ\) is an angle, whereas \(\alpha^{-1}\) is a dimensionless coupling. A numerical resemblance between the angle in degrees and 137 is not an equation: changing angular units changes the number. A physical relation would need a coordinate-independent construction, a specified normalization and a derivation from an action or observable.

## 3. What a genuine RGS–PSU derivation must supply

The electromagnetic coupling is not fixed by the topology of a graph alone. In a continuum gauge theory, it is the coefficient/normalization of the gauge kinetic term after the gauge field and charge conventions have been fixed. A discrete-to-continuum derivation therefore needs all of the following.

### Gate A — one declared electromagnetic action

Start with a convention-explicit action, for example

\[
S_{EM}=-\frac{1}{4g_0^2}\int d^4x\sqrt{-g}\,
F_{\mu\nu}F^{\mu\nu}
+S_{matter}[\psi,A_\mu],
\]

or its SI-normalized equivalent. State whether the gauge potential is normalized so that matter has unit charge or physical charge \(e\). Avoid counting the same coupling both in \(F=dA\) and in the covariant derivative.

For the RGS–PSU construction, specify the discrete link variables, gauge group, representation, plaquette action and continuum field-strength limit. Internal SU(2)/SU(3) links must not be confused with the spacetime Lorentz/Weitzenböck connection.

### Gate B — exact geometric and group-theoretic input

For any use of Q6, PRISM, \(2I\), the 600-cell, \(S^3/2I\), \(\varphi\), or the supergolden constant \(\psi\), define the exact object and map:
- group order, multiplication, inverses and representations;
- graph/complex incidence operators and orientation;
- invariant measure and normalization;
- the map/intertwiner connecting the discrete operator to the continuum gauge sector.

In particular, do not use the KILAW 'naressilaw' \(2I\) data as a source until its 216-element output is replaced by an exact 120-element binary icosahedral group and its group axioms are tested.

### Gate C — derive, do not fit, the gauge kinetic coefficient

Derive the coefficient of \(F_{\mu\nu}F^{\mu\nu}\) from the declared discrete action or a controlled continuum limit. Show that the normalization is independent of arbitrary rescalings of links, generators, trace conventions, cell volumes and lattice spacing. If the result depends on a regulator, demonstrate its removal.

### Gate D — renormalization and physical matching

Compute or justify the relation between the bare lattice coupling and the low-energy Thomson coupling. This requires a renormalization prescription, beta function/threshold treatment, matter content and a controlled continuum limit. A bare plaquette coefficient is not automatically \(\alpha^{-1}\).

### Gate E — parameter-free prediction

Freeze all geometric constants and conventions before comparing with CODATA. No tuning of \(W_{geo}\), powers of \(\pi\), integer factors, \(\varphi/\psi\) exponents, or spectral windows after observing the target. Include uncertainty propagation and a declared family of candidate expressions to control selection bias.

## 4. Required numerical audit and tests

Implement the following as a reproducible module, separate from production physics until the gates pass:

1. **Reference constants:** store the CODATA reference, uncertainty, source, adjustment year and convention in a machine-readable record.
2. **Candidate registry:** record every previously proposed formula, its exact inputs, derivation status, numerical output, absolute/relative error and whether any parameter was fitted.
3. **Dimensional/invariance tests:** prove the candidate is dimensionless and invariant under unit changes and allowed field/generator rescalings.
4. **Sensitivity:** vary every non-exact input within its justified uncertainty; report propagated uncertainty.
5. **Selection-bias control:** compare against a predeclared set of alternatives, not just the closest formula selected after the fact.
6. **Action-to-observable bridge:** compute the effective Maxwell term and match its coefficient to the Thomson coupling.
7. **Independent reproduction:** a second implementation should reproduce the exact numerical candidate from only the frozen input record.

Minimum tests:
- CODATA reference value and uncertainty are stored exactly;
- candidate formula is recomputed rather than hard-coded;
- no reference alpha is imported into candidate construction;
- no hidden fitted coefficient is allowed;
- changing units does not change the dimensionless output;
- the low-energy matching convention is explicit;
- candidate status cannot be upgraded to a prediction solely by a small numerical residual.

## 5. Integration points with current RGS work

| Existing subsystem | Legitimate role in alpha programme | Missing bridge |
|---|---|---|
| Clifford/Dirac and EC–Dirac | Matter representation and current coupling | Consistent U(1) gauge covariant derivative and action normalization |
| Q6 cochains | Discrete differential/field-strength construction | A continuum Maxwell limit and coupling matching |
| PRISM/Krist-All | Candidate spectral/geometric substrate | Proven map from spectral data to the Maxwell kinetic coefficient |
| \(2I\), 600-cell, \(S^3/2I\) | Exact symmetry/geometry input | Validated group action and gauge-sector intertwiner |
| Vierbein/LC/TEGR | Spacetime geometry and measure | Demonstrate how geometry affects, or does not affect, gauge normalization |
| Q-ball/QED/MHD | Charged matter and electromagnetic response | One common action, current normalization and conserved charge |
| PSU/Anti-Twin | Candidate two-sheet geometry | Derive how the sheets determine the U(1) coupling without fitted weights |

## 6. Evidence ledger

| Claim | Status | Reason |
|---|---|---|
| CODATA \(\alpha^{-1}=137.035999177(21)\) | REFERENCE VALUE | NIST/CODATA 2022 |
| Feigenbaum/\(\pi\) candidate gives 137.0359996989 | NUMERICALLY VERIFIED ARITHMETIC | Recomputed; mismatch is about 24.9 CODATA standard uncertainties |
| Feigenbaum/\(\pi\) formula derives electromagnetic coupling | HYPOTHESIS / OPEN | No action, group or renormalization derivation |
| Wyler expression | HISTORICAL NUMERICAL COMPARATOR | Not derived in RGS–PSU |
| 'scriptanal27.py' first/second outputs | NUMERICALLY VERIFIED ARITHMETIC; NOT A DERIVATION | 197.3919 and 1973.2254, respectively |
| Golden angle equals inverse fine-structure constant | REJECTED AS STATED | Degrees are unit-dependent; no invariant physical equation |
| RGS–PSU predicts low-energy alpha | OPEN | Gates A–E not passed |

## 7. Sources

- NIST, *CODATA Recommended Values of the Fundamental Physical Constants: 2022*, published 30 April 2025: https://www.nist.gov/publications/codata-recommended-values-fundamental-physical-constants-2022
- NIST CODATA database and version history: https://physics.nist.gov/cuu/Reference/versioncon.html
- Main repository exploratory script: 'experiments/analysis_sequence/scriptanal27.py'

**Conclusion:** alpha is now an explicit research target in the RGS–PSU–Q6 integration programme. The near-match involving \(F_a\) is worth retaining as a candidate to explain or falsify, but its numerical closeness is not a physical derivation. The immediate task is to derive the Maxwell kinetic normalization from the same declared discrete/continuum action used for the gauge, matter and vierbein sectors, then match it to the low-energy coupling without fitting.
