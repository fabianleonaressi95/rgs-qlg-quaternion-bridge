# Wyler D5: invariant measure, Euclidean volume, and the scale obstruction

## 1. The domain and the measure must not be conflated

The Wyler expression uses the Euclidean volume factor

\[
V_E(D_5;1)=\frac{\pi^5}{2^4 5!}=\frac{\pi^5}{1920},
\qquad
\alpha_W(1)=\frac{9}{8\pi^4}V_E(D_5;1)^{1/4}.
\]

Here (D_5) is the five-complex-dimensional bounded symmetric domain in the normalization used in the Wyler expression. The real dimension is ten.

There are two distinct measures:

- (dV_E): Euclidean Lebesgue measure in a chosen bounded realization. Its total volume is finite, but it is not invariant under the full noncompact automorphism group.
- (dV_B): the volume measure of the invariant Bergman metric. This is intrinsic up to the normalization convention for the metric, but the complete invariant metric on a noncompact bounded symmetric domain has infinite total volume.

Thus one cannot silently call the finite Euclidean volume a fully invariant homogeneous-space volume. The measure and its normalization are part of the mathematical input.

## 2. Theorem W3 — exact scale law

**Statement.** Let the Euclidean realization be dilated by (z\mapsto Rz), (R>0), and let (D_5(R)=R D_5(1)). Then

\[
V_E(D_5;R)=R^{10}V_E(D_5;1),
\qquad
\alpha_W(R)=R^{5/2}\alpha_W(1).
\]

**Proof.** Lebesgue measure in real dimension ten transforms under dilation by the Jacobian determinant (R^{10}). The Wyler expression contains the fourth root of this volume, so its scale factor is ((R^{10})^{1/4}=R^{5/2}). \(\square\)

This proves that the finite Euclidean volume expression is not independent of the radius normalization. It is the explicit scaling obstruction behind the issue raised by Robertson.

## 3. Theorem W4 — invariant-volume obstruction

**Statement.** A complete, connected, noncompact homogeneous Riemannian manifold with its invariant Riemannian volume cannot have finite total volume.

**Proof sketch.** Homogeneity gives a uniform positive lower bound on the volume of sufficiently small metric balls. Completeness and noncompactness imply the existence of an infinite sequence of pairwise disjoint balls of a fixed sufficiently small radius (choose a sequence of points escaping every compact set). Their volumes have the same positive lower bound, so the total volume is infinite. \(\square\)

The Bergman metric on an irreducible bounded symmetric domain is complete and invariant under its automorphism group. Its total Bergman volume is therefore infinite. The finite (V_E(D_5)) in Wyler's formula is not that invariant Bergman volume.

## 4. Numerical audit and the fitted-radius diagnostic

The code in `src/rgs_wyler_domain_audit.py` computes the formula and the exact scale law. With the unit-radius Euclidean normalization it obtains

\[
\alpha_W^{-1}=137.0360824481643\ldots
\]

Using CODATA 2022, (\alpha(0)^{-1}=137.035999177), the inverse-coupling difference is about (0.608) ppm. Solving the scale law for a radius that forces agreement gives

\[
R_{\rm fit}=\left(\frac{\alpha_{\rm CODATA}}{\alpha_W(1)}\right)^{2/5}
\approx 1.0000002431.
\]

This is a diagnostic only. Since (R_{\rm fit}) is obtained from the measured target, it is not an independent prediction and cannot resolve the normalization objection.

Reproduce:

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_wyler_domain_audit.py -q
PYTHONPATH=. python -c "from src.rgs_wyler_domain_audit import audit_wyler_scaling; print(audit_wyler_scaling())"
```

## 5. What RGS would have to prove

A viable RGS route cannot simply declare (R=1). It must construct an independent map from the RGS geometric data to the ten-real-dimensional tangent metric of the candidate (D_5) domain and prove that the induced Euclidean realization has a fixed scale. In particular:

1. Identify the actual (SO(5,2)/(SO(5)\times SO(2))) homogeneous geometry (or justify a different domain) from the RGS variables.
2. Construct its invariant metric from a specified normalization, for example a Killing/Bergman convention, and distinguish that metric from the bounded-coordinate Euclidean metric.
3. Derive the coordinate Jacobian connecting the RGS metric volume form to (dV_E). This is the missing scale map; it must not be fitted to (alpha).
4. Construct the electromagnetic (U(1)) kinetic term and matter charge normalization. A Q6 (SU(2)) link field or Dirac–Kähler operator alone does not supply this.
5. Prove the appropriate continuum and renormalization limit, then state the scale at which the predicted coupling is compared with (\alpha(0)).

The existing six-dimensional Q6 cube has a 64-dimensional exterior algebra and 729 cubical cells; neither fact by itself identifies a ten-real-dimensional Cartan domain or fixes its Euclidean radius. Any such map is a new theorem to be constructed, not an identification by dimension-counting.

## 6. Evidence ledger

- **Exact:** the dilation law (V_E\sim R^{10}), (\alpha_W\sim R^{5/2}).
- **Numerically audited:** Wyler's unit-radius value and the radius required to fit CODATA 2022.
- **Geometric fact:** finite Euclidean coordinate volume and invariant Bergman volume are different measures; the latter is infinite on this complete noncompact domain.
- **Open:** derivation of the domain, scale map, U(1) action, charge normalization, and QED matching from RGS/PSU/Q6.
- **Not claimed:** derivation of the measured fine-structure constant.

## References

- A. Wyler, “L'espace symétrique du groupe des équations de Maxwell,” *C. R. Acad. Sci. Paris* **269**, A743–A745 (1969).
- B. Robertson, “Wyler's Expression for the Fine-Structure Constant (\alpha),” *Phys. Rev. Lett.* **27**, 1545–1547 (1971), DOI: 10.1103/PhysRevLett.27.1545.
- S. L. Adler, “Theories of the Fine Structure Constant,” Fermilab report 72/059-T (1972).
- NIST/CODATA 2022 recommended value of (\alpha).
