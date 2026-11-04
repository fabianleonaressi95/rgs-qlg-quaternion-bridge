# Q6–D5 Bergman/Liouville pullback: exact local measure theorem

**Status:** exact measure identity on the D5 preimage of the five-mode Q6 tangent phase space. This is not a global identification of the full Q6 quotient with D5.

## 1. Maxwell-normalized Q6 phase-space sector

The Q6 coexact Hodge eigenspace at eigenvalue \(\lambda=12\) has dimension five. Its cotangent lift has canonical coordinates \((q,p)\in\mathbb R^5\times\mathbb R^5\). In the Maxwell convention already used by the repository, the quadratic Hamiltonian is

\[
H_5(q,p)=\frac12\left(\beta\lambda |q|^2+\beta^{-1}|p|^2\right),
\qquad \omega=\sqrt{\lambda}.
\]

The frequency is independent of \(\beta\). Define

\[
z=\Phi_\beta(q,p)
=\frac{1}{\sqrt{20}}
\left(\sqrt{\beta\omega}\,q
+i\frac{p}{\sqrt{\beta\omega}}\right).
\]

Then

\[
H_5=10\omega |z|^2,
\qquad
\left|\det\frac{\partial(\Re z,\Im z)}{\partial(q,p)}\right|
=20^{-5}.
\]

Thus \(d\mu_L=d^5q\,d^5p=20^5d^{10}z\). The Jacobian is independent of the bare coefficient \(\beta\). The chart is a linear isomorphism of the ten-dimensional tangent phase space with \(\mathbb C^5\); its restriction to \(\Phi_\beta^{-1}(D^{IV}_5)\) is the domain on which the Bergman pullback below is defined.

## 2. Type-IV generic norm and Bergman metric

For the bounded type-IV domain,

\[
N(z)=1-2z^\dagger z+|z^Tz|^2,
\qquad
D^{IV}_5=\{z:N(z)>0,\ |z^Tz|<1\}.
\]

With the repository's Bergman-kernel convention \(K(z,z)=C N(z)^{-5}\), the Hermitian metric is

\[
g_{i\bar j}(z)=\partial_i\partial_{\bar j}\log K(z,z)
=-5\partial_i\partial_{\bar j}\log N(z).
\]

At the origin, \(g(0)=10I_5\).

The determinant identity for the type-IV Bergman metric in this normalization is

\[
\boxed{\det g(z)=\frac{10^5}{N(z)^5}}.
\]

It can be obtained by differentiating \(N\) and applying the matrix determinant lemma to the resulting identity-plus-rank-two update; the normalization constant is fixed by \(g(0)=10I_5\).

## 3. Theorem W5 — exact pullback of Bergman volume

Use the real-coordinate convention

\[
dV_B=2^5\det(g(z))\,d^{10}z.
\]

Combining the determinant identity with the chart Jacobian gives

\[
\begin{aligned}
\Phi_\beta^*(dV_B)
&=2^5\frac{10^5}{N(\Phi_\beta(q,p))^5}
  20^{-5}\,d\mu_L\\
&=\boxed{N(\Phi_\beta(q,p))^{-5}\,d\mu_L}.
\end{aligned}
\]

The prefactor cancels exactly because \(2^5 10^5/20^5=1\). Therefore the D5 Bergman volume is the Q6 five-mode Liouville measure weighted by the fifth inverse power of the D5 generic norm, on the preimage of the bounded domain.

At the origin, \(N(0)=1\), so the pullback density equals one. Away from the origin it is nonconstant. This is a genuine local phase-space-to-invariant-measure identity, not merely a match of dimensions or a fitted volume ratio.

## 4. What this does and does not establish

**Established by the theorem and its implementation**
- the Maxwell-normalized five-mode chart has frequency \(\sqrt{12}\), not \(\sqrt{12\beta}\);
- its real Liouville Jacobian is exactly \(20^{-5}\);
- the D5 Bergman determinant has the factor \(10^5N^{-5}\);
- the Bergman volume pulls back to \(N^{-5}d\mu_L\);
- the identity is independent of the bare \(\beta\) in its constant Jacobian, although the domain coordinate \(\Phi_\beta\) and its preimage depend on \(\beta\).

**Still open**
- whether the full Q6 Wilson measure, including the rank-124 integer fibre sum, induces this five-mode measure;
- whether the induced measure has a canonical finite Euclidean D5 realization and radius;
- whether a microscopic RGS/PSU action fixes \(\beta\), the U(1) charge normalization and the continuum limit;
- whether any resulting renormalized coupling predicts \(\alpha(0)\).

The complete Bergman volume of the noncompact D5 domain is not the finite Euclidean volume used in Wyler's expression. W5 constructs an invariant Bergman measure, not the missing finite-volume prescription or the scale normalization. It therefore does not yet derive Wyler's number.

## 5. Relation to Gödel and Weyl

The result gives a precise geometric bridge from the Q6 five-mode sector to a D5 invariant measure. It does not turn Gödel's modal positivity into a probability measure, nor does it establish necessary existence. Hermann Weyl's gauge-geometric ideas, Armand Wyler's bounded-domain volume expression, and Gödel's modal argument remain distinct mathematical inputs. Their connection is now constrained by an explicit measure identity rather than a numerical analogy.

## 6. Reproducibility

Run:

    pytest -q tests/test_rgs_q6_d5_liouville_measure.py

The tests check the Bergman determinant, the pullback density, the origin normalization, and the Maxwell Hamiltonian identity over multiple positive \(\beta\) and eigenvalue choices.
