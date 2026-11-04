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

At the origin, \(g(0)=10I_5\). The determinant identity in this normalization is

\[
\boxed{\det g(z)=\frac{10^5}{N(z)^5}}.
\]

It follows by differentiating \(N\) and applying the matrix determinant lemma to the resulting identity-plus-rank-two update; the normalization constant is fixed by \(g(0)=10I_5\).

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

At the origin, \(N(0)=1\), so the pullback density equals one. Away from the origin it is nonconstant.

## 4. Theorem W6 — exact Hua cutoff volume in Q6 Liouville units

Let \(V_E=\pi^5/1920\) be the Euclidean volume of the unit-radius bounded realization, and let \(d\nu_E=d^{10}z/V_E\) be normalized Euclidean measure. The exact D5 Hua density of \(t=N(z)\) is

\[
f_N(t)=\frac{20}{3}(1-\sqrt t)^4(1+4\sqrt t+t),
\qquad 0\le t\le1.
\]

The existing Hua audit proves that

\[
Z_H(\epsilon)=\int_\epsilon^1 t^{-5}f_N(t)\,dt
\]

is exactly

\[
\begin{aligned}
Z_H(\epsilon)={}&\frac53(\epsilon^{-4}-1)
-20(\epsilon^{-3}-1)-30(\epsilon^{-2}-1)\\
&+\frac{128}{3}(\epsilon^{-5/2}-1)
+\frac{20}{3}(\epsilon^{-1}-1).
\end{aligned}
\]

Since \(d\mu_L=20^5d^{10}z=20^5V_Ed\nu_E\), W5 gives

\[
\boxed{
\operatorname{Vol}_B\{N\ge\epsilon\}
=20^5\frac{\pi^5}{1920}Z_H(\epsilon)
}.
\]

In particular,

\[
\boxed{
\operatorname{Vol}_B\{N\ge\epsilon\}
\sim\frac{25000\pi^5}{9}\epsilon^{-4}
\quad(\epsilon\downarrow0).
}
\]

This is the exact bridge between the Q6 five-mode Liouville measure, the D5 Bergman density, and the already established Hua cutoff. It also proves explicitly that the uncut Bergman volume diverges. The cutoff is mathematically controlled, but \(\epsilon\) is not yet a physical scale.

## 5. What this does and does not establish

**Established**
- the Maxwell-normalized five-mode chart has frequency \(\sqrt{12}\), independent of \(\beta\);
- its real Liouville Jacobian is exactly \(20^{-5}\);
- the D5 Bergman determinant is \(10^5N^{-5}\);
- the Bergman volume pulls back to \(N^{-5}d\mu_L\);
- the Hua cutoff volume is exactly the stated \(20^5V_E Z_H(\epsilon)\), with leading divergence \((25000\pi^5/9)\epsilon^{-4}\).

**Still open**
- whether the full Q6 Wilson measure, including the rank-124 integer fibre sum, induces this five-mode measure;
- whether a Wilson-derived regulator selects \(\epsilon\);
- whether the induced measure yields a canonical finite Euclidean D5 realization and radius;
- whether a microscopic RGS/PSU action fixes the U(1) charge and its continuum renormalization;
- whether any resulting coupling predicts \(\alpha(0)\).

The complete Bergman volume of noncompact D5 is not the finite Euclidean volume in Wyler's expression. W5–W6 derive the invariant Bergman measure and its cutoff divergence, not Wyler's finite-volume prescription or its radius normalization. They do not yet derive Wyler's number.

## 6. Relation to Gödel and Weyl

The result gives a precise geometric bridge from the Q6 five-mode sector to a D5 invariant measure. It does not turn Gödel's modal positivity into a probability measure, nor does it establish necessary existence. Hermann Weyl's gauge-geometric ideas, Armand Wyler's bounded-domain volume expression, and Gödel's modal argument remain distinct mathematical inputs. Their connection is now constrained by an explicit measure identity rather than a numerical analogy.

## 7. Reproducibility

Run:

    pytest -q tests/test_rgs_q6_d5_liouville_measure.py

The tests check the Bergman determinant, the pullback density, the origin normalization, the Maxwell Hamiltonian identity, and the exact Hua cutoff-volume factor.
