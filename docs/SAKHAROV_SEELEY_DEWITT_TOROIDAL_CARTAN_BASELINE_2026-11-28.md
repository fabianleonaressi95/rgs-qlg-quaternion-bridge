# Seeley–DeWitt calculation baseline — implementation contract

Primary derivation: theorem/research/SAKHAROV_SEELEY_DEWITT_TOROIDAL_CARTAN_BASELINE_2026-11-28.md.

## Geometry and conventions
- Euclidean toroidal coframe \(e^0=Nc\,d\tau,\ e^1=dl,\ e^2=r\,d\theta,\ e^3=(R+r\cos\theta)d\varphi\).
- Volume density \(\sqrt g=Ncr(R+r\cos\theta)\).
- Independent connection \(\omega=\omega_{\rm LC}+K\); do not set \(K=0\) in the target calculation.
- For every fluctuation operator store \(P=-(\mathcal D^2+E)\), \(\Omega=[\mathcal D,\mathcal D]\), bundle rank, domain, positivity and boundary/interface data.

## Exact coefficient engine
Use the universal four-dimensional Laplace-type formulas:
\[
a_0=(4\pi)^{-2}\operatorname{tr}I,\quad
a_2=(4\pi)^{-2}\operatorname{tr}(E+R_{\rm LC}/6),
\]
\[
a_4=(4\pi)^{-2}\operatorname{tr}\left[
\frac12E^2+\frac16R_{\rm LC}E+\frac16\mathcal D^2E+
\frac1{30}\Box R_{\rm LC}+\frac1{72}R_{\rm LC}^2
-\frac1{180}R_{\mu\nu}^2+\frac1{180}R_{\mu\nu\rho\sigma}^2+
\frac1{12}\Omega_{\mu\nu}^2\right].
\]
The complete scalar \(m,\xi\) specialization and torsion-free Dirac \(a_2\) check are in the primary derivation.

## Baseline result, not a final model prediction
For \(N_s\) real scalars and \(N_D\) four-component Dirac fermions, no vectors/ghosts yet, massless scalar+torsion-free Dirac check:
\[
C_R=\sum_s(1/6-\xi_s)+N_D/3,\qquad C_0=N_s-4N_D.
\]
With the declared proper-time cutoff \(\epsilon=\Lambda_{\rm UV}^{-2}\) and Euclidean action convention \(-M_{\rm ind}^2\int\sqrt g R/2\):
\[
M_{\rm ind}^2=C_R\Lambda_{\rm UV}^2/(4\pi)^2.
\]
This is explicitly a minimal baseline, regulator-dependent, and not the PSU model's final Newton constant.

## Blocking items before a physical number
1. Fix the actual \(K^{ab}(l,\theta)\) profile and compute \(E_D(K),\Omega_D(K)\).
2. Fix physical matter multiplets and all multiplicities/statistics.
3. Derive gauge/gravitational gauge-fixing and ghost operators.
4. Specify smooth throat gluing or shell/interface matching and domains.
5. Compute the full metric curvature and torsion invariants for the selected \(N(l),r(l),R(l)\).
6. State \(\Lambda_{\rm UV}\), renormalization conditions and analytic continuation.
7. Only then report \(G_{\rm ind}\), \(\Lambda_{\rm ind}\), induced torsion couplings and solve semiclassical backreaction.

Boundary/interface heat coefficients are not to be dropped when the throat is treated as a shell. The field supertrace is boson minus fermion with ghost weights; it is not the 64/124 matrix grading and not Anti-Twin sheet exchange.
