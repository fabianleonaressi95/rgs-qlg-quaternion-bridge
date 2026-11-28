# Sakharov induced Cartan gravity — implementation map

This companion record connects the unified PSU–RGS–Anti-Twin wormhole implementation to the research specification in theorem/research/SAKHAROV_SUPERTRACE_SUSY_PSU_ANTITWIN_INDUCED_CARTAN_GRAVITY_2026-11-28.md.

## Operator separation (mandatory)

Keep three data products separate:
1. Q6 graph/covariant operator and its heat trace \(K_L(s)=\mathrm{Tr}(e^{-sL})\);
2. the geometric 188D RGS operator \(K_0\), which has negative eigenvalues and does not commute with the proposed 64/124 grading;
3. the algebraic SUSY candidate \(H=Q^2\), with \(\Gamma=\mathrm{diag}(I_{64},-I_{124})\), finite index \(-60\), and \(\mathrm{Str}_\Gamma(e^{-sH})=-60\).

Do not identify any pair of these operators. Do not equate Anti-Twin sheet exchange \(\tau_x\) with fermion parity.

## Required software interfaces

- A typed operator record must include Hilbert-space dimension, inner product, domain, self-adjointness/positivity evidence, zero modes, grading (if any), geometric background, and provenance/checksum.
- A heat-kernel audit must report \(K(s)\), zero-mode multiplicity, large-\(s\) limit, spectral dimension \(d_s=-2\,d\log K/d\log s\), cutoff dependence and refinement scaling.
- A field-supertrace registry must list every boson, fermion, gauge mode and ghost, with statistics, multiplicity, operator, normalization and boundary conditions.
- The continuum operator factory must accept the same independent coframe \(e\) and Lorentz connection \(\omega\) as the Palatini–Cartan wormhole model. Torsion \(T=de+\omega\wedge e\) must remain available to the operator coefficients.
- The two-sheet operator must use throat transport derived from the PSU/Anti-Twin gluing data; sheet-even/odd decomposition is not a boson/fermion decomposition.
- The Euclidean heat-kernel calculation must include throat/interface and asymptotic-boundary contributions. Lorentzian results require a stated analytic continuation.

## Sakharov output contract

For the gauge-fixed Euclidean fluctuation operators \(\mathcal O_i[e,\omega]\), compute the field supertrace of their heat kernels and the Seeley–DeWitt coefficients \(a_0,a_2,a_4\). The induced action is to be reported in the form
\[
\Gamma_{\mathrm{eff}}=\int e\left[-\Lambda_{\mathrm{ind}}+
\frac{M_{\mathrm{ind}}^2}{2}\mathcal R(\omega,e)
+c_T\mathcal I_T(T)+c_1\mathcal R^2+c_2R_{\mu\nu}R^{\mu\nu}+\cdots\right]
+\Gamma_\Sigma.
\]
Every coefficient must include its regulator, field-content and renormalization-condition dependence. No \(G_{\mathrm{ind}}\) is to be reported until the actual \(a_2\) calculation and normalization are complete.

## Workstreams

- WQ6-09: Q6 spectral-to-continuum bridge.
- WQ6-10: Sakharov heat-kernel coefficients and induced Newton/cosmological couplings.
- WQ6-11: two-sheet spectral matching and throat boundary conditions.
- WQ6-12: physical SUSY realization of the 188D algebraic candidate.
- WQ6-13: torsionful induced Cartan terms and throat support.
- WQ6-14: semiclassical wormhole backreaction from \(S_{\mathrm{PC}}+\Gamma_1+S_\Sigma\).

All six are open workstreams. Existing finite-dimensional audits are inputs for reproducibility, not completed continuum calculations. The companion research specification gives the numerical historical values and evidence classifications.
