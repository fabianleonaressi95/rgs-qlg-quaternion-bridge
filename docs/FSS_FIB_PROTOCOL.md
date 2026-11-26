# Fibonacci Aubry–André finite-size scaling protocol

## Scope and evidence status

This protocol specifies reproducible finite-chain calculations. It does **not**
establish a thermodynamic-limit collective mode, intrinsic damping, or an exact
interacting many-body result. The direct density-channel RPA is an approximation
built on the non-interacting single-particle spectrum.

The CSV in data/legacy_fss_2026-08-15.csv transcribes 12 values from a report
dated 15 August 2026. They are reported historical values, not independently
reproduced by this implementation.

## Model and size sequence

For sites j=0,...,N-1,
\[
H_0=-t\sum_j(c^\dagger_{j+1}c_j+\mathrm{h.c.})
+V\sum_j\cos(2\pi\alpha j+\theta)n_j.
\]
The critical choice is V=2t. Record boundary conditions, phase, filling,
temperature and momentum for every run. Open boundaries are the code default.

For rational Fibonacci approximants, N=F_k and alpha_k=F_{k-1}/F_k.
Do not silently mix this with irrational alpha=1/phi at Fibonacci N. The
historical report states N=55,144,377,987 and alpha=1/phi; its complete
numerical conventions are not available for independent reconstruction.

## Density response and RPA

The one-particle density vertex is rho_q=diag(exp(i q j)). The retarded
non-interacting response is
\[
\chi_0^R(q,z)=\frac1N\sum_{mn}
\frac{(f_m-f_n)|\langle m|\rho_q|n\rangle|^2}
{z+E_m-E_n},\quad \operatorname{Im}z>0.
\]
The occupation difference is essential. At q=0, rho_0 is the conserved total
particle number and the finite-frequency response vanishes.

For spinless nearest-neighbour interaction U sum_j n_j n_(j+1), the direct
density-channel Fourier kernel is V(q)=2U cos(q), with each bond counted once.
The direct-channel RPA is
\[
\chi_{\rm RPA}=\frac{\chi_0}{1-V(q)\chi_0}.
\]
This is not the exact interacting response and includes no separately derived
exchange or self-energy correction.

At a simple pole z_L satisfying 1-V(q)chi_0(z_L)=0, the residue is
\[
\operatorname{Res}_{z=z_L}\chi_{\rm RPA}
=-\frac{1}{V(q)^2\,\partial_z\chi_0(z_L)}.
\]
A real-axis sign change of Re[1-V chi_0] alone is not a reliable complex-pole
locator. Pole tracking must distinguish a zero of the RPA denominator from
a pole of chi_0.

## Broadening and linewidth

The historical report used eta=0.02 and tabulated gamma_L approximately
0.0200–0.0202. This coincidence makes the reported linewidth potentially
resolution-dominated; it is not evidence of intrinsic damping.

A physical linewidth claim requires an eta-to-zero sweep, consistent mode
tracking, finite-size scaling at each broadening, and a stated thermodynamic
limit or bath/self-energy prescription. A finite closed chain has a discrete
spectrum; a fixed numerical broadening is not an intrinsic linewidth.

## Finite-size scaling

A descriptive candidate fit is y_N=y_inf+a N^{-p}. Fix p from independent
theory or use explicit model comparison and sensitivity analysis. Four sizes
do not justify a robust three-parameter asymptotic fit. Report residuals,
fit-window dependence, exponent sensitivity and uncertainty.

The historical values appear close to omega_L approximately U/3 and Z_L
approximately U/3 for the three reported couplings. These are candidate
finite-size trends only, not theorem statements, independently reproduced
limits, or established physical laws.

## Reproduction checklist

Archive N,t,V,alpha,theta; boundary condition; filling/chemical potential,
temperature and occupation rule; q and interaction kernel; eta and frequency
convention; eigensolver and precision; tracked pole/peak identity, root
residual and residue; FSS fit window, exponent, residuals and sensitivity.

## Evidence labels

- Implemented: finite-chain Hamiltonian, occupations, density vertex,
  non-interacting retarded response, direct RPA kernel and fixed-p descriptive fit.
- Unit-testable: algebraic invariants and synthetic fit recovery.
- Not independently reproduced: the 12 historical rows.
- Candidate: omega_L -> U/3 and Z_L -> U/3.
- Open: robust pole tracking, eta -> 0, larger-size FSS, physical
  thermodynamic-limit interpretation.
