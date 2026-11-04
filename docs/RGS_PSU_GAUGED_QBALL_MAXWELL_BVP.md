# RGS–PSU: stationary gauged Q-ball coupled to Maxwell

## Scope

`src/rgs_psu_qball_maxwell_bvp.py` solves the stationary, spherically
symmetric, electrostatic Maxwell–complex-scalar system in flat spacetime.
It is the next layer after the global-`U(1)` BVP in
`src/rgs_psu_qball_bvp.py`. The gauge field is dynamical and sourced by
the scalar Noether current; it is not an externally prescribed potential.

The model uses signature `(-,+,+,+)`,
`D_mu = partial_mu - i q A_mu`,
`Phi=f(r) exp(i omega t)`, and `A=A_0(r) dt`. Define
`Omega(r)=omega-q A_0(r)`. With
`V(f)=m^2 f^2-lambda f^4+g f^6`, the radial equations are

```text
f'' + 2 f'/r = (m^2-Omega^2) f - 2 lambda f^3 + 3 g f^5
A0'' + 2 A0'/r = -2 q Omega f^2
```

The Maxwell sign follows from `nabla_mu F^{mu nu}=j^nu` and
`j^0=2 q Omega f^2` in these conventions. Reversing the convention for
`D_mu` or the sign of the scalar phase reverses the corresponding charge
signs; the equations and charge diagnostics must be changed together.

## Boundary-value problem

The solver imposes
`f(0)=f0`, `f'(0)=0`, `f(Rmax)=0`,
`A0'(0)=0`, and `A0(Rmax)=0`. The last condition fixes the residual
constant-gauge freedom by choosing the electrostatic potential to vanish
at the finite outer cutoff. The frequency `omega` is the BVP eigenparameter.
The two radial (2/r) terms are handled with the regular-singular matrix in
the collocation solver; the code does not divide by (r=0).

The phase/particle Noether charge and electric charge are distinct:
```text
Q_particle = 8 pi integral_0^Rmax r^2 Omega(r) f(r)^2 dr
Q_electric = q Q_particle
```
The total energy includes both the scalar and electrostatic field:
```text
E = 4 pi integral_0^Rmax r^2 [
      Omega(r)^2 f^2 + f'^2 + V(f) + (1/2) A0'^2
    ] dr
```
The implementation reports scalar and Maxwell energy separately, both
equations' independently evaluated residuals, and a Gauss-law check between
the electric flux at the boundary and (q Q_{particle}).

## Validation strategy

The tests cover:
- the (q=0) decoupling limit against the existing global-(U(1)) BVP;
- a nonzero-coupling solution with both scalar and Maxwell residual checks;
- origin and outer boundary conditions;
- Gauss-law flux/charge consistency and positive field energy;
- invalid input rejection.

A nonzero-(q) result is a finite-radius numerical stationary solution, not
a proof of existence on all of (mathbb R^3), and not a stability result.
The necessary asymptotic frequency window is checked using (A_0(Rmax)=0);
it is not a sufficient existence or stability theorem.

## RGS vierbein/TEGR extension: next layer, not included here

For a static spherical metric
`ds^2=-N(r)^2dt^2+A(r)^2dr^2+r^2dOmega_2^2`, the scalar equation to be
coupled to the curved Maxwell equation is
```text
(1/(N A r^2)) d/dr[(N r^2/A) f']
  + Omega^2 f/N^2 - V'(f)/2 = 0
```
and the Maxwell equation must use the same metric volume form and inverse
metric in `nabla_mu F^{mu nu}=j^nu`. In TEGR, (N,A) must ultimately be
obtained from the chosen tetrad/field equations (or a controlled fixed
background) and the same action. The present code does not insert a
phenomenological torsion coupling and does not claim teleparallel
backreaction.

## Epistemic status

**Implemented target:** flat, spherical, electrostatic gauged-(U(1))
Q-ball with self-consistent Maxwell (A_0).

**Not implemented:** magnetic Maxwell components, rotating/non-spherical
solutions, Q-ball/fluxon bound states, curved vierbein/TEGR backreaction,
MHD evolution, fluctuation spectrum, linear or nonlinear stability, and a
continuum/Mosco limit.
