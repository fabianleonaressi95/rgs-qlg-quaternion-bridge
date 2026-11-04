# Explicit tetrad variation of the EC axial contact operator

## Scope and conventions

Use eta_ab=diag(-1,+1,+1,+1), kappa=8 pi G (c=1), and the real axial
current convention already fixed by the EC Dirac module:
J5^a = i bar(psi) gamma5 gamma^a psi.
The effective operator is
L_4f = c4 J5_a J5^a,  c4=3 kappa/16.

This audit varies the contact term in the tetrad formalism while holding
the local Lorentz spinor components fixed. It is an explicit operator
variation, not yet the complete on-shell Hilbert tensor of the interacting
Dirac system.

## Tetrad variation

Let e^a_mu be an oriented tetrad, e=det(e^a_mu)>0, and
g_mu_nu=eta_ab e^a_mu e^b_nu. In a local orthonormal frame the components
J5^a are independent of the tetrad under the declared fixed-spinor
variation. Therefore L_4f=c4 eta_ab J5^a J5^b has no explicit tetrad
dependence, and
  delta(e L_4f) = e L_4f e_a^mu delta e^a_mu.

Thus
  (1/e) delta S_4f/delta e^a_mu = L_4f e_a^mu.

Equivalently, since delta e/e = -(1/2) g_mu_nu delta g^mu_nu and with
T_mu_nu=-2/sqrt(-g) delta S/delta g^mu_nu,
  T_mu_nu^(4f, explicit) = g_mu_nu L_4f.

The code returns this tensor and checks the tetrad derivative against a
central finite difference of det(e)*L_4f. It also reports the metric and
inverse metric, allowing independent contraction and frame checks.

## Critical limitation: this is not yet a cosmological fluid

It would be incorrect to identify this isolated explicit contribution
with the complete on-shell stress tensor or to infer a cosmological
equation of state from it. The full effective Dirac action also contains
the tetrad-dependent kinetic term, and the spinor state/current changes
along solutions. The total Hilbert tensor must be varied before imposing
the homogeneous state and its equations of motion. Kinetic and contact
contributions can exchange energy; the contact term alone need not be a
separately conserved perfect fluid.

In particular, the often-used thermodynamic expression
p=n d rho/dn-rho requires a justified conserved number density n and a
state/ensemble closure. It is not obtained by simply relabelling
T_mu_nu^(4f, explicit). No rho_spin, w_spin, beta, or Friedmann correction
is assigned here.

## Required next derivation

1. Vary the full torsion-eliminated Dirac action, including the symmetrized
   kinetic term and the four-fermion operator, with respect to e^a_mu.
2. Derive the spinor equations and use them only after obtaining the
   Hilbert tensor; track possible boundary terms.
3. Construct a homogeneous/isotropic fermionic state and compute both
   <J5^a> and <J5^a J5^b>, including connected correlations.
4. Derive the total rho and p from the full tensor, and verify
   covariant conservation/Bianchi identities.
5. Match conventions and variables against the exact TorsionLambda action.
   Only then can a beta-like coefficient be derived.

## Evidence status

- ACTION-DERIVED: the contact coefficient and current convention inherited
  from the prior EC derivation.
- NUMERICALLY VERIFIED: the explicit tetrad derivative against central
  finite differences; Lorentz scalar invariance under spatial frame rotations.
- OPEN: complete on-shell Dirac Hilbert tensor, state closure, equation of
  state, conservation reduction, and TorsionLambda matching.

Minimal Einstein-Cartan torsion remains algebraic in this construction.
This module introduces no propagating torsion modes and makes no
identification of axial spin current with electric/Q-ball charge or MHD
current.
