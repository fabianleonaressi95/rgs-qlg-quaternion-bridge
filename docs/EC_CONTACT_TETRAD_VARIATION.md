# Fixed-state tetrad variation of the EC axial contact term

## Scope and conventions

Use eta_ab=diag(-1,1,1,1), kappa=8 pi G (c=1), and the repository's real
axial tangent-frame current J5^a=i bar(psi) gamma5 gamma^a psi. The minimal
EC contact term is

    L4f = alpha J5_a J5^a,   alpha=3 kappa/16.

Let M^{ab}=<J5^a J5^b> be the local-frame second moment. The ensemble
average is taken before tetrad variation, and the spinor components and
M^{ab} are held fixed in the local orthonormal frame during the explicit
coframe variation.

## Explicit variation

The interaction action is S4f=int d4x e alpha eta_ab M^{ab}. At fixed local
state, M^{ab} and eta_ab have no explicit coframe dependence; only
e=det(e^a_mu) varies:

    delta e = e e_a^mu delta e^a_mu.

With the convention delta S_m = int d4x e T^mu_a delta e^a_mu, this gives

    T^mu_a|4f,explicit = e_a^mu <L4f>,
    T_munu|4f,explicit = g_munu <L4f>.

In the local orthonormal frame,

    <L4f> = (3 kappa/16) eta_ab M^{ab},
    rho_4f,explicit = T_00 = -<L4f>,
    p_4f,explicit = <L4f>.

Thus the *explicit fixed-state contact contribution* has vacuum form
p=-rho. This is not yet the equation of state of the complete fermion
ensemble or of cosmology: if M^{ab}, occupation numbers, spinor amplitudes,
or the state itself respond to the tetrad/scale factor, their variation
contributes additional terms. Kinetic Dirac stress, constraints, and
conservation must be included before interpreting this as a cosmological
fluid. In particular, this result does not determine TorsionLambda beta.

## Isotropy condition

For an FLRW-compatible local-frame ensemble, it is not enough that the
covariance C^{ij} be isotropic. The full second moment must satisfy
M^{0i}=0 and M^{ij}=q delta^{ij}; a nonzero mean spatial polarization
bar(J)^i also selects a direction. The audit therefore reports separately
the second-moment anisotropy, mixed time-space second moment, and mean
spatial polarization.

## Evidence status

- PROVED: the displayed explicit coframe variation under the declared
  fixed-local-state prescription.
- NUMERICALLY VERIFIED: tensor construction, Lorentzian contraction,
  and validation tests.
- OPEN: state-response contribution, full Dirac-plus-contact stress tensor,
  thermodynamic equation of state, covariant conservation/Bianchi reduction,
  and matching to the exact TorsionLambda action and beta definition.

This is a controlled variational component, not a derivation of a complete
spin-fluid cosmology.
