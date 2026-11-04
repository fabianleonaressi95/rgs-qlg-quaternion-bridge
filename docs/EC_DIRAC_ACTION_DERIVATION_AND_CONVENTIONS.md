# Einstein–Cartan–Dirac: action, axial current and Cartan coefficient

## Fixed conventions

This implementation fixes the conventions of Bonder and Corral,
*Unimodular Einstein–Cartan gravity: Dynamics and conservation laws*,
Phys. Rev. D 97, 084001 (2018), arXiv:1712.09903. The first-order setup uses
eta_ab=diag(-1,+1,+1,+1), independent vierbein e^a and Lorentz connection
omega^ab, with
    T^a = de^a + omega^a_b wedge e^b
    R^ab = d omega^ab + omega^a_c wedge omega^cb.
For zero cosmological constant:
    S_g = (1/(4 kappa)) integral epsilon_abcd R^ab wedge e^c wedge e^d,
    kappa = 8 pi G (c=1).
The symmetrized Dirac action is
    S_D = 1/2 integral [bar(psi) *gamma wedge D psi
                         + D bar(psi) wedge *gamma psi],
    D psi = d psi + (1/4) omega^ab gamma_ab psi.
The bilinear conventions are
    bar(psi)=i psi^dagger gamma^0,
    {gamma^a,gamma^b}=2 eta^ab,
    gamma5=i gamma^0 gamma^1 gamma^2 gamma^3,
    J5^a=i bar(psi) gamma5 gamma^a psi.
The alternating internal symbol is epsilon_0123=+1, not the metric-raised
Levi-Civita tensor.

## Connection variation and coefficient

The Lorentz-connection equation is
    E_ab = kappa sigma_ab,
    E_ab = epsilon_abcd T^c wedge e^d,
    sigma_ab = (1/2) J5 wedge e_a wedge e_b.
Its algebraic solution is
    T^a = (kappa/2) * (J5 wedge e^a).
Writing T^a=(1/2)T^a_bc e^b wedge e^c and
*(e^d wedge e^a)=(1/2)epsilon^da_bc e^b wedge e^c, lowering the first
Lorentz index with eta=(-+++) gives
    T_abc = -(kappa/2) epsilon_abcd J5^d.
Therefore the coefficient in T_abc=c_T*kappa*epsilon_abcd*J5^d is
    c_T = -1/2.
This replaces the former provisional -1/4. Factors differing in other
sources must be traced to their definitions of spin density, torsion,
axial current, epsilon orientation and gravitational coupling.

## Executable spinor current

The module supplies explicit 4x4 matrices for the (-+++) Clifford algebra,
the stated Dirac adjoint, gamma5, and the bilinear J5^a. The function
axial_torsion_from_dirac_spinor computes the current from a spinor and inserts
it into the Cartan solution. Tests check the Clifford algebra, gamma5 squared,
current reality, zero-spinor limit and the Cartan residual.

## Status and limits

- PROVED within the declared first-order action: the Cartan equation is
  algebraic; minimal EC has no propagating torsion degrees of freedom.
- NUMERICALLY VERIFIED: Clifford algebra and spinor-current-to-torsion tests.
- CANDIDATE: bridge to the RGS spinor/Kirchhoff-current construction.
- OPEN: derive the full contorsion with identical index ordering; substitute
  it into the complete Dirac action and independently verify the induced
  axial-axial four-fermion term; then construct a fermionic spin-fluid average
  and compare with TorsionLambda. Q-ball U(1) charge is not spin density.
  No EC-to-TEGR connection map is asserted.

## Reference

Y. Bonder and C. Corral, Phys. Rev. D 97, 084001 (2018), arXiv:1712.09903.
See the first-order action, connection field equation, Dirac spin density
and axial current, and the Cartan solution T^a=kappa * (J5 wedge e^a)/2.
