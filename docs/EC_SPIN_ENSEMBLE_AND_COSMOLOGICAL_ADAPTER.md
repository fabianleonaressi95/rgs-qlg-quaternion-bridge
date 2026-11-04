# Einstein–Cartan spin ensemble: finite-density adapter

## 1. What is being averaged

The Dirac module fixes eta=(-,+,+,+) and the real axial current
J5^a=i bar(psi) gamma5 gamma^a psi. Its action-derived Cartan solution is
T_abc=-(kappa/2) epsilon_abcd J5^d.

For interoperability, this module defines the **Cartan-normalized source**
Sigma_abc := -(1/2) epsilon_abcd J5^d, so T_abc=kappa Sigma_abc.
Sigma is a rank-three tensor with the normalization needed by the Cartan
equation. It is intentionally not asserted to equal a canonical Noether spin
current without a separate variation and normalization audit.

## 2. Ensemble moments: average before closing the fluid

For samples J_(n)^a with normalized weights w_n, the module computes
< J^a > = sum_n w_n J_(n)^a,
C^{ab}=< (J^a-<J^a>)(J^b-<J^b>) >,
and M^{ab}=<J^a J^b> = <J^a><J^b>+C^{ab}.

The Lorentz scalar entering the local contact operator is
< J_a J^a > = eta_ab M^{ab}, not generally
< J >_a < J >^a. An unpolarized ensemble can have <J>=0 while retaining
nonzero axial-current fluctuations and therefore a nonzero mean contact
operator. This is the key distinction between a vanishing macroscopic
polarization and vanishing microscopic spin correlations.

The code also evaluates the spatial isotropy residual in a declared local
orthonormal frame. Statistical isotropy is a condition on the spatial
second moment, not merely on the mean current.

## 3. What the cosmological adapter does not assume

The adapter records a, H, matter rho and p, units, conventions, derivation
identifier, and the spin ensemble. It deliberately leaves rho_spin,
effective rho/p, and TorsionLambda beta unset. The contact Lagrangian alone
does not determine an effective perfect-fluid energy density: one must vary
the full effective matter action with respect to the tetrad/metric, specify
the quantum or statistical state, and derive its equation of state and
conservation law. Nor can a phenomenological beta be read off without
matching the precise TorsionLambda action and definitions.

Minimal EC has algebraic torsion and no propagating torsion modes. A
propagating torsion cosmology requires additional kinetic terms and a new
action/constraint analysis.

## 4. Validation and evidence status

- ACTION-DERIVED INPUT: Cartan relation and current convention inherited
  from rgs_ec_algebraic_torsion.
- NUMERICALLY VERIFIED: weighted moments, covariance decomposition,
  Cartan source reconstruction, torsion average, contact average, and
  isotropy diagnostics.
- OPEN: canonical spin-current normalization from the complete Noether
  variation; tetrad variation of the effective four-fermion action;
  isotropic fermion-state equation of state; Bianchi/conservation reduction;
  comparison with TorsionLambda beta; mapping to RGS/TEGR variables.

No identification is made between the Dirac spin density and a Q-ball U(1)
Noether charge, electric current, or MHD magnetic field. Those require
separate currents and coupling terms in a common action.

## 5. Units

Every adapter record requires explicit units_id and derivation_id.
The numerical routines do not convert units. The user must supply kappa and
J5 in a mutually consistent system; otherwise dimensional comparisons are
invalid. The frame is local orthonormal Lorentz frame with signature
(-,+,+,+).


## Isotropy refinement

The FLRW test is applied to the full second moment M^{ab}=<J5^a J5^b>,
not to the covariance alone. The audit also checks the coherent spatial
mean. A vanishing covariance anisotropy cannot hide a nonzero mean spatial
polarization, which would select a preferred direction.
