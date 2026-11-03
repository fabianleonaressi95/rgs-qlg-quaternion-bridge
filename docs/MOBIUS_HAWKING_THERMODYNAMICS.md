# RGS-Q6 Mobius boundary: fermionic and modular thermodynamics

## Scope

This phase implements reproducible finite-spectrum thermodynamics from the
24 naive and 24 Wilson boundary levels and the corresponding 24 eigenvalues
of the restricted Gaussian correlation matrix. It does not recompute the
Dirac-Wilson matrix, establish a Hawking horizon, or derive a continuum
Hawking flux. The temperature T_H=1/L_b is an input scale from the model.

## Fermion determinant and vacuum normalization

For independent fermionic levels E_k at chemical potential mu=0,

    Z_raw(beta) = product_k (1 + exp(-beta E_k)).

The spectrum includes negative and positive Dirac levels. With the filled
negative-energy sea as reference,

    E_vac = sum_{E_k<0} E_k
    Z_rel = exp(beta E_vac) Z_raw
    U_rel = sum_k E_k/(exp(beta E_k)+1) - E_vac
    S_th = log(Z_raw) + beta sum_k E_k/(exp(beta E_k)+1)
    C_V = beta^2 sum_k E_k^2 f_k(1-f_k)
    f_k = 1/(exp(beta E_k)+1).

The implementation returns raw and vacuum-normalized partition functions.
Repeated eigenvalues already encode degeneracy; no additional multiplicity
is applied.

## Modular partition and entanglement entropy

For correlation eigenvalues 0<=zeta_k<=1,

    epsilon_k = log((1-zeta_k)/zeta_k)
    Z_A_mod = product_k (1 + exp(-epsilon_k))
            = product_k (1-zeta_k)^(-1)
    S_EE = -sum_k [zeta_k log(zeta_k)
                  +(1-zeta_k) log(1-zeta_k)].

Only roundoff-scale excursions at endpoints are projected. Exact zeta=0
modes have zero entropy and epsilon=+infinity; exact zeta=1 modes have zero
entropy and epsilon=-infinity, while the unshifted modular product diverges.
Infinite modular energies are encoded as JSON null plus explicit endpoint
counts, avoiding non-standard JSON Infinity.

The modular partition represents the reduced Gaussian state. It is not an
independent thermal sector and is not multiplied by Z_rel without a separate
derivation.

## Reference calculation

Input: L_b=37.50208453459409, T_H=0.026665184413350203, 24 levels in each
spectrum, and 24 correlation eigenvalues in each sector.

| Quantity | Naive | Wilson |
|---|---:|---:|
| E_vac | -5.17341132985 | -1.37921562633 |
| log Z_rel | 0.10186763291 | 0.21726633396 |
| Z_rel | 1.10723690030 | 1.24267502462 |
| U_rel | 0.01176147366 | 0.01627531430 |
| S_th | 0.54294741231 | 0.82762454669 |
| C_V | 1.89809918785 | 1.66973717275 |
| F_rel | -0.00271631922 | -0.00579344686 |
| S_EE | 9.23887876499 | 7.53353795513 |

The small Wilson difference between reconstructed S_EE and the reported
7.53353795536 is due to the rounded eigenvalues supplied in the reference
data.

The supplied Hawking-power values are carried into the report as imported
diagnostics. They are not derived from this partition function: a flux
calculation additionally requires mode normalization, transmission
coefficients/greybody factors, and a precise definition of the radiative
observable.

## Reproduction

    PYTHONPATH=. python examples/mobius_hawking_thermodynamics.py
    PYTHONPATH=. python -m pytest tests/test_thermodynamics.py -v

Default input: results/json/mobius_hawking_spectrum_reference.json
Default output: results/json/mobius_hawking_thermodynamics.json

## Evidence boundary

- NUMERICALLY VERIFIED: finite-spectrum fermion partition functions,
  vacuum-normalized thermodynamic observables, modular partition and
  correlation-matrix entanglement entropy for the supplied 24-mode arrays.
- NOT ESTABLISHED: derivation of the Hawking temperature from a genuine
  Killing horizon / Bisognano-Wichmann setting; continuum Hawking radiation;
  greybody-corrected flux; Einstein-Cartan backreaction; a coupled
  gravitational path integral; and a continuum limit.
