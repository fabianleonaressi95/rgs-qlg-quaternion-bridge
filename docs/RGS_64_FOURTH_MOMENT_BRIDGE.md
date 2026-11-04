# RGS 64-site fourth-moment bridge

## Why this is the bridge

The recovered RGS seed is the actual hash-verified (64\times3) coordinate array
from `load_rgs_seed_data()`, not the 120-point 600-cell coordinate set.
Its unnormalised second moment is

\[
M^{(2)}_{ij}=\sum_{n=1}^{64}X_{n,i}X_{n,j}=80\,\delta_{ij}.
\]

The bridge is the **rank-four moment** of this same seed:

\[
M^{(4)}_{ijkl}=\sum_{n=1}^{64}
X_{n,i}X_{n,j}X_{n,k}X_{n,l}.
\]

For (d=3), the orthogonal projection onto the SO(3)-invariant rank-four
subspace is

\[
M^{(4),\mathrm{iso}}_{ijkl}
=c(\delta_{ij}\delta_{kl}+\delta_{ik}\delta_{jl}
+\delta_{il}\delta_{jk}),\qquad
c=\frac{M^{(4)}_{iikk}}{3(3+2)}.
\]

The anisotropic part is

\[
R_{ijkl}=M^{(4)}_{ijkl}-M^{(4),\mathrm{iso}}_{ijkl}.
\]

This is not obtained by simply setting the cross component (M_{1122})
as the isotropic coefficient. The Frobenius-orthogonal SO(3) projection
is fixed by the double trace. Consequently (R_{iikl}=0), and (R) is
orthogonal to the isotropic tensor.

## The linear constitutive bridge

The fourth-rank residual defines a canonical map on symmetric traceless
rank-two tensors:

\[
\mathcal{A}:\operatorname{Sym}_0(3)\to\operatorname{Sym}_0(3),
\qquad
(\mathcal{A}Q)_{ij}=R_{ijkl}Q_{kl}.
\]

It is equivariant under the recovered cubic group (O_h). Under the
standard cubic decomposition

\[
\operatorname{Sym}_0(3)=E_g\oplus T_{2g},
\]

the response splits into a two-dimensional (E_g) sector and a
three-dimensional (T_{2g}) sector. The audit computes the full (5\times5)
matrix in a Frobenius-orthonormal basis, its eigenvalues, and the
off-block mixing residual.

This gives a mathematically defined bridge from the discrete fourth moment
to a rank-two anisotropic constitutive response. It is a candidate input
for an RGS metric/coframe constitutive law; it is **not yet** a vierbein
equation or a torsion tensor.

## Important tensor-index obstruction

(R_{ijkl}) is fully symmetric in all four indices. Teleparallel torsion
(T^a{}_{\mu\nu}) is antisymmetric in (mu,\nu). Therefore no direct
identification (T=R) is admissible. A torsion construction would require
a separately derived differential map involving a coframe, its derivatives,
and a compatible connection, with the antisymmetry and local Lorentz/gauge
transformation laws checked explicitly.

## Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_fourth_moment_bridge.py -q
PYTHONPATH=. python -m src.rgs_fourth_moment_bridge --output results/json/rgs_64_fourth_moment_bridge.json
```

## Evidence status

- **NUMERICALLY VERIFIED:** the recovered 64-site seed has (M^{(2)}=80I);
  the fourth-moment SO(3) projection, trace-free residual, (O_h) invariance,
  and (E_g/T_{2g}) response split are independently tested.
- **CANDIDATE:** the map (Q\mapsto R:Q) as an anisotropic constitutive
  response feeding a future coframe/metric law.
- **OPEN:** derive the actual coframe field equation, construct torsion with
  correct index symmetries, couple the response to SU(2)/Yang--Mills/Q6,
  and establish a continuum/refinement or Mosco limit.

The 600-cell audit in `RGS_PSU_Equidistribution_Proof.py` is a distinct
quaternionic/Haar finite-moment check. It is not a substitute for this
64-site RGS fourth-moment bridge.
