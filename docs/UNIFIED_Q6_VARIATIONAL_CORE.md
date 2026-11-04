# Q6 unified variational core — status and limits

This implementation adds a finite-dimensional Euclidean action on the canonical
64-vertex Q6 / 192-edge / 240-face complex.

## Fields

- `e_v`: a real, invertible 4x4 coframe at each vertex;
- `U_uv in SU(2)`: one link variable per oriented Q6 edge;
- `psi_+, psi_- in C^2`: two fundamental SU(2) fields at every vertex.

The link convention is the existing repository convention:
`U_uv -> G_v U_uv G_u^dagger`, `psi_v -> G_v psi_v`.

## Action

With `V_v = |det(e_v)|`, the prototype is

\[
S_h = S_e + S_W + S_{\rm hop} + S_{\rm twin} + S_4,
\]
\[
S_e = \frac{\kappa_e}{2}\sum_{(u,v)\in E}\|e_u-e_v\|_F^2,
\]
\[
S_W = \beta\sum_{p\in F}\left(1-\frac12\operatorname{ReTr}U_p\right),
\]
\[
S_{\rm hop} = \frac{\kappa_\psi}{2}\sum_{(u,v)\in E}
\sqrt{V_uV_v}\sum_{\sigma=\pm}
\|\psi_{\sigma,v}-U_{uv}\psi_{\sigma,u}\|^2,
\]
\[
S_{\rm twin} = \frac{\mu}{2}\sum_v V_v\|\psi_{+,v}-\psi_{-,v}\|^2,
\quad
S_4 = \frac{\lambda}{4}\sum_v V_v
(\|\psi_{+,v}\|^2+\|\psi_{-,v}\|^2)^2.
\]

All couplings are non-negative. The action is real and non-negative. Its
Wilson and matter terms are invariant under the local SU(2) transformations
above; the coframe is an SU(2)-singlet in this prototype.

## What is actually established

The code evaluates the action and reports its five terms separately. Tests
check the zero-field reference configuration, non-negativity, exact graph
support, rejection of degenerate coframes, and local SU(2) gauge invariance
to floating-point tolerance.

## Deliberate limits

1. `S_e` is a fixed-frame coframe smoothness regulariser, not
   Einstein-Hilbert or TEGR gravity.
2. `S_hop` is a covariant graph energy for SU(2) doublets, not
   a Dirac operator or a fermionic kinetic term.
3. The model is Euclidean and finite; it does not define Lorentzian dynamics.
4. No spin connection `omega` or Clifford-valued transport is yet
   derived from the coframe.
5. The implementation evaluates the action but does not yet solve its coupled
   Euler-Lagrange equations or prove a continuum/Mosco limit.
6. PSU sheets are represented as two fields with a local locking interaction;
   this is a candidate encoding, not evidence of a physical anti-twin sector.

## Next mathematical gate

Derive and independently test the constrained variations with respect to
`e_v`, `U_uv`, `psi_+`, `psi_-`. Then introduce a coframe-derived
spin connection and Clifford transport without breaking the local gauge
covariance already tested here. Only after that should the Euclidean prototype
be extended to a Lorentzian TEGR/Dirac system.
