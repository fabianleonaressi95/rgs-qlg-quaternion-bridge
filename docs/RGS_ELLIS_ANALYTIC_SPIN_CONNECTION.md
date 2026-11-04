# Analytic LC connection on the Ellis-RGS refinement chart

This increment supplies the missing *analytic first coframe jet* for the
continuum Ellis reference metric used by `rgs_ellis_refinement.py`:

`ds² = -dt² + du² + (u²+r0²)(dθ² + sin²θ dφ²)`.

With coordinate order `(t,u,theta,phi)`, the orthonormal coframe is
`e^a_mu = diag(1,1,R,R sin(theta))`, where `R=sqrt(u²+r0²)`.
Its only nonzero derivatives are

- `partial_u e^2_theta = u/R`;
- `partial_u e^3_phi = (u/R) sin(theta)`;
- `partial_theta e^3_phi = R cos(theta)`.

`src/rgs_ellis_spin_connection.py` builds this jet and passes it to the
independent local Levi-Civita kernel. Tests audit Lorentz antisymmetry, the
tetrad postulate and zero torsion.

**Scope boundary:** the angular refinement uses 56 distinct RGS directions
and a piecewise-linear simplicial embedding in R4. This module evaluates the
exact smooth Ellis chart at continuum coordinates; it does *not* yet
differentiate the 64 site-dependent covariance coframes
`E_i = principal_sqrt(C_i/c0)`, nor does it define a unique interpolation
of their local frame gauge across tetrahedra. That discrete-to-continuum
coframe map remains OPEN. It also does not provide the inertial connection
for TEGR or the independent connection of Einstein-Cartan theory.
