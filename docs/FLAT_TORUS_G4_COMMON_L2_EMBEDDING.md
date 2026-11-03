# Flat-torus G4: common-Hilbert-space P1 embedding

**Status:** implemented and numerically audited on the flat periodic benchmark.  
**Scope limit:** this is not yet a curved RGS-Q6 embedding or a Mosco theorem.

## 1. Discrete and continuum spaces

For a periodic Freudenthal mesh K_h of the flat torus T_L^d, let
H_h = R^{V_h} and H = L^2(T_L^d).
The embedding J_h: H_h -> H is the continuous piecewise-linear nodal
interpolant. On every simplex sigma=[i_0,...,i_d],
(J_h u_h)|_sigma = sum_a u_{i_a} lambda_a^sigma.
Periodic simplices are evaluated in locally unwrapped coordinates and
identified modulo L. The implementation is J_h(...) in
src/flat_torus_embedding.py.

The map is injective: if J_h u_h=0 in L^2, its continuous P1
representative vanishes everywhere, hence at every mesh vertex, so
u_h=0.

## 2. Exact discrete Hilbert metric

The inner product on H_h is not the lumped G3 quadrature. It is the
exact L^2 inner product of the P1 reconstructions:
  <u_h,v_h>_H_h = <J_h u_h,J_h v_h>_L2 = u_h^* M_h v_h.
On a d-simplex of Euclidean volume V_sigma, the local consistent
P1 mass matrix is
  (M_sigma)_ab = V_sigma / ((d+1)(d+2)) * (1 + delta_ab).
The global matrix is assembled by summing local contributions. It is
symmetric positive definite on the nodal space. The G3 lumped diagonal
mass remains a separate object and is not silently substituted for M_h.

The common-space extended form is
  Q_h(v) = u_h^* K_h u_h if v=J_h u_h, and +infinity otherwise,
where K_h is the P1 stiffness matrix already assembled in G0-G3.

## 3. Reproducible audit

Run:

    PYTHONPATH=. python -m pytest tests/test_flat_torus_embedding.py -q

The tests check:
- positive definiteness and symmetry of the consistent mass matrix in 2D/3D;
- exact total volume from the constant function;
- nodal interpolation and periodicity of J_h;
- affine barycentric reconstruction on every tested simplex;
- equality between the matrix inner product and the P1 L2 product.

The evaluator is an explicit pointwise realization of the P1 map. The
matrix formula supplies the exact L2 norm and inner product without
numerical quadrature.

## 4. What G4 does not establish

G4 now has a precise common-Hilbert-space embedding and a compatible
discrete Hilbert norm for the flat benchmark. It does not prove:
- uniform coercivity or compactness across refinement levels (G5);
- Mosco liminf or recovery on curved RGS meshes (G6);
- convergence of the RGS vierbein, metric, dual cells or SU(2) spin links;
- convergence of Dirac, overlap or Yang-Mills operators.

The next gate is G5: uniform geometric estimates and energy/norm controls
across a refinement sequence. Numerical condition numbers and finite-level
residuals are diagnostics; they are not a proof of uniform coercivity.
