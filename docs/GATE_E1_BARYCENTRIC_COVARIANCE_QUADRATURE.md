# Gate E.1 — Barycentric–covariance quadrature on a planar triangle

**Status: PROVED (local cell theorem); NUMERICALLY VERIFIED (implementation tests).**

This gate records a local Euclidean statement only. It does not assert that the
Q6 2-skeleton is a closed embedded surface, that the recovered RGS coordinates
define a refinement family, or that the Q6 vertex covariance is a triangle
covariance.

## 1. Geometric data

Let F be a nondegenerate planar triangle in R^3 with vertices v_0,v_1,v_2,
area A_F and diameter h_F. Its areal barycentre is

    c_F = (v_0 + v_1 + v_2)/3.

Its normalized ambient covariance is

    Sigma_F = (1/A_F) integral_F (x-c_F)(x-c_F)^T dA.

Uniform barycentric coordinates satisfy E[lambda_i^2]=1/6 and
E[lambda_i lambda_j]=1/12 for i != j. Consequently

    Sigma_F = (1/12) sum_{i=0}^2 (v_i-c_F)(v_i-c_F)^T.

The covariance is a 3x3 positive-semidefinite matrix of rank at most two.

## 2. Theorem E.1

For g in C^3(U), where U is an open Euclidean neighbourhood of F, Taylor
formula at c_F gives, with y=x-c_F,

    g(x) = g(c_F) + Dg(c_F)[y] + (1/2)D^2g(c_F)[y,y] + r_3(x),

and

    |r_3(x)| <= (1/6) sup_F ||D^3 g|| ||y||^3.

The linear term integrates to zero by the barycentre property. Therefore

    integral_F g dA = A_F g(c_F) + (A_F/2) D^2g(c_F):Sigma_F + R_F,

with

    |R_F| <= (A_F h_F^3 / 6) sup_F ||D^3g||.

Here the third derivative norm is the operator norm of the trilinear form.
The quadrature error for Q_F(g)=A_F g(c_F) has the sign

    Q_F(g)-integral_F g dA = -(A_F/2)D^2g(c_F):Sigma_F - R_F.

For every polynomial of degree at most two, R_F=0. In particular, for

    g(x)=a+b^T x+(1/2)x^T H x,

with symmetric H,

    integral_F g dA = A_F [g(c_F) + (1/2) H:Sigma_F].

## 3. Exact implementation

The module src/barycentric_triangle_quadrature.py provides:
- triangle_barycentric_data(vertices): area, centroid and covariance;
- quadratic_triangle_integral(...): exact quadratic integral;
- triangle_cubic_remainder_bound(...): A h^3 sup||D^3g||/6.

The implementation rejects malformed, non-finite and degenerate input.
The dedicated test file covers centroid, covariance, constant/affine exactness,
quadratic exactness against an independent degree-three Dunavant rule, the
cubic remainder bound and h^5 scaling under similarity refinement, and
degenerate-input rejection.

Run:

    PYTHONPATH=. python -m pytest tests/test_barycentric_triangle_quadrature.py -q

## 4. Relation to existing Q6 results

The existing Q6 FEM record concerns the 64 hypercube vertices and 192 graph
edges. Its vertex barycentre (1/2,...,1/2) and covariance (1/4)I_6 are exact
finite combinatorial facts. The coordinate-field edge-gradient and
Dirichlet-energy tests are a separate discrete FEM result.

Neither the 6D covariance of the Q6 vertex measure nor the Q6 edge-gradient
identity is identified here with Sigma_F or with a bound on third derivatives.
Such a bridge requires a specified geometric map and a coherent family of
surface cells. That is Gate E.2, not part of this theorem.

## 5. Surface-only functions and limitations

The Hessian in this theorem is the ambient Euclidean Hessian of a function
defined on a neighbourhood of the planar cell. If a function is defined only
on a curved limiting surface, an ambient extension or an intrinsic
tangential/curvature formulation must first be specified.

The theorem is local to a planar triangle. It does not prove global quadrature
convergence, shape regularity of a refinement family, surface closure,
orientability, absence of self-intersections, Mosco convergence, or a
Bernoulli-number expansion. The Q6 square complex and the Gate D integral
chain filling do not by themselves supply an embedded closed surface.

## 6. Evidence ledger

- **PROVED:** barycentric covariance formula and local Taylor identity/bound.
- **PROVED:** exactness for all polynomials of degree at most two.
- **NUMERICALLY VERIFIED:** implementation checks, conditional on running the
  test suite in the target environment.
- **OPEN:** compatibility with RGS maps, refined surface geometry and any
  Euler–Maclaurin/Bernoulli asymptotics.