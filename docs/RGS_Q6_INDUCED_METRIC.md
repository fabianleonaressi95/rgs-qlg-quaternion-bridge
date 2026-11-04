# RGS–Q6 induced metric and graph P1 audit

## Purpose and provenance

This phase takes the 64-point RGS seed from
`load_rgs_seed_data()`. That routine reconstructs the coordinates and checks
the archived SHA-256 hashes for the 64x3 seed, Q6 edge list and normalized
four-component edge coframe. The canonical Q6 incidence ordering is obtained
from `build_q6_complex()`; the archived and canonical edge endpoint sets are
compared explicitly, not assumed to have identical enumeration order.

The map under study is
[
X:{0,1}^6longrightarrowmathbb R^3,qquad vmapsto X_v.
]
The six-bit vertex order is big-endian. The six finite differences at a
vertex define (J_vinmathbb R^{3	imes6}), and
[
G_v=J_v^{mathsf T}J_v.
]
Since (operatorname{rank}J_vle3), (G_v) is necessarily singular as a
6x6 matrix. This is an exact rank bound, not a numerical defect. The
construction is a degenerate pullback semimetric on the six-dimensional
parameter directions; it is not a nondegenerate 6D Riemannian metric.

## Face geometry

For each canonical square ((v,v_a,v_{ab},v_b)), define
[
a=X_{v_a}-X_v,quad b=X_{v_b}-X_v,quad
G_f=egin{pmatrix}acdot a&acdot b\acdot b&bcdot bend{pmatrix}.
]
The parallelogram area diagnostic is (sqrt{det G_f}). The actual
piecewise-linear quadrilateral area is computed as the sum of the areas of
triangles ((v,v_a,v_{ab})) and ((v,v_{ab},v_b)). The closure defect
(|X_{v_{ab}}-X_{v_a}-X_{v_b}+X_v|) measures non-affinity of the
image quadrilateral. It is not silently set to zero.

## Graph P1 finite elements

Each Q6 edge is treated as a one-dimensional affine P1 element of its
measured Euclidean length (L_e=|X_v-X_u|):
[
M_e=rac{L_e}{6}egin{pmatrix}2&1\1&2end{pmatrix},qquad
K_e=rac1{L_e}egin{pmatrix}1&-1\-1&1end{pmatrix}.
]
These matrices are assembled into 64x64 vertex mass and stiffness matrices.
The mass is positive definite for the connected Q6 graph; stiffness is
positive semidefinite and annihilates constants. This is a metric-weighted
graph P1 FEM, not a six-dimensional volume FEM and not a Whitney 1-form
Hodge star.

## Symmetry gate

The 48 signed coordinate permutations are checked on the recovered point
cloud. Their induced permutations are separately tested against the Q6 edge
and square sets. Coordinate (O_h) symmetry does not imply Q6-label
equivariance. The audit records the number of group elements preserving
edges and faces; no symmetry is imposed by hand, and Q6 equivariance is a
diagnostic rather than a precondition for the graph metric construction.

## Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_induced_metric.py -q
PYTHONPATH=. python -m src.rgs_q6_induced_metric --output results/rgs_q6_induced_metric
```

The JSON report and compressed NPZ contain all finite arrays, residual
diagnostics and SHA-256 of the exported archive.

## Evidence ledger

- **NUMERICALLY VERIFIED:** seed reconstruction hashes, finite edge/face
  geometry, graph P1 matrix properties and the separate symmetry gates.
- **ANALYTICALLY FORCED:** the 6x6 pullback Gram has rank at most 3 because
  the target is (mathbb R^3).
- **OPEN:** a nondegenerate six-dimensional metric, a curved RGS Whitney
  Hodge star, Mosco convergence, continuum Laplace–Beltrami/Yang–Mills,
  TEGR closure and physical interpretation.
