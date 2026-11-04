# Equivariant icosahedral-axis / Q6-direction intertwiner

## Construction

The twelve normalized icosahedral vertices are reduced to six antipodal axes
using the deterministic ordering of `build_icosahedral_c6()`. The 60 proper
rotations are reconstructed from the vertex set: three linearly independent
vertex images determine each rotation, and candidates are retained only if
they belong to SO(3) and preserve all twelve vertices.

For each rotation R and representative axis a_j, find the unique target a_k
and sign s such that R a_j = s a_k. This gives a signed permutation S(R):
\[
S(R)\in O(6),\qquad S(R_1R_2)=S(R_1)S(R_2).
\]

## Q6 realization

Use the big-endian Q6 convention in `src/q6_complex.py`. If
S(R)_{kj}=+1, copy bit j to bit k; if it is -1, complement bit j before
placing it at k:
\[
x'_k=x_j\quad(+1),\qquad x'_k=1-x_j\quad(-1).
\]
This is a Q6 graph automorphism: it permutes all 64 vertices and sends every
one-bit edge to a one-bit edge. The sign is retained as bit complementation.

## Intertwining statement

The signed six-axis representation of the proper icosahedral group is
realized on the six Q6 coordinate directions. The C6 tensor is equivariant:
\[
S(R)^T C_6 S(R)=C_6
\]
for all 60 rotations. The code also checks the group homomorphism law and
all 192 edge directions at all 64 vertices.

## Evidence and limits

- **PROVED from the construction:** proper rotations induce signed
  permutations; the bit action is a Q6 automorphism.
- **NUMERICALLY VERIFIED:** order 60, representation law, C6 covariance,
  and Q6 edge preservation.
- **NOT CLAIMED:** a metric isometry between the icosahedron and Q6, a unique
  physical frame identification, spacetime geometry, or fermionic dynamics.
  The initial axis-to-bit labeling is a declared coordinate convention.
- **OPEN:** lift the intertwiner to Q6 cochains with connection transport and
  test a Clifford-valued Dirac operator, gauge covariance and continuum limit.

Implementation: `src/rgs_icosahedral_q6_intertwiner.py`; tests:
`tests/test_rgs_icosahedral_q6_intertwiner.py`.
