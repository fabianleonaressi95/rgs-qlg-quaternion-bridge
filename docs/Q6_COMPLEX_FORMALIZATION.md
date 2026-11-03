# Q6 canonical square complex: formal construction and audit

## 1. Object and conventions

Let V = {0,1}^6, identified with integer labels 0,...,63.
The repository convention is big-endian: axis a in {0,...,5} flips
the mask m_a = 2^(5-a). Vertices are joined when they differ in exactly
one bit. Each edge is oriented from its smaller integer label to its larger.

For every pair of axes a < b, and every base vertex whose bits on those
two axes vanish, include the square
(v, v xor m_a, v xor m_a xor m_b, v xor m_b).
Its orientation is the displayed cyclic order.

The counts follow directly:
- |V| = 2^6 = 64;
- |E| = 6 * 2^5 = 192;
- |F2| = binomial(6,2) * 2^4 = 240.

This is the canonical 2-skeleton of the six-cube. The combinatorial
construction is exact and does not require geometric coordinates.

## 2. Incidence operators and chain identity

Let B1 be the 64 x 192 oriented vertex-edge incidence matrix and B2
the 192 x 240 oriented edge-face incidence matrix. The sign in each
column of B2 follows the chosen face orientation and canonical edge
orientation.

Every square boundary has four edges and its vertex boundary cancels
pairwise. Therefore B1 B2 = 0. The implementation checks this as an
exact integer identity, not a floating-point tolerance test.

## 3. Graph Laplacian and exact spectrum

The 0-form graph Laplacian is L0 = B1 B1^T = 6 I - A(Q6).
The six-cube is the Cartesian product of six copies of K2. The eigenvalues
of the K2 Laplacian are {0,2}; product eigenvalues are sums of six such
terms. Hence:
- eigenvalue 2j;
- multiplicity binomial(6,j);
- j = 0,...,6.

The multiplicities are 1, 6, 15, 20, 15, 6, 1 and sum to 64.
The code diagonalizes the matrix numerically and clusters values by
nearest integer before checking multiplicities. This avoids treating
roundoff-level eigenvalue splitting as distinct eigenvalues.

## 4. Graph cycle rank versus cellular homology

These invariants must not share the same label.

The cycle rank (first Betti number of the 1-dimensional graph) is
|E| - |V| + 1 = 192 - 64 + 1 = 129.

After attaching all 240 square 2-cells, the cellular chain complex is
C2 --B2--> C1 --B1--> C0. The graph is connected, so rank(B1)=63.
The chain identity gives rank_Q(B2) <= 192 - 63 = 129.

The implementation computes rank(B2) modulo the prime 1,000,003.
A rank of 129 over this finite field is a lower bound for rational rank;
combined with the upper bound, it certifies rank_Q(B2)=129. Consequently
the square 2-skeleton has cellular Betti numbers:
- b0 = 1;
- b1 = 0;
- b2 = 240 - 129 = 111.

Its Euler characteristic is 64 - 192 + 240 = 112 = 1 - 0 + 111.

Thus 129 is the graph cycle rank, not the first Betti number after the
square faces have been attached.

## 5. Antipodal complement

On six-bit labels define iota(v) = v xor 63. Then iota^2 = id, and bit
flips commute with complementation, so iota preserves vertices, edges
and square faces. This is an exact combinatorial antipodal involution
of the abstract cube.

It is not automatically the same as a coordinate antipodal map on a
separately reconstructed point cloud X_Oh. The latter requires an
explicit vertex correspondence and a coordinate residual audit.

## 6. Evidence classification

### PROVED (finite combinatorics)
- Counts 64, 192 and 240 follow directly from the cube construction.
- B1 B2 = 0 follows from the oriented boundary-of-boundary identity.
- The Laplacian spectrum and multiplicities follow from the Cartesian
  product spectrum of K2^6.
- The graph cycle rank is 129.
- Given the verified modular rank and chain-identity upper bound,
  the cellular Betti numbers of this 2-skeleton are (1,0,111).

### NUMERICALLY VERIFIED BY THE AUDIT
- Matrix construction and exact integer residual B1 B2 = 0.
- Floating-point eigenspectrum against the exact multiplicity table.
- Finite-field rank witness for B2.

### OPEN
- Identification of this abstract Q6 labeling with the actual 64-point
  RGS X_Oh embedding, including the induced action of all 48 signed
  coordinate permutations.
- Whether the geometric nearest-neighbor graph equals this canonical Q6
  edge set.
- SU(2) lattice gauge dynamics on these 192 edges and 240 faces.
- Any continuum limit, Green-Kubo coefficient or physical interpretation.

## 7. Reproduction

From the repository root run:
- python -m pytest tests/test_q6_complex.py tests/test_q6_spectral_audit.py -q
- python examples/q6_formalization_audit.py

The example writes a JSON report and NumPy arrays to the selected output
directory. It does not load or fabricate the RGS coordinate file.
