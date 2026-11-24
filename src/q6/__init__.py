"""Q6 hypercube structure: 64 vertices, 192 edges, 240 plaquettes.

Q6 = {0,1}^6 is the 6-dimensional hypercube.
- 64 vertices (2^6)
- 192 edges (64 * 6 / 2)
- Degree 6 (each vertex connected to 6 neighbors)
- Bipartite (even/odd parity)
- Automorphism group: S_6 ⋊ C_2
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List, Dict
from itertools import product


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def q6_states() -> IntArray:
    """Generate all 64 states of Q6 = {0,1}^6.

    Convention: big-endian binary representation
    Row i corresponds to binary representation of i

    Returns:
        states: (64, 6) array of binary states
    """
    states = np.zeros((64, 6), dtype=np.int64)

    for i in range(64):
        # Convert i to 6-bit binary (big-endian)
        binary = format(i, "06b")
        states[i] = np.array([int(bit) for bit in binary], dtype=np.int64)

    return states


def q6_adjacency() -> IntArray:
    """Generate adjacency matrix of Q6.

    Two vertices are adjacent if they differ in exactly one bit.

    Returns:
        adjacency: (64, 64) symmetric adjacency matrix
    """
    states = q6_states()
    N = 64
    adjacency = np.zeros((N, N), dtype=np.int64)

    for i in range(N):
        for j in range(i + 1, N):
            # Count differing bits (Hamming distance)
            diff = np.sum(states[i] != states[j])
            if diff == 1:
                adjacency[i, j] = 1
                adjacency[j, i] = 1

    return adjacency


def q6_laplacian() -> FloatArray:
    """Generate Laplacian matrix of Q6.

    L = D - A, where D is degree matrix, A is adjacency

    Returns:
        laplacian: (64, 64) Laplacian matrix
    """
    adjacency = q6_adjacency()
    degree = np.sum(adjacency, axis=1)
    laplacian = np.diag(degree) - adjacency

    return laplacian.astype(np.float64)


def q6_edges() -> List[Tuple[int, int]]:
    """Generate list of all 192 edges.

    Returns:
        edges: List of (vertex1, vertex2) tuples
    """
    states = q6_states()
    edges = []

    for i in range(64):
        for j in range(i + 1, 64):
            diff = np.sum(states[i] != states[j])
            if diff == 1:
                edges.append((i, j))

    return edges


def q6_plaquettes() -> List[List[int]]:
    """Generate all 240 square plaquettes (faces) of Q6.

    Each plaquette is a cycle of 4 vertices.
    In 6D hypercube, there are C(6,2) * 2^4 = 15 * 16 = 240 squares.

    Returns:
        plaquettes: List of 4-vertex cycles
    """
    states = q6_states()
    plaquettes = []

    # For each pair of dimensions (choose 2 out of 6)
    for dim1 in range(6):
        for dim2 in range(dim1 + 1, 6):
            # For each combination of the other 4 dimensions
            fixed_dims = [d for d in range(6) if d not in (dim1, dim2)]

            for fixed_values in product([0, 1], repeat=4):
                # Build the 4 vertices of the square
                # Base vertex
                base = np.zeros(6, dtype=np.int64)
                for d, val in zip(fixed_dims, fixed_values):
                    base[d] = val

                # Four vertices: (0,0), (1,0), (1,1), (0,1) in (dim1, dim2)
                v1 = base.copy()
                v2 = base.copy()
                v2[dim1] = 1
                v3 = base.copy()
                v3[dim1] = 1
                v3[dim2] = 1
                v4 = base.copy()
                v4[dim2] = 1

                # Convert to vertex indices
                def to_index(state):
                    binary = "".join(str(bit) for bit in state)
                    return int(binary, 2)

                indices = [to_index(v) for v in [v1, v2, v3, v4]]
                plaquettes.append(indices)

    return plaquettes


def q6_incidence_matrices() -> Tuple[IntArray, IntArray]:
    """Generate incidence matrices B1 (vertices-edges) and B2 (edges-faces).

    B1: (64, 192) - vertex-edge incidence
    B2: (192, 240) - edge-face incidence

    Convention: B1[i,e] = ±1 if vertex i is endpoint of edge e
              B2[e,f] = ±1 if edge e is in face f

    Returns:
        B1: (64, 192) vertex-edge incidence matrix
        B2: (192, 240) edge-face incidence matrix
    """
    edges = q6_edges()
    plaquettes = q6_plaquettes()

    N_vertices = 64
    N_edges = len(edges)
    N_faces = len(plaquettes)

    # Build B1
    B1 = np.zeros((N_vertices, N_edges), dtype=np.int8)
    for e, (v1, v2) in enumerate(edges):
        B1[v1, e] = 1
        B1[v2, e] = -1  # Orientation convention

    # Build B2
    B2 = np.zeros((N_edges, N_faces), dtype=np.int8)
    edge_to_index = {edges[i]: i for i in range(len(edges))}

    for f, face in enumerate(plaquettes):
        # Face is a cycle of 4 vertices
        for i in range(4):
            v1 = face[i]
            v2 = face[(i + 1) % 4]
            edge = tuple(sorted([v1, v2]))
            if edge in edge_to_index:
                e_idx = edge_to_index[edge]
                # Orientation: +1 if going v1->v2 matches cycle direction
                # This is a simplified convention
                B2[e_idx, f] = 1 if v1 < v2 else -1

    return B1, B2


def q6_antipodal() -> IntArray:
    """Generate antipodal pairing matrix.

    For each vertex, identify its antipodal vertex
    (the vertex opposite in the hypercube, complement of all bits).

    Returns:
        antipodal: (64, 64) matrix where A[i,j] = 1 if j is antipodal to i
    """
    states = q6_states()
    N = 64
    antipodal = np.zeros((N, N), dtype=np.int64)

    for i in range(N):
        # Antipodal is bitwise complement
        complement = 1 - states[i]
        # Find index of complement
        for j in range(N):
            if np.array_equal(states[j], complement):
                antipodal[i, j] = 1
                break

    return antipodal


def verify_q6_structure(
    tol: float = 1e-10,
) -> Dict[str, any]:
    """Verify Q6 geometric and combinatorial properties.

    Args:
        tol: Tolerance for numerical checks

    Returns:
        diagnostics: Dictionary with verification results
    """
    states = q6_states()
    adjacency = q6_adjacency()
    laplacian = q6_laplacian()
    edges = q6_edges()
    plaquettes = q6_plaquettes()
    B1, B2 = q6_incidence_matrices()

    # Check vertex count
    vertex_count_correct = states.shape[0] == 64

    # Check edge count
    edge_count_correct = len(edges) == 192

    # Check face count
    face_count_correct = len(plaquettes) == 240

    # Check degree
    degrees = np.sum(adjacency, axis=1)
    degree_6 = np.all(degrees == 6)

    # Check bipartition (even/odd parity)
    parities = np.sum(states, axis=1) % 2
    # Adjacent vertices should have opposite parity
    bipartite = True
    for i, j in edges:
        if parities[i] == parities[j]:
            bipartite = False
            break

    # Check B1 @ B2 = 0
    B1B2 = B1 @ B2
    B1B2_zero = np.allclose(B1B2, 0, atol=tol)

    # Check ranks
    rank_B1 = np.linalg.matrix_rank(B1.astype(float))
    rank_B2 = np.linalg.matrix_rank(B2.astype(float))

    # Check Laplacian symmetry
    laplacian_symmetric = np.allclose(laplacian, laplacian.T, atol=tol)

    # Check non-zero counts
    nonzero_B1 = np.count_nonzero(B1)
    nonzero_B2 = np.count_nonzero(B2)

    # Check entries are in {-1, 0, 1}
    B1_entries_valid = np.all(np.isin(B1, [-1, 0, 1]))
    B2_entries_valid = np.all(np.isin(B2, [-1, 0, 1]))

    return {
        "vertex_count": int(states.shape[0]),
        "edge_count": int(len(edges)),
        "face_count": int(len(plaquettes)),
        "vertex_count_correct": vertex_count_correct,
        "edge_count_correct": edge_count_correct,
        "face_count_correct": face_count_correct,
        "degree_all_6": bool(degree_6),
        "bipartite": bipartite,
        "B1B2_zero": B1B2_zero,
        "B1B2_norm": float(np.linalg.norm(B1B2)),
        "rank_B1": int(rank_B1),
        "rank_B2": int(rank_B2),
        "expected_rank_B1": 63,
        "expected_rank_B2": 129,
        "rank_B1_correct": rank_B1 == 63,
        "rank_B2_correct": rank_B2 == 129,
        "laplacian_symmetric": laplacian_symmetric,
        "nonzero_B1": int(nonzero_B1),
        "nonzero_B2": int(nonzero_B2),
        "B1_entries_valid": B1_entries_valid,
        "B2_entries_valid": B2_entries_valid,
    }
