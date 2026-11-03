"""Star tetrahedron (stellated octahedron) geometry.

The star tetrahedron is formed by two interpenetrating tetrahedra,
creating an octahedron at their intersection. This structure is
fundamental in sacred geometry and appears in the RGS golden geometry.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List

from .constants import PHI


FloatArray = NDArray[np.float64]


def tetrahedron_vertices(
    scale: float = 1.0,
    centered: bool = True,
) -> FloatArray:
    """Generate vertices of a regular tetrahedron.

    The tetrahedron has 4 vertices. Using coordinates:
        (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)

    Args:
        scale: Scaling factor
        centered: If True, center at origin

    Returns:
        vertices: (4, 3) array of vertex coordinates
    """
    coords = [
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
    ]

    vertices = np.array(coords, dtype=np.float64) * scale

    if centered:
        vertices = vertices - np.mean(vertices, axis=0)

    return vertices


def star_tetrahedron_vertices(
    scale: float = 1.0,
    centered: bool = True,
) -> Tuple[FloatArray, FloatArray]:
    """Generate vertices of a star tetrahedron (two interpenetrating tetrahedra).

    The star tetrahedron consists of:
    - Tetrahedron 1 (pointing up)
    - Tetrahedron 2 (pointing down, inverted)

    Together they form 8 vertices that create an octahedron at their center.

    Args:
        scale: Scaling factor
        centered: If True, center at origin

    Returns:
        tet1_vertices: (4, 3) vertices of first tetrahedron
        tet2_vertices: (4, 3) vertices of second (inverted) tetrahedron
    """
    # First tetrahedron
    tet1 = tetrahedron_vertices(scale, centered=False)

    # Second tetrahedron (inverted)
    tet2 = -tetrahedron_vertices(scale, centered=False)

    if centered:
        tet1 = tet1 - np.mean(tet1, axis=0)
        tet2 = tet2 - np.mean(tet2, axis=0)

    return tet1, tet2


def star_tetrahedron_combined(
    scale: float = 1.0,
    centered: bool = True,
) -> FloatArray:
    """Return all 8 vertices of the star tetrahedron.

    Args:
        scale: Scaling factor
        centered: If True, center at origin

    Returns:
        vertices: (8, 3) array of all vertex coordinates
    """
    tet1, tet2 = star_tetrahedron_vertices(scale, centered)
    return np.vstack([tet1, tet2])


def star_tetrahedron_edges() -> List[Tuple[int, int]]:
    """Return the edge connectivity of the star tetrahedron.

    Returns:
        edges: List of (vertex1, vertex2) tuples
    """
    # Tetrahedron 1 edges (vertices 0-3)
    tet1_edges = [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    ]

    # Tetrahedron 2 edges (vertices 4-7)
    tet2_edges = [
        (4, 5),
        (4, 6),
        (4, 7),
        (5, 6),
        (5, 7),
        (6, 7),
    ]

    # Inter-tetrahedron edges (intersection octahedron)
    # These connect the two tetrahedra through their intersection
    inter_edges = [
        (0, 5),
        (0, 6),
        (0, 7),
        (1, 4),
        (1, 6),
        (1, 7),
        (2, 4),
        (2, 5),
        (2, 7),
        (3, 4),
        (3, 5),
        (3, 6),
    ]

    return tet1_edges + tet2_edges + inter_edges


def verify_star_tetrahedron(
    vertices: FloatArray,
    tol: float = 1e-10,
) -> dict:
    """Verify star tetrahedron geometric properties.

    Args:
        vertices: (8, 3) vertex coordinates
        tol: Tolerance for checks

    Returns:
        diagnostics: Dictionary with verification results
    """
    N = vertices.shape[0]

    # Check if we have 8 vertices
    vertex_count_correct = N == 8

    # Check centroid
    centroid = np.mean(vertices, axis=0)
    centered = np.allclose(centroid, np.zeros(3), atol=tol)

    # Check symmetry: vertices should come in antipodal pairs
    antipodal_pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            if np.allclose(vertices[i] + vertices[j], np.zeros(3), atol=tol):
                antipodal_pairs.append((i, j))

    # Should have 4 antipodal pairs
    correct_antipodal_structure = len(antipodal_pairs) == 4

    return {
        "num_vertices": int(N),
        "expected_vertices": 8,
        "vertex_count_correct": vertex_count_correct,
        "centered": centered,
        "centroid": centroid.tolist(),
        "antipodal_pairs": antipodal_pairs,
        "antipodal_structure_correct": correct_antipodal_structure,
    }
