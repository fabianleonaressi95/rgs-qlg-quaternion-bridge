"""Icosahedron and dodecahedron geometry with golden ratio coordinates.

This module implements:
- Icosahedron vertices using golden ratio φ
- Dodecahedron vertices using golden ratio φ
- Dual relationship between icosahedron and dodecahedron
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple

from .constants import PHI, PHI_INVERSE


FloatArray = NDArray[np.float64]


def icosahedron_vertices(
    scale: float = 1.0,
    centered: bool = True,
) -> FloatArray:
    """Generate vertices of a regular icosahedron using golden ratio.

    The icosahedron has 12 vertices. Using golden ratio φ:
        - 2 vertices at (0, ±1, ±φ)
        - 2 vertices at (±1, ±φ, 0)
        - 2 vertices at (±φ, 0, ±1)
    (all sign combinations)

    Args:
        scale: Scaling factor for the vertices
        centered: If True, center the icosahedron at origin

    Returns:
        vertices: (12, 3) array of vertex coordinates
    """
    # Base coordinates using golden ratio
    coords = [
        [0, 1, PHI],
        [0, -1, PHI],
        [0, 1, -PHI],
        [0, -1, -PHI],
        [1, PHI, 0],
        [-1, PHI, 0],
        [1, -PHI, 0],
        [-1, -PHI, 0],
        [PHI, 0, 1],
        [PHI, 0, -1],
        [-PHI, 0, 1],
        [-PHI, 0, -1],
    ]

    vertices = np.array(coords, dtype=np.float64) * scale

    if centered:
        vertices = vertices - np.mean(vertices, axis=0)

    return vertices


def dodecahedron_vertices(
    scale: float = 1.0,
    centered: bool = True,
) -> FloatArray:
    """Generate vertices of a regular dodecahedron using golden ratio.

    The dodecahedron has 20 vertices. Using golden ratio φ:
        - 8 vertices at (±1, ±1, ±1) (cube)
        - 12 vertices at (0, ±1/φ, ±φ)
        - 12 vertices at (±1/φ, ±φ, 0)
        - 12 vertices at (±φ, 0, ±1/φ)
    (all sign combinations, with appropriate permutations)

    Args:
        scale: Scaling factor for the vertices
        centered: If True, center the dodecahedron at origin

    Returns:
        vertices: (20, 3) array of vertex coordinates
    """
    inv_phi = PHI_INVERSE

    # Cube vertices
    cube_vertices = []
    for x in [-1, 1]:
        for y in [-1, 1]:
            for z in [-1, 1]:
                cube_vertices.append([x, y, z])

    # Rectangle vertices on each face
    rect_vertices = [
        [0, inv_phi, PHI],
        [0, -inv_phi, PHI],
        [0, inv_phi, -PHI],
        [0, -inv_phi, -PHI],
        [inv_phi, PHI, 0],
        [-inv_phi, PHI, 0],
        [inv_phi, -PHI, 0],
        [-inv_phi, -PHI, 0],
        [PHI, 0, inv_phi],
        [PHI, 0, -inv_phi],
        [-PHI, 0, inv_phi],
        [-PHI, 0, -inv_phi],
    ]

    vertices = np.array(cube_vertices + rect_vertices, dtype=np.float64) * scale

    if centered:
        vertices = vertices - np.mean(vertices, axis=0)

    return vertices


def icosahedron_faces() -> NDArray[np.int64]:
    """Return the 20 triangular faces of the icosahedron.

    Each face is a triangle of 3 vertex indices (0-11).

    Returns:
        faces: (20, 3) array of vertex indices
    """
    # Standard icosahedron face indexing
    faces = np.array(
        [
            [0, 2, 8],
            [0, 8, 4],
            [0, 4, 6],
            [0, 6, 10],
            [0, 10, 2],
            [1, 3, 11],
            [1, 11, 5],
            [1, 5, 7],
            [1, 7, 9],
            [1, 9, 3],
            [2, 5, 8],
            [3, 2, 11],
            [4, 9, 8],
            [5, 4, 7],
            [6, 11, 10],
            [7, 6, 9],
            [8, 5, 4],
            [9, 4, 7],
            [10, 11, 6],
            [11, 3, 2],
        ],
        dtype=np.int64,
    )

    return faces


def dodecahedron_faces() -> NDArray[np.int64]:
    """Return the 12 pentagonal faces of the dodecahedron.

    Each face is a pentagon of 5 vertex indices (0-19).

    Returns:
        faces: (12, 5) array of vertex indices
    """
    # Standard dodecahedron face indexing
    faces = np.array(
        [
            [0, 1, 2, 3, 4],
            [0, 5, 6, 7, 1],
            [0, 8, 9, 10, 5],
            [1, 7, 11, 12, 2],
            [2, 12, 13, 14, 3],
            [3, 14, 15, 8, 0],
            [4, 3, 14, 16, 17],
            [4, 17, 18, 19, 9],
            [4, 9, 8, 15, 16],
            [5, 10, 18, 17, 6],
            [6, 17, 19, 11, 7],
            [10, 9, 19, 18, 13],
        ],
        dtype=np.int64,
    )

    return faces


def verify_icosahedron_geometry(
    vertices: FloatArray,
    faces: NDArray[np.int64],
    tol: float = 1e-10,
) -> dict:
    """Verify icosahedron geometric properties.

    Args:
        vertices: (12, 3) vertex coordinates
        faces: (20, 3) face indices
        tol: Tolerance for geometric checks

    Returns:
        diagnostics: Dictionary with verification results
    """
    N_vertices = vertices.shape[0]
    N_faces = faces.shape[0]

    # Check edge lengths (should all be equal)
    edge_lengths = []
    for face in faces:
        for i in range(3):
            v1 = vertices[face[i]]
            v2 = vertices[face[(i + 1) % 3]]
            edge_lengths.append(np.linalg.norm(v2 - v1))

    edge_lengths = np.array(edge_lengths)
    edge_std = np.std(edge_lengths)
    edge_mean = np.mean(edge_lengths)

    # Check face areas (should all be equal)
    face_areas = []
    for face in faces:
        v0 = vertices[face[0]]
        v1 = vertices[face[1]]
        v2 = vertices[face[2]]
        cross = np.cross(v1 - v0, v2 - v0)
        area = 0.5 * np.linalg.norm(cross)
        face_areas.append(area)

    face_areas = np.array(face_areas)
    area_std = np.std(face_areas)

    return {
        "num_vertices": int(N_vertices),
        "num_faces": int(N_faces),
        "expected_vertices": 12,
        "expected_faces": 20,
        "vertex_count_correct": N_vertices == 12,
        "face_count_correct": N_faces == 20,
        "edge_length_mean": float(edge_mean),
        "edge_length_std": float(edge_std),
        "edges_equal": edge_std < tol,
        "face_area_std": float(area_std),
        "faces_equal": area_std < tol,
    }
