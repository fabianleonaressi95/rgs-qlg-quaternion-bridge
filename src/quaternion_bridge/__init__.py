"""Quaternion bridge: connection from vierbein to SU(2) and 2I group.

This module implements:
- Quaternionic connection M_ij = E_j * E_i^{-1}
- Polar decomposition M_ij = R_ij * S_ij
- SU(2) lift of rotations
- Distance from binary icosahedral group 2I
- Holonomy calculation on Q6 plaquettes
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List, Dict
from scipy.linalg import polar

try:
    from ..quaternions import Quaternion
    from ..rotations import quaternion_to_matrix, axis_angle_to_quaternion
    from ..binary_icosahedral import binary_icosahedral_elements
    from ..q6 import q6_edges, q6_plaquettes
except ImportError:
    from quaternions import Quaternion
    from rotations import quaternion_to_matrix, axis_angle_to_quaternion
    from binary_icosahedral import binary_icosahedral_elements
    from q6 import q6_edges, q6_plaquettes


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


def quaternionic_connection(
    E: FloatArray,
    adjacency: List[Tuple[int, int]],
) -> Dict[Tuple[int, int], FloatArray]:
    """Compute quaternionic connection M_ij = E_j * E_i^{-1}.

    Args:
        E: (N, 3, 3) array of vierbein matrices
        adjacency: List of (i, j) edge pairs

    Returns:
        connection: Dictionary mapping (i, j) to M_ij (3, 3) matrix
    """
    connection = {}

    for i, j in adjacency:
        # Compute M_ij = E_j * E_i^{-1}
        E_i = E[i]
        E_j = E[j]

        # Inverse of E_i
        E_i_inv = np.linalg.inv(E_i)

        # Connection matrix
        M_ij = E_j @ E_i_inv
        connection[(i, j)] = M_ij

    return connection


def polar_decomposition_connection(
    connection: Dict[Tuple[int, int], FloatArray],
) -> Dict[Tuple[int, int], Tuple[FloatArray, FloatArray]]:
    """Compute polar decomposition M_ij = R_ij * S_ij.

    Args:
        connection: Dictionary of M_ij matrices

    Returns:
        decomposed: Dictionary mapping (i, j) to (R_ij, S_ij)
            where R_ij ∈ SO(3) and S_ij > 0
    """
    decomposed = {}

    for edge, M in connection.items():
        # Polar decomposition: M = R * S
        # R is orthogonal, S is symmetric positive semidefinite
        R, S = polar(M, side="right")

        # Ensure R ∈ SO(3) (det = +1)
        if np.linalg.det(R) < 0:
            # Flip sign if det = -1
            R = -R

        decomposed[edge] = (R, S)

    return decomposed


def so3_to_su2_quaternion(
    R: FloatArray,
) -> Quaternion:
    """Lift SO(3) rotation matrix to SU(2) quaternion.

    Args:
        R: (3, 3) rotation matrix

    Returns:
        q: Unit quaternion representing the rotation
    """
    # Extract axis and angle from rotation matrix
    trace = np.trace(R)
    angle = np.arccos((trace - 1) / 2)

    if np.abs(angle) < 1e-10:
        # Identity rotation
        return Quaternion.identity()

    # Axis of rotation
    axis = np.array([
        R[2, 1] - R[1, 2],
        R[0, 2] - R[2, 0],
        R[1, 0] - R[0, 1],
    ])
    axis = axis / (2 * np.sin(angle))

    # Create quaternion
    q = axis_angle_to_quaternion(axis, angle)
    return q


def su2_connection_from_vierbein(
    E: FloatArray,
    adjacency: List[Tuple[int, int]],
) -> Dict[Tuple[int, int], Quaternion]:
    """Compute SU(2) connection from vierbein.

    Pipeline:
        1. M_ij = E_j * E_i^{-1}
        2. Polar decomposition: M_ij = R_ij * S_ij
        3. Lift R_ij to quaternion q_ij ∈ SU(2)

    Args:
        E: (N, 3, 3) vierbein matrices
        adjacency: List of (i, j) edges

    Returns:
        su2_connection: Dictionary mapping (i, j) to quaternion q_ij
    """
    # Step 1: Quaternionic connection
    connection = quaternionic_connection(E, adjacency)

    # Step 2: Polar decomposition
    decomposed = polar_decomposition_connection(connection)

    # Step 3: SU(2) lift
    su2_connection = {}

    for edge, (R, S) in decomposed.items():
        q = so3_to_su2_quaternion(R)
        su2_connection[edge] = q

    return su2_connection


def distance_from_2I(
    q: Quaternion,
    elements_2I: FloatArray,
) -> float:
    """Compute distance from quaternion to binary icosahedral group 2I.

    d_2I(q) = min_{g∈2I} min(||q - g||, ||q + g||)

    The min over q+g accounts for the double cover SU(2) → SO(3).

    Args:
        q: Quaternion to test
        elements_2I: (120, 4) array of 2I elements (scalar-first)

    Returns:
        distance: Minimum distance to 2I
    """
    q_arr = np.array(q.as_tuple(), dtype=np.float64)

    min_distance = np.inf

    for g in elements_2I:
        # Distance to g
        dist_g = np.linalg.norm(q_arr - g)
        # Distance to -g (same SO(3) rotation)
        dist_minus_g = np.linalg.norm(q_arr + g)

        min_distance = min(min_distance, dist_g, dist_minus_g)

    return float(min_distance)


def analyze_2I_proximity(
    su2_connection: Dict[Tuple[int, int], Quaternion],
    elements_2I: FloatArray = None,
) -> Dict[str, any]:
    """Analyze proximity of SU(2) connection to 2I group.

    Args:
        su2_connection: Dictionary of quaternion link variables
        elements_2I: (120, 4) array of 2I elements (computed if None)

    Returns:
        analysis: Dictionary with proximity statistics
    """
    if elements_2I is None:
        elements_2I = binary_icosahedral_elements()

    distances = []

    for edge, q in su2_connection.items():
        dist = distance_from_2I(q, elements_2I)
        distances.append(dist)

    distances = np.array(distances)

    return {
        "num_links": len(distances),
        "mean_distance": float(np.mean(distances)),
        "std_distance": float(np.std(distances)),
        "min_distance": float(np.min(distances)),
        "max_distance": float(np.max(distances)),
        "median_distance": float(np.median(distances)),
        "fraction_close": float(np.mean(distances < 0.1)),  # Arbitrary threshold
        "all_distances": distances.tolist(),
    }


def compute_plaquette_holonomy(
    su2_connection: Dict[Tuple[int, int], Quaternion],
    plaquette: List[int],
    orientation: List[int] = None,
) -> Quaternion:
    """Compute holonomy around a plaquette.

    U_f = ∏_{e∈∂f} q_e^{ε_fe}

    Args:
        su2_connection: Dictionary of quaternion link variables
        plaquette: List of 4 vertex indices
        orientation: List of ±1 for each edge (default: assume cyclic)

    Returns:
        holonomy: Quaternion representing the holonomy
    """
    if orientation is None:
        # Assume cyclic orientation: 0→1→2→3→0
        orientation = [1, 1, 1, 1]

    # Start with identity
    holonomy = Quaternion.identity()

    for i in range(4):
        v1 = plaquette[i]
        v2 = plaquette[(i + 1) % 4]
        epsilon = orientation[i]

        # Get link variable
        if (v1, v2) in su2_connection:
            q = su2_connection[(v1, v2)]
        elif (v2, v1) in su2_connection:
            # Reverse direction: use inverse
            q = su2_connection[(v2, v1)].inverse()
        else:
            raise ValueError(f"Link ({v1}, {v2}) not in connection")

        # Apply with orientation
        if epsilon > 0:
            holonomy = q * holonomy
        else:
            holonomy = q.inverse() * holonomy

    return holonomy


def compute_all_holonomies(
    su2_connection: Dict[Tuple[int, int], Quaternion],
    plaquettes: List[List[int]] = None,
) -> List[Quaternion]:
    """Compute holonomies for all 240 Q6 plaquettes.

    Args:
        su2_connection: Dictionary of quaternion link variables
        plaquettes: List of plaquettes (computed from Q6 if None)

    Returns:
        holonomies: List of 240 holonomy quaternions
    """
    if plaquettes is None:
        plaquettes = q6_plaquettes()

    holonomies = []

    for plaquette in plaquettes:
        h = compute_plaquette_holonomy(su2_connection, plaquette)
        holonomies.append(h)

    return holonomies


def analyze_holonomy_distribution(
    holonomies: List[Quaternion],
    elements_2I: FloatArray = None,
) -> Dict[str, any]:
    """Analyze distribution of plaquette holonomies.

    Args:
        holonomies: List of holonomy quaternions
        elements_2I: (120, 4) array of 2I elements

    Returns:
        analysis: Dictionary with statistical analysis
    """
    if elements_2I is None:
        elements_2I = binary_icosahedral_elements()

    # Compute curvature angles
    try:
        from ..rotations import rotation_distance
    except ImportError:
        from rotations import rotation_distance

    angles = []
    distances_2I = []

    for h in holonomies:
        # Curvature angle
        angle = rotation_distance(h, Quaternion.identity())
        angles.append(angle)

        # Distance from 2I
        dist = distance_from_2I(h, elements_2I)
        distances_2I.append(dist)

    angles = np.array(angles)
    distances_2I = np.array(distances_2I)

    return {
        "num_plaquettes": len(holonomies),
        "curvature_angles": {
            "mean": float(np.mean(angles)),
            "std": float(np.std(angles)),
            "min": float(np.min(angles)),
            "max": float(np.max(angles)),
            "median": float(np.median(angles)),
        },
        "2I_distances": {
            "mean": float(np.mean(distances_2I)),
            "std": float(np.std(distances_2I)),
            "min": float(np.min(distances_2I)),
            "max": float(np.max(distances_2I)),
            "median": float(np.median(distances_2I)),
        },
        "flat_plaquettes": int(np.sum(angles < 1e-6)),
        "max_curvature": float(np.max(angles)),
    }
