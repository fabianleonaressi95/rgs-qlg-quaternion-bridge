import numpy as np

def generate_icosahedral_lattice(num_points=100, projection_width=1.0):
    """
    Generates a set of points for an icosahedral quasicrystalline lattice
    using the cut-and-project method from 6D to 3D.
    
    1. 6D hypercubic grid points.
    2. Project into 3D using the golden ratio-based transformation matrix.
    3. Apply the 'cut' filter based on the projection window.
    """
    # Golden ratio
    phi = (1 + np.sqrt(5)) / 2
    
    # Basis vectors for 3D projection from 6D
    # These represent the projection of 6 orthogonal axes in 6D to 3D
    # optimized for icosahedral symmetry
    v = np.array([
        [1, 0, 1/phi],
        [-1, 0, 1/phi],
        [0, 1/phi, 1],
        [0, -1/phi, 1],
        [1/phi, 1, 0],
        [-1/phi, 1, 0]
    ])
    v = v / np.linalg.norm(v, axis=1)[:, np.newaxis]
    
    # Generate a smaller sample for manageable FEA
    grid_range = 1  # Reduced from 3 to limit point count
    points_6d = np.array(np.meshgrid(*[np.arange(-grid_range, grid_range+1)]*6)).T.reshape(-1, 6)
    
    # Project to 3D
    points_3d = np.dot(points_6d, v)
    
    # Cut: Select points that project into the 'physical' space
    # (Simplified window: points with low norm in the perpendicular space)
    # The perpendicular space basis is the orthogonal complement of the 3D projection basis
    # Here we simplify by using a spherical window in 3D
    mask = np.linalg.norm(points_3d, axis=1) < (grid_range * projection_width)
    return points_3d[mask]

def create_truss_from_points(points, threshold=1.5):
    """
    Connect points to form a truss network (edges).
    Simple distance-based connection.
    """
    from scipy.spatial import KDTree
    tree = KDTree(points)
    pairs = tree.query_pairs(r=threshold)
    return list(pairs)

if __name__ == "__main__":
    # Generate geometry
    qc_points = generate_icosahedral_lattice(projection_width=1.5)
    truss_edges = create_truss_from_points(qc_points)
    
    print(f"Generated {len(qc_points)} nodes and {len(truss_edges)} edges.")
    
    # Save for visualization/FEA
    np.savetxt("/workspace/qc_nodes.csv", qc_points, delimiter=",")
    with open("/workspace/qc_edges.csv", "w") as f:
        for edge in truss_edges:
            f.write(f"{edge[0]},{edge[1]}\n")
    print("Geometry saved to /workspace/qc_nodes.csv and /workspace/qc_edges.csv")
