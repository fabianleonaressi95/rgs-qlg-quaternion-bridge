import numpy as np
import pandas as pd
from scipy.linalg import solve

def run_implicit_newmark():
    # Load geometry
    nodes = pd.read_csv("/workspace/qc_nodes.csv", header=None).values
    edges = pd.read_csv("/workspace/qc_edges.csv", header=None).values
    num_nodes = len(nodes)
    
    # Material / Simulation Parameters
    mass_per_node = 0.01
    k = 1e3
    M = np.eye(3 * num_nodes) * mass_per_node
    
    # Newmark-Beta constants (Average acceleration method)
    gamma = 0.5
    beta = 0.25
    dt = 1e-3
    steps = 50
    
    # State vectors
    u = np.zeros(3 * num_nodes)
    v = np.zeros(3 * num_nodes)
    a = np.zeros(3 * num_nodes)
    
    # Initial impact
    z_max = np.max(nodes[:, 2])
    impact_mask = nodes[:, 2] > (z_max - 0.2)
    v[impact_mask.repeat(3)] = -5.0
    
    # Precompute structural Stiffness Matrix (Simplified)
    K = np.zeros((3 * num_nodes, 3 * num_nodes))
    for edge in edges:
        i, j = edge
        # Local stiffness matrix components (simplified 1D spring in 3D)
        dir_vec = (nodes[i] - nodes[j])
        dist = np.linalg.norm(dir_vec)
        n = dir_vec / dist
        k_mat = k * np.outer(n, n)
        
        # Add to global K
        for dim in range(3):
            K[3*i+dim, 3*i+dim] += k_mat[dim, dim]
            K[3*i+dim, 3*j+dim] -= k_mat[dim, dim]
            K[3*j+dim, 3*i+dim] -= k_mat[dim, dim]
            K[3*j+dim, 3*j+dim] += k_mat[dim, dim]

    # Effective Stiffness
    K_eff = K + (1 / (beta * dt**2)) * M
    
    for _ in range(steps):
        # Predicted displacement/velocity
        u_pred = u + dt * v + (0.5 - beta) * dt**2 * a
        v_pred = v + (1 - gamma) * dt * a
        
        # RHS (Load vector)
        F = -K @ u_pred
        
        # Solve for new acceleration
        a = solve(K_eff, F)
        
        # Update
        u = u_pred + beta * dt**2 * a
        v = v_pred + gamma * dt * a
        
    print("Simulation Complete (Implicit Integration).")
    print(f"Final System Energy (Approx): {0.5 * np.dot(v, M @ v) + 0.5 * np.dot(u, K @ u):.2f} J")

if __name__ == "__main__":
    run_implicit_newmark()
