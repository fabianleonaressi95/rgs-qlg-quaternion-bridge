import numpy as np
import pandas as pd
from scipy.linalg import solve

def generate_cubic_truss(size=5):
    """Generates a simple cubic lattice."""
    x = np.linspace(0, size, size)
    y = np.linspace(0, size, size)
    z = np.linspace(0, size, size)
    X, Y, Z = np.meshgrid(x, y, z)
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    
    # Simple connectivity
    edges = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if np.linalg.norm(nodes[i] - nodes[j]) <= 1.0:
                edges.append([i, j])
    return nodes, np.array(edges)

def run_implicit_newmark_periodic():
    nodes, edges = generate_cubic_truss(size=5)
    num_nodes = len(nodes)
    
    # Same material/simulation parameters for comparison
    mass_per_node = 0.01
    k = 1e3
    M = np.eye(3 * num_nodes) * mass_per_node
    gamma, beta, dt, steps = 0.5, 0.25, 1e-3, 50
    
    u, v, a = np.zeros(3 * num_nodes), np.zeros(3 * num_nodes), np.zeros(3 * num_nodes)
    
    # Impact
    z_max = np.max(nodes[:, 2])
    impact_mask = nodes[:, 2] > (z_max - 0.2)
    v[impact_mask.repeat(3)] = -5.0
    
    # Stiffness Matrix
    K = np.zeros((3 * num_nodes, 3 * num_nodes))
    for edge in edges:
        i, j = edge
        n = (nodes[i] - nodes[j]) / np.linalg.norm(nodes[i] - nodes[j])
        k_mat = k * np.outer(n, n)
        for dim in range(3):
            K[3*i+dim, 3*i+dim] += k_mat[dim, dim]
            K[3*i+dim, 3*j+dim] -= k_mat[dim, dim]
            K[3*j+dim, 3*i+dim] -= k_mat[dim, dim]
            K[3*j+dim, 3*j+dim] += k_mat[dim, dim]

    K_eff = K + (1 / (beta * dt**2)) * M
    
    for _ in range(steps):
        u_pred = u + dt * v + (0.5 - beta) * dt**2 * a
        v_pred = v + (1 - gamma) * dt * a
        F = -K @ u_pred
        a = solve(K_eff, F)
        u = u_pred + beta * dt**2 * a
        v = v_pred + gamma * dt * a
        
    print("Periodic Cubic Truss Simulation Complete.")
    print(f"Final System Energy (Approx): {0.5 * np.dot(v, M @ v) + 0.5 * np.dot(u, K @ u):.2f} J")

if __name__ == "__main__":
    run_implicit_newmark_periodic()
