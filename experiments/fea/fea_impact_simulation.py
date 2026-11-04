import numpy as np
import pandas as pd

def run_impact_simulation():
    # Load geometry
    nodes = pd.read_csv("/workspace/qc_nodes.csv", header=None).values
    edges = pd.read_csv("/workspace/qc_edges.csv", header=None).values
    
    # Simulation Parameters
    num_nodes = len(nodes)
    mass = 0.01  # kg per node
    k = 1e1      # Drastically reduced stiffness to stabilize
    damping_alpha = 5.0 # Increased damping
    dt = 1e-6    # Even smaller timestep
    steps = 500
    
    # Initialize velocities/positions
    vel = np.zeros((num_nodes, 3))
    pos = nodes.copy()
    
    # Apply impact at top nodes (z-max)
    z_max = np.max(pos[:, 2])
    impact_nodes = np.where(pos[:, 2] > (z_max - 0.2))[0]
    vel[impact_nodes, 2] = -10.0  # Impact velocity
    
    # Damping parameters
    damping_alpha = 0.5 # Velocity-proportional damping
    
    # Simulation Loop
    energy_history = []
    
    for _ in range(steps):
        accel = np.zeros((num_nodes, 3))
        for edge in edges:
            i, j = edge
            diff = pos[i] - pos[j]
            dist = np.linalg.norm(diff)
            if dist < 1e-6: continue
            
            # Simple linear spring
            force = -k * (dist - 0.2) * (diff / dist)
            accel[i] += force / mass
            accel[j] -= force / mass
            
        # Add damping force
        accel -= damping_alpha * vel
            
        vel += accel * dt
        pos += vel * dt
        
        ke = 0.5 * mass * np.sum(vel**2)
        energy_history.append(ke)
        
    print("Simulation Complete.")
    print(f"Initial Energy: {energy_history[0]:.2f} J")
    print(f"Final Energy: {energy_history[-1]:.2f} J")
    print(f"Energy Dissipated: {energy_history[0] - energy_history[-1]:.2f} J")
    
    # Save deformation
    np.savetxt("/workspace/qc_deformed_nodes.csv", pos, delimiter=",")

if __name__ == "__main__":
    run_impact_simulation()
