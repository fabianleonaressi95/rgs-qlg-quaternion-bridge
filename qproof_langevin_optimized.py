import numpy as np

class QProofLangevinSimulator:
    def __init__(self, n_particles, box_size=10.0, T=0.1, gamma=1.0, cutoff=2.5):
        self.n_particles = n_particles
        self.box_size = box_size
        self.positions = np.random.rand(n_particles, 3) * box_size
        self.T = T
        self.gamma = gamma
        self.cutoff = cutoff
        # Grid setup for Cell Lists
        self.n_cells_per_dim = int(box_size / cutoff)
        if self.n_cells_per_dim < 3: self.n_cells_per_dim = 3

    def get_cell_lists(self, positions):
        """Build cell lists O(N)."""
        cells = np.floor(positions / self.cutoff).astype(int)
        # Flatten 3D cell coords to 1D
        cell_ids = cells[:, 0] + cells[:, 1] * self.n_cells_per_dim + cells[:, 2] * self.n_cells_per_dim**2
        
        # Sort particles by cell_id
        order = np.argsort(cell_ids)
        sorted_pos = positions[order]
        unique_cells, cell_starts = np.unique(cell_ids[order], return_index=True)
        cell_starts = np.append(cell_starts, self.n_particles)
        
        return sorted_pos, unique_cells, cell_starts, order

    def compute_forces_on(self, i, pos, sorted_pos, unique_cells, cell_starts):
        """Compute force for one particle using neighbors in nearby cells."""
        # This is a conceptual snippet; O(N) optimized implementation 
        # usually involves vectorized neighbor searching across grid cells.
        return np.zeros(3) # Placeholder for actual force logic

    def step(self, dt=0.01):
        sorted_pos, unique_cells, cell_starts, order = self.get_cell_lists(self.positions)
        
        # In a real O(N) implementation, we iterate over particles 
        # and only check neighbors in current and adjacent cells.
        forces = np.zeros_like(self.positions)
        # [Force calc logic using Cell Lists here]
        
        # ... (rest of the Langevin step logic as before)
        return self.positions

if __name__ == "__main__":
    sim = QProofLangevinSimulator(n_particles=1000)
    print("Cell list optimized simulator ready.")
