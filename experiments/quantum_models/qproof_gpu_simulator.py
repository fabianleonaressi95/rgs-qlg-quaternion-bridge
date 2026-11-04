import cupy as cp
import numpy as np

class QProofGPUSimulator:
    def __init__(self, n_particles, box_size=10.0, T=0.1, gamma=1.0):
        self.n_particles = n_particles
        self.box_size = box_size
        # Move state to GPU
        self.positions = cp.random.rand(n_particles, 3, dtype=cp.float32) * box_size
        self.T = T
        self.gamma = gamma

    def compute_forces(self):
        """
        Vectorized force calculation on GPU.
        For large N, this outperforms CPU by orders of magnitude.
        """
        # Broadcasting difference: (N, 1, 3) - (1, N, 3) -> (N, N, 3)
        diff = self.positions[:, cp.newaxis, :] - self.positions[cp.newaxis, :, :]
        dist = cp.linalg.norm(diff, axis=2)
        
        # Avoid self-interaction and define cutoff
        mask = (dist < 2.5) & (dist > 0.1)
        dist = cp.where(mask, dist, cp.inf)
        
        # Lennard-Jones simplified
        forces = cp.sum(diff * (1.0 / dist**3)[:, :, cp.newaxis], axis=1)
        return forces

    def step(self, dt=0.01):
        # Everything happens on GPU
        forces = self.compute_forces()
        
        # Stress-based noise (Z-Score calculation on GPU)
        # Simplified: using neighbor count as stress
        dist_sq = cp.sum((self.positions[:, cp.newaxis, :] - self.positions[cp.newaxis, :, :])**2, axis=2)
        neighbor_count = cp.sum(dist_sq < 2.5**2, axis=1)
        z_scores = (neighbor_count - cp.mean(neighbor_count)) / (cp.std(neighbor_count) + 1e-12)
        
        # Adaptive noise
        noise = cp.sqrt(self.T * (1.0 + self.gamma * cp.maximum(z_scores, 0)))[:, cp.newaxis] \
                * cp.random.normal(size=(self.n_particles, 3), dtype=cp.float32) * cp.sqrt(dt)
        
        self.positions += forces * dt + noise
        self.positions = cp.mod(self.positions, self.box_size)
        
        return self.positions

# Example Usage
if __name__ == "__main__":
    # Ensure cupy is available in the environment
    try:
        sim = QProofGPUSimulator(n_particles=1000)
        pos = sim.step()
        print(f"GPU accelerated simulation: {pos.shape[0]} particles processed on {cp.cuda.Device().name}")
    except Exception as e:
        print(f"CuPy initialization failed: {e}")
