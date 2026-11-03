# -*- coding: utf-8 -*-
"""
Naressi Law Spin Chain Simulation Framework
Author: Fabian Leo Naressi
Description: Comprehensive simulation suite for analyzing scaling properties and 
             Von Neumann Entropy in 1D spin chains with Golden Ratio geometric coupling.
             Includes Dense QuTiP, Sparse SciPy, TeNPy MPS implementations, and 
             formal mathematical verification suites.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
# Forza matplotlib a non usare una GUI (essenziale per Docker/Code Ocean)
import matplotlib
matplotlib.use('Agg')
# Try importing QuTiP and TeNPy for completeness, but handle exceptions if not installed
try:
    import qutip as qt
except ImportError:
    qt = None

try:
    from tenpy.networks.mps import MPS
    from tenpy.models.spins import SpinChain
    from tenpy.algorithms import dmrg
except ImportError:
    MPS, SpinChain, dmrg = None, None, None

def get_naressi_coupling(N, alpha, phi=(1 + np.sqrt(5)) / 2):
    """
    Computes the geometric coupling constants according to the Naressi Law:
    J_n = alpha * phi^-n * exp(-n / phi^2)
    """
    n = np.arange(N - 1)
    return alpha * (phi**(-n)) * np.exp(-n / (phi**2))

# =============================================================================
# 1. DENSE QUTIP IMPLEMENTATION (Optimal for N <= 16)
# =============================================================================
def build_hamiltonian_qutip(N, alpha, g, phi=(1 + np.sqrt(5)) / 2):
    if qt is None:
        raise ImportError("QuTiP is required for this function.")
    J = get_naressi_coupling(N, alpha, phi)
    si = qt.qeye(2)
    sx = qt.sigmax()
    sz = qt.sigmaz()
    
    H = 0
    # Ising interaction term
    for n in range(N - 1):
        op_list = [si] * N
        op_list[n] = sx
        op_list[n+1] = sx
        H -= J[n] * qt.tensor(op_list)
        
    # Transverse field term
    for n in range(N):
        op_list = [si] * N
        op_list[n] = sz
        H -= g * qt.tensor(op_list)
        
    return H

def calculate_vn_entropy_qutip(N, alpha, g, phi=(1 + np.sqrt(5)) / 2):
    if qt is None:
        raise ImportError("QuTiP is required for this function.")
    
    # 1. Calcoliamo l'autovettore fondamentale usando SciPy Sparse (veloce e robusto)
    H_sparse = build_hamiltonian_sparse(N, alpha, g, phi)
    _, evec_sparse = eigsh(H_sparse, k=1, which='SA')
    
    # 2. Convertiamo l'autovettore in un Qobj di QuTiP (Ket dello stato fondamentale)
    dims = [[2]*N, [1]*N]
    gs = qt.Qobj(evec_sparse[:, 0], dims=dims)
    
    # 3. Calcolo dell'entropia di Von Neumann a metà catena
    rho_red = gs.ptrace(range(N // 2))
    return qt.entropy_vn(rho_red)
# =============================================================================
# 2. SPARSE SCIPY IMPLEMENTATION (Optimal for 16 < N <= 24)
# =============================================================================
def build_hamiltonian_sparse(N, alpha, g, phi=(1 + np.sqrt(5)) / 2):
    """
    Efficiently constructs the Transverse-Field Ising model Hamiltonian 
    using Compressed Sparse Row (CSR) matrices.
    """
    J = get_naressi_coupling(N, alpha, phi)
    I = sp.eye(2, format='csr')
    X = sp.csr_matrix([[0, 1], [1, 0]])
    Z = sp.csr_matrix([[1, 0], [0, -1]])
    
    H = sp.csr_matrix((2**N, 2**N))
    
    # Ising Interaction: Sum J_n * X_n * X_{n+1}
    for n in range(N - 1):
        op = 1
        for i in range(N):
            if i == n or i == n + 1:
                op = sp.kron(op, X, format='csr')
            else:
                op = sp.kron(op, I, format='csr')
        H -= J[n] * op
        
    # Transverse Field: Sum g * Z_n
    for n in range(N):
        op = 1
        for i in range(N):
            if i == n:
                op = sp.kron(op, Z, format='csr')
            else:
                op = sp.kron(op, I, format='csr')
        H -= g * op
        
    return H

def run_sparse_eigenanalysis(N, alpha, g, phi=(1 + np.sqrt(5)) / 2):
    H = build_hamiltonian_sparse(N, alpha, g, phi)
    evals, evecs = eigsh(H, k=1, which='SA') 
    return evals[0], evecs[:, 0]

# =============================================================================
# 3. MATRIX PRODUCT STATES / DMRG IMPLEMENTATION (Optimal for N >= 50)
# =============================================================================
def run_tenpy_dmrg(N=100, alpha=1.618, g=1.0, chi_max=64, phi=(1 + np.sqrt(5)) / 2):
    if dmrg is None:
        print("TeNPy is not installed. Returning conceptual placeholder.")
        return None
    
    J_list = list(get_naressi_coupling(N, alpha, phi))
    
    model_params = {
        'L': N,
        'S': 0.5,
        'J': J_list, 
        'g': g, 
        'bc_MPS': 'finite'
    }
    
    model = SpinChain(model_params)
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*N, bc=model.lat.bc_MPS)
    
    dmrg_params = {
        'mixer': True,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1.e-10}
    }
    
    info = dmrg.run(model, psi, dmrg_params)
    entropy = psi.entanglement_entropy()
    
    return info['E'], entropy

# =============================================================================
# 4. SCIENTIFIC VERIFICATION SUITE (Scale-Invariance & Falsifiability)
# =============================================================================
def verify_scale_invariance(N_sizes=[4, 8, 12, 16], alpha=1.61803398875, g=1.0):
    """
    Executes a thermodynamic finite-size analysis to map the scaling behavior
    of the ground-state energies across progressive architectural layers.
    """
    print("\n" + "="*70)
    print(" VERIFICATION SUITE I: FINITE-SIZE SCALING & RATIO ANALYSIS")
    print("="*70)
    energies = {}
    
    for N in N_sizes:
        E0, _ = run_sparse_eigenanalysis(N, alpha, g)
        energies[N] = E0
        print(f"Lattice Size N={N:2d} | Ground State Energy E0: {E0:.6f}")
        
    print("\n--- Structural Scaling Ratios ---")
    sizes = list(energies.keys())
    for idx in range(1, len(sizes)):
        n_curr = sizes[idx]
        n_prev = sizes[idx-1]
        ratio = energies[n_curr] / energies[n_prev]
        print(f"Energy Ratio E0(N={n_curr}) / E0(N={n_prev}): {ratio:.6f}")

def verify_falsifiability_perturbation(N=12, alpha=1.61803398875, g=1.0, epsilon=0.05):
    """
    Tests system resilience against geometric configuration breakdown.
    Compares ideal golden coupling vs a perturbed geometry (phi + epsilon).
    """
    print("\n" + "="*70)
    print(" VERIFICATION SUITE II: GEOMETRIC BREAKDOWN & FALSIFIABILITY")
    print("="*70)
    if qt is None:
        print("QuTiP is required to evaluate reduced density matrix entropy. Skipping.")
        return
        
    phi_ideal = (1 + np.sqrt(5)) / 2
    phi_perturbed = phi_ideal + epsilon
    
    S_ideal = calculate_vn_entropy_qutip(N, alpha, g, phi=phi_ideal)
    S_perturbed = calculate_vn_entropy_qutip(N, alpha, g, phi=phi_perturbed)
    
    print(f"Lattice Setup (N={N}, alpha={alpha}, g={g})")
    print(f"-> Ideal System Von Neumann Entropy (phi = {phi_ideal:.5f}): {S_ideal:.6f}")
    print(f"-> Perturbed System Von Neumann Entropy (phi + epsilon = {phi_perturbed:.5f}): {S_perturbed:.6f}")
    print(f"-> Net Structural Entropy Shift (Delta S): {S_perturbed - S_ideal:.6f}")

# =============================================================================
# AUTOMATION AND VISUALIZATION LOOP (QuTiP Example)
# =============================================================================
def plot_scaling_analysis():
    if qt is None:
        print("QuTiP not found. Cannot run standard automated plot.")
        return
        
    N_values = [4, 8, 12, 16]
    g = 1.0
    alpha_range = np.linspace(0.5, 2.5, 20)
    results = {n: [] for n in N_values}

    print("\nStarting automated scaling plot visualization...")
    for n in N_values:
        for a in alpha_range:
            entropy = calculate_vn_entropy_qutip(n, a, g)
            results[n].append(entropy)

    plt.figure(figsize=(8, 5))
    for n in N_values:
        plt.plot(alpha_range, results[n], label=f'N={n}')
    plt.xlabel(r'$\alpha$ (Scale Parameter)')
    plt.ylabel('Von Neumann Entropy')
    plt.title('Naressi Law: Scaling Analysis')
    plt.legend()
    plt.grid(True)
    plt.savefig('naressi_scaling_plot.png', dpi=300)
    print("Plot saved as 'naressi_scaling_plot.png'.")

if __name__ == '__main__':
    print("Naressi Law Framework Initialized successfully.")
    
    # Run the added predictive and formal check routines
    verify_scale_invariance(N_sizes=[4, 8, 12, 16])
    verify_falsifiability_perturbation(N=12, epsilon=0.05)
    
    # Optional plot execution
    plot_scaling_analysis()
    import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt

def analyze_correlations_and_potential(N, alpha, g, phi=(1 + np.sqrt(5)) / 2):
    """
    Calcola le correlazioni spaziali a lungo raggio <Z_i Z_{i+r}> 
    e stima il profilo del potenziale chimico locale sul Ground State.
    """
    # 1. Estrazione del Ground State (sfruttando il solutore sparse stabile)
    H_sparse = build_hamiltonian_sparse(N, alpha, g, phi)
    E0, evec_sparse = eigsh(H_sparse, k=2, which='SA') # Estraiamo i primi due livelli per il gap
    gs = evec_sparse[:, 0]
    E1 = E0[1]
    E0_gs = E0[0]
    
    # Il gap di eccitazione locale come proxy del potenziale chimico strutturale
    chemical_potential_proxy = E1 - E0_gs
    
    # 2. Calcolo delle correlazioni spaziali <Z_i Z_{i+r}> partendo dal centro della catena
    mid_site = N // 2
    correlations = []
    distances = list(range(1, N - mid_site))
    
    # Costruiamo l'operatore Z per il sito centrale
    # (Adatta questa logica a seconda di come generi i tuoi operatori di spin nel codice)
    # Esempio concettuale tramite prodotti tensoriali:
    def get_Z_op(site):
        # Genera l'operatore sigma_z sul sito specificato
        ops = [sp.eye(2) for _ in range(N)]
        ops[site] = sp.csr_matrix([[1, 0], [0, -1]])
        Z = ops[0]
        for op in ops[1:]:
            Z = sp.kron(Z, op, format='csr')
        return Z

    Z_mid = get_Z_op(mid_site)
    
    for r in distances:
        Z_dist = get_Z_op(mid_site + r)
        # Operatore di correlazione <Z_i Z_{i+r}>
        Corr_op = Z_mid.dot(Z_dist)
        val_corr = np.dot(gs.conj(), Corr_op.dot(gs))
        correlations.append(np.real(val_corr))
        
    return distances, correlations, chemical_potential_proxy

def plot_advanced_metrics(distances, correlations, chem_pot, N):
    """Genera e salva il grafico a due pannelli per le metriche avanzate."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Pannello 1: Correlazioni spaziali a lungo raggio
    ax1.plot(distances, correlations, 'o-', color='#1f77b4', linewidth=2, label=r'$\langle \sigma_c^z \sigma_{c+r}^z \rangle$')
    ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_title(f"Funzione di Correlazione a Lungo Raggio (N={N})")
    ax1.set_xlabel("Distanza spaziale (r)")
    ax1.set_ylabel("Ampiezza di Correlazione")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Pannello 2: Spettro/Potenziale chimico efficace
    ax2.bar(["E0 (Ground State)", "E1 (1st Excited)", "Effective Gap ($\mu$)"], 
            [E0_gs, E1, chem_pot], color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8)
    ax2.set_title("Profilo Energetico & Potenziale Chimico Efficace")
    ax2.set_ylabel("Energia / Gap")
    ax2.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plot_path = "naressi_advanced_metrics.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Advanced metrics plot saved successfully as '{plot_path}'.")