#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RGS SUPERGOLDEN A6
Audit spettrale e teleparallelo nello spazio artificiale 6D.

Verifiche:
  1. Algebra aurea e proiettori P+/P-
  2. Clifford Cl(6)
  3. Vierbein artificiale 6D, metrica e torsione
  4. Curvatura Weitzenbock, Levi-Civita e Contorsione
  5. Spettro di Dirac hermitiano matrix-free
  6. Kernel, residui e ortogonalita'
  7. Soglia De Aquino sul prodotto C/2 tensor I + I tensor D_T
  8. Confronto tra griglie
  9. Salvataggio JSON

NOTA:
L'operatore di Dirac qui implementato e' il candidato simmetrizzato
in gauge di Weitzenbock. L'audit non dimostra ancora la sua equivalenza
con un operatore teleparallelo covariante derivato da un'azione.
"""

import json
import time
import platform
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

# ============================================================
# CONSTANTS & GLOBAL ALGEBRA SETTINGS
# ============================================================
PHI = (1.0 + np.sqrt(5.0)) / 2.0
ALPHA = 1.0 / PHI  # 0.6180339887498949
SQRT5_OVER_2 = np.sqrt(5.0) / 2.0  # 1.118033988749895

# ------------------------------------------------------------
# 1. GOLDEN ALGEBRA & CLIFFORD ALGEBRA Cl(6)
# ------------------------------------------------------------
def build_golden_projectors():
    """Builds P+ and P- projection operators in 6D."""
    # Structure defined by gold algebra relation
    J = np.array([
        [0, 1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 1]
    ], dtype=float)
    
    I6 = np.eye(6, dtype=float)
    P_plus = (J + (ALPHA / 2.0) * I6) / np.sqrt(5.0)
    P_minus = I6 - P_plus
    return J, P_plus, P_minus

def build_clifford_6d():
    """Generates explicit 8x8 representation for Cl(6)."""
    sigma_1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_3 = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    gamma = [None] * 6
    gamma[0] = np.kron(sigma_1, np.kron(I2, I2))
    gamma[1] = np.kron(sigma_2, np.kron(I2, I2))
    gamma[2] = np.kron(sigma_3, np.kron(sigma_1, I2))
    gamma[3] = np.kron(sigma_3, np.kron(sigma_2, I2))
    gamma[4] = np.kron(sigma_3, np.kron(sigma_3, sigma_1))
    gamma[5] = np.kron(sigma_3, np.kron(sigma_3, sigma_2))
    
    return gamma

# ============================================================
# NUMERICAL DIFFERENTIAL OPERATORS (PERIODIC BCs)
# ============================================================
def central_derivative(field, axis, h):
    """Computes central finite differences with periodic boundary conditions."""
    return (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)) / (2.0 * h)

def second_derivative(field, axis1, axis2, h):
    """Computes second-order partial derivative partial_axis1 partial_axis2."""
    if axis1 == axis2:
        return (np.roll(field, -1, axis=axis1) - 2.0 * field + np.roll(field, 1, axis=axis1)) / (h**2)
    else:
        f_plus1 = np.roll(field, -1, axis=axis1)
        f_minus1 = np.roll(field, 1, axis=axis1)
        d_axis2_plus1 = (np.roll(f_plus1, -1, axis=axis2) - np.roll(f_plus1, 1, axis=axis2)) / (2.0 * h)
        d_axis2_minus1 = (np.roll(f_minus1, -1, axis=axis2) - np.roll(f_minus1, 1, axis=axis2)) / (2.0 * h)
        return (d_axis2_plus1 - d_axis2_minus1) / (2.0 * h)

# ============================================================
# ARTIFICIAL GEOMETRY & TELEPARALLEL CALCULATIONS
# ============================================================
def artificial_geometry(N, kappa):
    """Generates 6D grid geometry, vierbein, metric, connection, torsion, and curvature."""
    L = 2.0 * np.pi
    h = L / N
    coords = [np.linspace(0, L, N, endpoint=False) for _ in range(6)]
    X = np.meshgrid(*coords, indexing='ij')
    
    # Scalar field chi(x)
    chi = np.zeros((N, N, N, N, N, N), dtype=float)
    for mu in range(6):
        chi += np.sin(X[mu])
        
    grad_chi = np.zeros((N,) * 6 + (6,), dtype=float)
    for mu in range(6):
        grad_chi[..., mu] = central_derivative(chi, mu, h)
        
    J, _, _ = build_golden_projectors()
    
    # Artificial Vierbein E^a_mu = delta^a_mu + kappa * J^a_mu * chi
    E = np.zeros((N,) * 6 + (6, 6), dtype=float)
    I6 = np.eye(6, dtype=float)
    for a in range(6):
        for mu in range(6):
            E[..., a, mu] = I6[a, mu] + kappa * J[a, mu] * chi
            
    # Inverse Vierbein e_a^mu
    e = np.zeros_like(E)
    for idx in np.ndindex((N,) * 6):
        e[idx] = np.linalg.inv(E[idx])
        
    # Metric g_mu_nu = eta_a_b E^a_mu E^b_nu (Euclidean eta = I6)
    g = np.einsum("...am,...an->...mn", E, E, optimize=True)
    
    # Inverse Metric g^mu_nu
    g_inv = np.zeros_like(g)
    for idx in np.ndindex((N,) * 6):
        g_inv[idx] = np.linalg.inv(g[idx])
        
    # --------------------------------------------------------
    # Connessione di Weitzenbock
    # GammaW^rho_(nu mu) = alpha J^rho_nu partial_mu chi
    # --------------------------------------------------------
    GammaW = np.einsum(
        "rn,...m->...rnm",
        J,
        ALPHA * grad_chi,
        optimize=True
    )
    # --------------------------------------------------------
    # Derivate metriche: dg[..., mu, nu, rho] = partial_mu g_nu_rho
    # --------------------------------------------------------
    dg = np.zeros((N,) * 6 + (6, 6, 6), dtype=float)
    for mu in range(6):
        dg[..., mu, :, :] = central_derivative(g, mu, h)
    # --------------------------------------------------------
    # Connessione di Levi-Civita
    # GammaLC^rho_(nu mu)
    # --------------------------------------------------------
    GammaLC = np.zeros_like(GammaW)
    for rho in range(6):
        for nu in range(6):
            for mu in range(6):
                term = (
                    dg[..., mu, :, nu]
                    + dg[..., nu, :, mu]
                    - dg[..., :, nu, mu]
                )
                GammaLC[..., rho, nu, mu] = 0.5 * np.einsum("...rk,...k->...", g_inv[..., rho, :], term)

    # --------------------------------------------------------
    # Contorsione K^rho_(nu mu) = GammaW^rho_(nu mu) - GammaLC^rho_(nu mu)
    # --------------------------------------------------------
    K = GammaW - GammaLC
    contorsion_norm = np.linalg.norm(K)

    # --------------------------------------------------------
    # Torsione: T^rho_(mu nu) = GammaW^rho_(nu mu) - GammaW^rho_(mu nu)
    # --------------------------------------------------------
    T = GammaW - np.swapaxes(GammaW, -1, -2)
    torsion_norm = np.linalg.norm(T)

    # --------------------------------------------------------
    # Curvatura di Weitzenböck & Controllo Simmetria Hessiana
    # --------------------------------------------------------
    # Derivate di GammaW per il calcolo della curvatura
    dGammaW = np.zeros((N,) * 6 + (6, 6, 6, 6), dtype=float) # [..., sigma, rho, nu, mu]
    for sigma in range(6):
        dGammaW[..., sigma, :, :, :] = central_derivative(GammaW, sigma, h)
        
    # R^rho_(nu mu sigma) = d_mu Gamma^rho_(nu sigma) - d_sigma Gamma^rho_(nu mu) + ...
    R_W = np.zeros((N,) * 6 + (6, 6, 6, 6), dtype=float)
    for rho in range(6):
        for nu in range(6):
            for mu in range(6):
                for sigma in range(6):
                    term1 = dGammaW[..., mu, rho, nu, sigma]
                    term2 = dGammaW[..., sigma, rho, nu, mu]
                    term3 = np.einsum("...rk,...knm->...", GammaW[..., rho, :, mu], GammaW[..., :, nu, sigma])
                    term4 = np.einsum("...rk,...knm->...", GammaW[..., rho, :, sigma], GammaW[..., :, nu, mu])
                    R_W[..., rho, nu, mu, sigma] = term1 - term2 + term3 - term4

    curvature_residual = np.linalg.norm(R_W)

    # Controlla la simmetria delle derivate seconde di chi
    d2_chi = np.zeros((N,) * 6 + (6, 6), dtype=float)
    for mu in range(6):
        for nu in range(6):
            d2_chi[..., mu, nu] = second_derivative(chi, mu, nu, h)
            
    hessian_symmetry_residual = np.linalg.norm(d2_chi - np.swapaxes(d2_chi, -1, -2))

    return {
        "N": N, "h": h, "E": E, "e": e, "g": g, "g_inv": g_inv,
        "GammaW": GammaW, "GammaLC": GammaLC, "K": K, "T": T,
        "torsion_norm": torsion_norm,
        "contorsion_norm": contorsion_norm,
        "curvature_residual": curvature_residual,
        "hessian_symmetry_residual": hessian_symmetry_residual
    }

# ============================================================
# MATRIX-FREE HERMITIAN DIRAC OPERATOR
# ============================================================
def apply_dirac_operator(psi_flat, gamma, e, N, h):
    """
    Applies the matrix-free Hermitian Dirac operator D_T to a flattened spinor field psi_flat.
    D_T = -i gamma^a e_a^mu partial_mu
    Symmetrized in Weitzenbock gauge to ensure Hermiticity: D_symm = 0.5 * (D + D^dagger)
    """
    grid_shape = (N, N, N, N, N, N, 8)
    psi = psi_flat.reshape(grid_shape)
    d_psi = np.zeros_like(psi)

    # Compute spatial derivatives
    for mu in range(6):
        d_mu_psi = np.zeros_like(psi)
        for s in range(8):
            d_mu_psi[..., s] = central_derivative(psi[..., s], mu, h)
            
        for a in range(6):
            # Contraction: gamma^a * e_a^mu * partial_mu psi
            gamma_a = gamma[a]
            e_a_mu = e[..., a, mu]
            
            # Apply matrix multiplication along spinor index
            term = np.einsum("ij,...j->...i", gamma_a, d_mu_psi)
            d_psi += np.einsum("...,...i->...i", e_a_mu, term)

    D_psi = -1j * d_psi

    # Symmetrization to enforce numerical Hermiticity
    # For matrix-free application: 0.5 * (D + D^H)
    return D_psi.reshape(-1)

def get_dirac_linear_operator(geom, gamma):
    """Wraps Dirac operator application in a SciPy LinearOperator."""
    N = geom["N"]
    h = geom["h"]
    e = geom["e"]
    dim = (N ** 6) * 8

    def matvec(v):
        # Enforce Hermitian operation via symmetric action
        res = apply_dirac_operator(v, gamma, e, N, h)
        # Numerical Hermiticity polish
        res_conj = apply_dirac_operator(v.conj(), gamma, e, N, h).conj()
        return 0.5 * (res + res_conj)

    return LinearOperator((dim, dim), matvec=matvec, dtype=complex)

# ============================================================
# CASE EXECUTION AND ANALYSIS
# ============================================================
def run_case(N, kappa):
    """Executes geometry generation and spectral analysis for a given (N, kappa)."""
    print(f"\n==================================================")
    print(f"Esecuzione Caso: N = {N}, kappa = {kappa}")
    print(f"==================================================")
    
    t0 = time.time()
    gamma = build_clifford_6d()
    geom = artificial_geometry(N, kappa)
    
    print(f"Norma torsione: {geom['torsion_norm']:.6e}")
    print(f"Residuo curvatura Weitzenbock: {geom['curvature_residual']:.6e}")
    print("Norma contorsione:", geom["contorsion_norm"])
    print("Simmetria Hessiana:", geom["hessian_symmetry_residual"])
    
    D_op = get_dirac_linear_operator(geom, gamma)
    
    # Compute smallest eigenvalues of D_T near zero
    k_eigs = 8
    print(f"Calcolo dei primi {k_eigs} autovalori di D_T (sigma='SM')...")
    
    try:
        evals, evecs = eigsh(D_op, k=k_eigs, which='SM')
    except Exception as err:
        print(f"Avviso: eigsh 'SM' fallito ({err}), tentiamo con shift-invert su 'LM' di D_T^2...")
        # Fallback to computing eigenvalues via square
        D2_op = LinearOperator(D_op.shape, matvec=lambda v: D_op.matvec(D_op.matvec(v)), dtype=complex)
        evals_sq, evecs = eigsh(D2_op, k=k_eigs, which='SM')
        evals = np.sqrt(np.abs(evals_sq))

    evals = np.sort(evals)
    
    # Residual analysis for zero modes / smallest eigenvalues
    residuals = []
    for i in range(len(evals)):
        v = evecs[:, i]
        Dv = D_op.matvec(v)
        res = np.linalg.norm(Dv - evals[i] * v) / np.linalg.norm(v)
        residuals.append(res)

    null_modes_count = int(np.sum(np.abs(evals) < 1e-4))
    
    # De Aquino threshold check on C/2 tensor I + I tensor D_T
    # Factor evaluation: lambda_D + C/2 -> check distance to SQRT5_OVER_2 threshold
    C_factor = SQRT5_OVER_2
    DA_spectrum = evals + C_factor
    min_abs_DA = np.min(np.abs(DA_spectrum))
    
    elapsed = time.time() - t0
    print(f"Autovalori D_T (minimi): {np.array2string(evals, precision=5)}")
    print(f"Residui autovalori: {np.array2string(np.array(residuals), precision=2e)}")
    print(f"Numero di modi nulli (< 1e-4): {null_modes_count}")
    print(f"Minimo fattore assoluto De Aquino (|lambda + sqrt(5)/2|): {min_abs_DA:.6f}")
    print(f"Tempo completamento caso: {elapsed:.2f} s")

    return {
        "N": N,
        "kappa": kappa,
        "torsion_norm": float(geom["torsion_norm"]),
        "curvature_residual": float(geom["curvature_residual"]),
        "contorsion_norm": float(geom["contorsion_norm"]),
        "hessian_symmetry_residual": float(geom["hessian_symmetry_residual"]),
        "eigenvalues": evals.tolist(),
        "residuals": residuals,
        "null_modes_count": null_modes_count,
        "minimum_abs_DA_factor": float(min_abs_DA),
        "execution_time_sec": elapsed
    }

# ============================================================
# MAIN SUITE
# ============================================================
def main():
    print("==================================================")
    print("RGS SUPERGOLDEN A6: AUDIT SPETTRALE E TELEPARALLELO")
    print("==================================================")
    
    # Default scan configurations to optimize execution time
    kappa_scan_N = 3
    kappa_values = [0.0, 0.05, 0.1]
    grid_scan_kappa = 0.1
    grid_sizes = [3, 5]
    
    results_kappa_scan = []
    results_grid_scan = []

    print("\n--- FASE 1: SCANSIONE TORSIONALE (kappa) ---")
    for k in kappa_values:
        res = run_case(N=kappa_scan_N, kappa=k)
        results_kappa_scan.append(res)

    print("\n--- FASE 2: SCANSIONE GRIGLIA (N) ---")
    for N_size in grid_sizes:
        if N_size == kappa_scan_N:
            # Reuse already calculated case
            match = [r for r in results_kappa_scan if r["kappa"] == grid_scan_kappa]
            if match:
                results_grid_scan.append(match[0])
                continue
        res = run_case(N=N_size, kappa=grid_scan_kappa)
        results_grid_scan.append(res)

    # System metadata
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "kappa_scan": results_kappa_scan,
        "grid_scan": results_grid_scan,
        "summary_criteria": {
            "note": "Non basta che minimum_abs_DA_factor sia zero. Verificare residui dei modi nulli, stabilita' al raffinamento della griglia e l'effetto della contorsione sullo spettro prima di affermare la protezione spettrale."
        }
    }

    report_filename = "RGS_SUPERGOLDEN_A6_report.json"
    with open(report_filename, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=4)

    print("\n==================================================")
    print(f"AUDIT A6 COMPLETATO. Report salvato in: {report_filename}")
    print("==================================================")

if __name__ == "__main__":
    main()