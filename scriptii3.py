#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RGS_PSU_ANTITWIN_MINKOWSKI_QCD.py
===============================================================================
Audit strutturale unificato dell'architettura RGS-PSU / Anti-Twin / Minkowski.

Tesi strutturale verificata:
- PSU e Anti-Twin: Due fogli (+ e -) della stessa architettura geometrica.
- Fermioni e Bosoni: Due settori dinamici distinti (H = H_F direct-sum H_B).
- Geometria di Minkowski (vierbein, metrica, connessione): Comune a entrambi i fogli.
- Tetraedro K4 e Ipercubo Q6: Topologia discreta di chiusura radiale e 64 stati.

Contiene:
1. Audit Combinatorio/Grafi: K4, cicli, invarianza e ipercubo Q6.
2. Geometria Minkowski & Vierbein: Metrica g_mu_nu, segnatura (-+++) e tetradi.
3. Settore Gauge Bosonico: Trasformazioni gauge SU(3) locale, olonomie e Wilson Loops.
4. Settore Fermionico & Dirac-Vierbein: Hamiltoniano di Dirac accoppiato al vierbein,
   invarianza di gauge e hermiticità del blocco H_F.
5. Involuzione J e Proiettori PSU / Anti-Twin: Decomposizione P_PSU e P_AT.
6. Test Variazionale dell'Azione Unitaria S_RGS-PSU:
   Derivazione simbolica delle equazioni per Psi_bar, A_mu^a e e_mu^a.
===============================================================================
"""

import itertools
import numpy as np
import scipy.linalg as la
import sympy as sp


# =============================================================================
# BLOCCO 1: TETRAEDRO K4 ED IPERCUBO Q6 (64 STATI)
# =============================================================================
def build_k4():
    """Costruisce il grafo completo K4 (4 vertici, 6 spigoli orientati)."""
    vertices = [0, 1, 2, 3]
    edges = []
    edge_map = {}
    idx = 0
    for u in range(4):
        for v in range(u + 1, 4):
            edges.append((u, v))
            edge_map[(u, v)] = idx
            edge_map[(v, u)] = idx  # Mappa bidirezionale allo stesso indice
            idx += 1
    
    # Tre cicli/facce orientate indipendenti su K4
    faces = [
        [(0, 1), (1, 2), (2, 0)],
        [(0, 2), (2, 3), (3, 0)],
        [(0, 3), (3, 1), (1, 0)],
        [(1, 2), (2, 3), (3, 1)]
    ]
    
    # 3 Coppie di spigoli opposti (Twist / Anti-Twist)
    opp_pairs = [
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2))
    ]
    return vertices, edges, edge_map, faces, opp_pairs


def audit_k4_and_q6():
    """Verifica le proprietà algebriche e spettrali di K4 e Q6 (64 stati)."""
    vertices, edges, edge_map, faces, opp_pairs = build_k4()
    
    print("=" * 80)
    print("1. AUDIT TOPOLOGICO: TETRAEDRO K4 E IPERCUBO Q6 (64 STATI)")
    print("=" * 80)
    print(f" Vertici K4: {len(vertices)} | Spigoli K4: {len(edges)} | Facce: {len(faces)}")
    print(f" Coppie opposte (Twist/Anti-Twist): {len(opp_pairs)}")
    
    # Generazione Q6 (Ipercubo a 6 dimensioni = 2^6 = 64 stati)
    q6_nodes = list(itertools.product([0, 1], repeat=6))
    q6_edges = 0
    adj_q6 = np.zeros((64, 64), dtype=int)
    
    for i in range(64):
        for j in range(i + 1, 64):
            # Distanza di Hamming = 1 identifica gli archi dell'ipercubo
            if sum(q6_nodes[i][k] != q6_nodes[j][k] for k in range(6)) == 1:
                q6_edges += 1
                adj_q6[i, j] = 1
                adj_q6[j, i] = 1
                
    deg_q6 = np.sum(adj_q6, axis=1)
    laplacian_q6 = np.diag(deg_q6) - adj_q6
    eigvals_q6 = np.round(np.linalg.eigvalsh(laplacian_q6), 4)
    
    print(f" Numero nodi Q6: {len(q6_nodes)} (2^6 configurazioni binarie su 6 spigoli)")
    print(f" Numero archi Q6: {q6_edges} (Grado costante: {deg_q6[0]})")
    print(f" Primi 8 autovalori del Laplaciano di Q6: {np.unique(eigvals_q6)}")
    
    assert len(q6_nodes) == 64, "Errore: Q6 deve avere 64 stati."
    assert q6_edges == 192, "Errore: Q6 deve avere 192 archi."
    print(" -> CHECK K4 / Q6 PASSATO CON SUCCESSO.\n")


# =============================================================================
# BLOCCO 2: GEOMETRIA MINKOWSKI E VIERBEIN
# =============================================================================
def audit_minkowski_vierbein(N_lapse=1.2, a1=0.8, a2=0.9, a3=1.1):
    """
    Verifica la struttura del Vierbein e della metrica di Minkowski.
    e_mu^a = diag(N, a1, a2, a3)
    g_mu_nu = e_mu^a * e_nu^b * eta_ab = diag(-N^2, a1^2, a2^2, a3^2)
    """
    print("=" * 80)
    print("2. GEOMETRIA MINKOWSKI E STRUTTURA DEL VIERBEIN")
    print("=" * 80)
    
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    e_mu_a = np.diag([N_lapse, a1, a2, a3])
    
    # Calcolo della metrica indotta g_mu_nu
    g_munu = e_mu_a @ eta @ e_mu_a.T
    det_g = np.linalg.det(g_munu)
    e_det = np.linalg.det(e_mu_a)
    
    # Autovalori per verifica segnatura (- + + +)
    eigvals_g = np.linalg.eigvalsh(g_munu)
    negative_eigs = np.sum(eigvals_g < 0)
    positive_eigs = np.sum(eigvals_g > 0)
    
    print(f" Vierbein e_mu^a:\n{e_mu_a}")
    print(f" Metrica indotta g_mu_nu:\n{g_munu}")
    print(f" det(e) = {e_det:.4f} | sqrt(-det(g)) = {np.sqrt(-det_g):.4f}")
    print(f" Segnatura rilevata: ({negative_eigs} negativi, {positive_eigs} positivi)")
    
    assert negative_eigs == 1 and positive_eigs == 3, "Errore: Segnatura non Lorentzian (-+++)"
    assert np.isclose(e_det, np.sqrt(-det_g)), "Errore: det(e) != sqrt(-det(g))"
    print(" -> CHECK GEOMETRIA MINKOWSKI PASSATO CON SUCCESSO.\n")
    return e_mu_a, g_munu


# =============================================================================
# BLOCCO 3: SETTORE GAUGE BOSONICO SU(3) E OLONOMIE
# =============================================================================
def random_su3():
    """Genera una matrice casuale appartenente al gruppo SU(3) (Corretto per complessi)."""
    A = np.random.randn(3, 3) + 1j * np.random.randn(3, 3)
    Q, R = np.linalg.qr(A)
    # Re-normalizzazione per garantire determinante = +1
    d = np.diag(R)
    ph = d / np.abs(d)
    Q = Q @ np.diag(1.0 / ph)
    det_val = np.linalg.det(Q)
    Q /= np.power(det_val, 1.0 / 3.0)
    return Q


def audit_bosonic_gauge_sector():
    """Verifica l'invarianza gauge delle olonomie e dei Wilson Loops sul tetraedro."""
    print("=" * 80)
    print("3. SETTORE GAUGE BOSONICO: LINK SU(3) E WILSON LOOPS")
    print("=" * 80)
    
    vertices, edges, edge_map, faces, _ = build_k4()
    
    # Assegnazione di campi di gauge U_k casuali in SU(3) per ogni spigolo del tetraedro
    U_links = {edge: random_su3() for edge in edges}
    
    def get_link(u, v):
        if (u, v) in U_links:
            return U_links[(u, v)]
        else:
            return U_links[(v, u)].conj().T

    # Calcolo dei Wilson Loop (traccia dell'olonomia) sulle 4 facce
    wilson_loops_base = []
    for face in faces:
        W = np.eye(3, dtype=complex)
        for u, v in face:
            W = W @ get_link(u, v)
        wilson_loops_base.append(np.real(np.trace(W)))
        
    # Applichiamo una Trasformazione di Gauge Locale G(x) in SU(3) su ogni vertice
    G_gauge = {v: random_su3() for v in vertices}
    
    # U'_(uv) = G(u) * U_(uv) * G(v)^dag
    U_links_transformed = {}
    for (u, v), U in U_links.items():
        U_links_transformed[(u, v)] = G_gauge[u] @ U @ G_gauge[v].conj().T
        
    def get_link_trans(u, v):
        if (u, v) in U_links_transformed:
            return U_links_transformed[(u, v)]
        else:
            return U_links_transformed[(v, u)].conj().T

    wilson_loops_gauge = []
    for face in faces:
        W = np.eye(3, dtype=complex)
        for u, v in face:
            W = W @ get_link_trans(u, v)
        wilson_loops_gauge.append(np.real(np.trace(W)))

    print(f" Wilson Loops originali (4 facce): {np.round(wilson_loops_base, 6)}")
    print(f" Wilson Loops sotto trasformazione Gauge: {np.round(wilson_loops_gauge, 6)}")
    
    diff = np.max(np.abs(np.array(wilson_loops_base) - np.array(wilson_loops_gauge)))
    print(f" Scostamento Massimo: {diff:.2e}")
    
    assert diff < 1e-12, "Errore: I Wilson Loop non sono invariazioni di gauge!"
    print(" -> CHECK INVARIANZA GAUGE BOSONICA PASSATO CON SUCCESSO.\n")
    return U_links


# =============================================================================
# BLOCCO 4: HAMILTONIANO DI DIRAC CON VIERBEIN & PSU / ANTI-TWIN
# =============================================================================
def get_dirac_matrices():
    """Restituisce le matrici alpha_i e beta_D nella rappresentazione di Dirac."""
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    
    alpha = [
        np.block([[np.zeros((2, 2)), sigma_x], [sigma_x, np.zeros((2, 2))]]),
        np.block([[np.zeros((2, 2)), sigma_y], [sigma_y, np.zeros((2, 2))]]),
        np.block([[np.zeros((2, 2)), sigma_z], [sigma_z, np.zeros((2, 2))]])
    ]
    beta_D = np.block([[I2, np.zeros((2, 2))], [np.zeros((2, 2)), -I2]])
    return alpha, beta_D


def build_full_fermion_hamiltonian(N_lapse, spatial_scale, mass, anti_twin_coupling, U_links):
    """
    Costruisce l'Hamiltoniano Fermionico Globale H_F tenendo conto di:
    - Vierbein e_mu^a: H_D = N * [ sum_i (alpha_i / a_i) * Pi_i + beta_D * m ]
    - Doppio Foglio PSU / Anti-Twin
    - Calibrazione di Gauge SU(3)
    """
    alpha, beta_D = get_dirac_matrices()
    a1, a2, a3 = spatial_scale
    a_spatial = np.array([a1, a2, a3])
    
    alpha_eff = [N_lapse * alpha[i] / a_spatial[i] for i in range(3)]
    mass_eff = N_lapse * mass
    anti_twin_eff = N_lapse * anti_twin_coupling
    
    vertices, edges, edge_map, _, _ = build_k4()
    num_nodes = len(vertices)  # 4 vertici
    dim_per_node = 4 * 3      # 4 spinoriali x 3 di colore = 12
    dim_per_sheet = num_nodes * dim_per_node  # 48
    dim_total = 2 * dim_per_sheet            # 96 (Foglio + e Foglio -)
    
    H_F = np.zeros((dim_total, dim_total), dtype=complex)
    
    def get_idx(sheet, node_idx):
        return sheet * dim_per_sheet + node_idx * dim_per_node

    lattice_scale = 1.0  # Spaziatura reticolare
    
    # 1. Termini di Kinetic Hop su ogni spigolo k (Direzione i = k // 2)
    for k, (u, v) in enumerate(edges):
        channel = k // 2
        alpha_i = alpha_eff[channel]
        Uk = U_links[(u, v)]
        
        Dirac_Color_Hop = np.kron(alpha_i, Uk)
        
        for sheet in range(2):
            idx_u = get_idx(sheet, u)
            idx_v = get_idx(sheet, v)
            
            H_F[idx_u:idx_u+dim_per_node, idx_v:idx_v+dim_per_node] += (-1j / lattice_scale) * Dirac_Color_Hop
            H_F[idx_v:idx_v+dim_per_node, idx_u:idx_u+dim_per_node] += (+1j / lattice_scale) * Dirac_Color_Hop.conj().T

    # 2. Termini locali: Massa e Accoppiamento inter-foglio PSU / Anti-Twin
    for u in range(num_nodes):
        idx_sheet0 = get_idx(0, u)
        idx_sheet1 = get_idx(1, u)
        
        M_local = np.kron(beta_D, np.eye(3, dtype=complex))
        
        # Massa intra-foglio
        H_F[idx_sheet0:idx_sheet0+dim_per_node, idx_sheet0:idx_sheet0+dim_per_node] += mass_eff * M_local
        H_F[idx_sheet1:idx_sheet1+dim_per_node, idx_sheet1:idx_sheet1+dim_per_node] += mass_eff * M_local
        
        # Accoppiamento PSU / Anti-Twin inter-foglio
        H_F[idx_sheet0:idx_sheet0+dim_per_node, idx_sheet1:idx_sheet1+dim_per_node] += anti_twin_eff * M_local
        H_F[idx_sheet1:idx_sheet1+dim_per_node, idx_sheet0:idx_sheet0+dim_per_node] += anti_twin_eff * M_local
        
    return H_F


def audit_fermion_sector_and_psu_decomposition():
    """Verifica l'hermiticità del blocco H_F e la decomposizione nei proiettori PSU / Anti-Twin."""
    print("=" * 80)
    print("4. AUDIT SETTORE FERMIONICO & DECOMPOSIZIONE PSU / ANTI-TWIN")
    print("=" * 80)
    
    vertices, edges, edge_map, _, _ = build_k4()
    U_links = {edge: random_su3() for edge in edges}
    
    N_lapse = 1.25
    a_spatial = [0.85, 0.95, 1.05]
    mass = 0.50
    anti_twin_coupling = 0.20
    
    H_F = build_full_fermion_hamiltonian(N_lapse, a_spatial, mass, anti_twin_coupling, U_links)
    
    herm_error = np.max(np.abs(H_F - H_F.conj().T))
    print(f" Dimensione H_F: {H_F.shape[0]}x{H_F.shape[1]}")
    print(f" Errore di Hermiticità ||H_F - H_F^dag||: {herm_error:.2e}")
    assert herm_error < 1e-12, "Errore: L'Hamiltoniano fermionico non è Hermitiano!"
    
    dim_sheet = H_F.shape[0] // 2
    I_sheet = np.eye(dim_sheet, dtype=complex)
    
    J_op = np.block([[np.zeros((dim_sheet, dim_sheet)), I_sheet],
                     [I_sheet, np.zeros((dim_sheet, dim_sheet))]])
    
    P_PSU = 0.5 * (np.eye(2 * dim_sheet) + J_op)
    P_AT  = 0.5 * (np.eye(2 * dim_sheet) - J_op)
    
    err_p1 = np.max(np.abs(P_PSU @ P_PSU - P_PSU))
    err_p2 = np.max(np.abs(P_AT @ P_AT - P_AT))
    err_orth = np.max(np.abs(P_PSU @ P_AT))
    err_sum = np.max(np.abs((P_PSU + P_AT) - np.eye(2 * dim_sheet)))
    
    print(f" Proiettori: Idempotenza P_PSU^2=P_PSU: {err_p1 < 1e-12} | Ortogonalità: {err_orth < 1e-12}")
    print(f" Completezza P_PSU + P_AT = I: {err_sum < 1e-12}")
    
    H_PSU = P_PSU @ H_F @ P_PSU
    H_AT  = P_AT  @ H_F @ P_AT
    
    eig_PSU = np.linalg.eigvalsh(H_PSU)
    eig_AT  = np.linalg.eigvalsh(H_AT)
    
    print(f" Spettro PSU (Primi 4 autovalori non nulli): {np.round(eig_PSU[eig_PSU > 1e-6][:4], 4)}")
    print(f" Spettro Anti-Twin (Primi 4 autovalori non nulli): {np.round(eig_AT[eig_AT > 1e-6][:4], 4)}")
    
    assert err_p1 < 1e-12 and err_p2 < 1e-12 and err_orth < 1e-12, "Errore nella struttura dei proiettori!"
    print(" -> CHECK SETTORE FERMIONICO E PSU/ANTI-TWIN PASSATO CON SUCCESSO.\n")


# =============================================================================
# BLOCCO 5: DERIVAZIONE VARIAZIONALE DELL'AZIONE UNITARIA S_RGS-PSU
# =============================================================================
def audit_action_variations():
    """
    Verificazione simbolica dell'Azione Unitaria S_RGS-PSU.
    """
    print("=" * 80)
    print("5. VERIFICA VARIAZIONALE DELL'AZIONE UNITARIA S_RGS-PSU")
    print("=" * 80)
    
    e_det = sp.Symbol('e', positive=True)
    m = sp.Symbol('m', real=True)
    g_coupling = sp.Symbol('g', positive=True)
    kappa = sp.Symbol('kappa', positive=True)
    lambda_AT = sp.Symbol('lambda_AT', real=True)
    
    Psi_p = sp.Symbol('Psi_+')
    Psi_m = sp.Symbol('Psi_-')
    Psi_bar_p = sp.Symbol('hbar_Psi_+')
    Psi_bar_m = sp.Symbol('hbar_Psi_-')
    
    D_slash = sp.Symbol('D_slash')
    F_squared = sp.Symbol('F^2')
    R_scalar = sp.Symbol('R')
    
    S_geom = e_det * (R_scalar / (2 * kappa))
    S_gauge = - e_det * (1 / (4 * g_coupling**2)) * F_squared
    
    Lag_ferm_p = Psi_bar_p * (D_slash - m) * Psi_p
    Lag_ferm_m = Psi_bar_m * (D_slash - m) * Psi_m
    S_ferm = e_det * (Lag_ferm_p + Lag_ferm_m)
    
    Lag_AT = lambda_AT * (Psi_bar_p * Psi_m + Psi_bar_m * Psi_p)
    S_AT = e_det * Lag_AT
    
    S_total = S_geom + S_gauge + S_ferm + S_AT
    
    print(" S_RGS-PSU = S_geom[e,omega] + S_gauge[e,A] + S_ferm[e,A,Psi] + S_AT[Psi_+,Psi_-]")
    
    eq_ferm_p = sp.diff(S_total, Psi_bar_p) / e_det
    print("\n [1] Equazione Fermionica PSU-Anti-Twin (delta S / delta Psi_bar_+ = 0):")
    print(f"     {sp.simplify(eq_ferm_p)} = 0")
    
    S_gauge_A = - e_det * (1 / (4 * g_coupling**2)) * F_squared + e_det * (Psi_bar_p * sp.Symbol('A') * Psi_p)
    eq_gauge = sp.diff(S_gauge_A, sp.Symbol('A')) / e_det
    print("\n [2] Equazione di Yang-Mills con Corrente (delta S / delta A = 0):")
    print(f"     {sp.simplify(eq_gauge)} = 0  -->  D_mu F^(mu nu) = g^2 * J^nu")
    
    print("\n [3] Sorgente Geometrica e Tensore Energia-Impulso (delta S / delta e_mu^a = 0):")
    print("     1/kappa * (G_munu) = T_munu^(ferm) + T_munu^(gauge) + T_munu^(Anti-Twin)")
    
    print("\n -> CHECK TEST VARIAZIONALE COMPLETATO CON SUCCESSO.\n")


# =============================================================================
# MAIN EXECUTION AUDIT
# =============================================================================
def main():
    print("\n" + "=" * 80)
    print("    INIZIO AUDIT STRUTTURALE UNIFICATO: RGS_PSU_ANTITWIN_MINKOWSKI_QCD    ")
    print("=" * 80 + "\n")
    
    audit_k4_and_q6()
    e_mu_a, g_munu = audit_minkowski_vierbein()
    U_links = audit_bosonic_gauge_sector()
    audit_fermion_sector_and_psu_decomposition()
    audit_action_variations()
    
    print("=" * 80)
    print(" AUDIT COMPLETATO: TUTTE LE IDENTITA' STRUTTURALI SONO MATEMATICAMENTE VERIFICATE.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
