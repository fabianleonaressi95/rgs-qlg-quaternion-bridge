import numpy as np
import scipy.linalg as la

def generate_aperiodic_chain(N, t=1.0, V=1.0, alpha=(np.sqrt(5) - 1) / 2):
    """
    Genera l'Hamiltoniana del sistema aperiodico (Modello di Aubry-André).
    """
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = V * np.cos(2 * np.pi * alpha * i)
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t
    H[0, N-1] = t
    H[N-1, 0] = t
    return H

def compute_filtered_gap_ratio(evals, fraction=0.4):
    """
    Calcola l'r-statistic isolando la regione centrale dello spettro 
    per evitare la distorsione causata dai gap di banda globali.
    """
    E = np.sort(evals)
    N = len(E)
    
    # Seleziona il 40% centrale dello spettro (evitando bordi e gap esterni)
    start_idx = int(N * (0.5 - fraction / 2))
    end_idx = int(N * (0.5 + fraction / 2))
    E_center = E[start_idx:end_idx]
    
    deltas = np.diff(E_center)
    r_values = []
    
    for i in range(len(deltas) - 1):
        d1 = deltas[i]
        d2 = deltas[i+1]
        # Filtra eventuali micro-gap residui anomali
        if max(d1, d2) > 1e-12:
            r = min(d1, d2) / max(d1, d2)
            if r <= 1.0:
                r_values.append(r)
                
    return np.mean(r_values) if r_values else 0.0

def run_rmt_analysis_corrected(N=2000):
    print("=" * 70)
    print(f" ANALISI RMT CORRETTA: LEVEL-SPACING CENTRALE (N = {N} nodi)")
    print("=" * 70)
    
    H = generate_aperiodic_chain(N)
    evals, _ = la.eigh(H)
    
    mean_r = compute_filtered_gap_ratio(evals, fraction=0.5)
    
    r_poisson = 2.0 * np.log(2) - 1.0  # ~0.3863
    r_goe = 0.5359                     # Wigner-Dyson
    
    print(f" <r> Corretto (Centro Spettro)     : {mean_r:.4f}")
    print(f" Riferimento Teorico - Poisson     : {r_poisson:.4f} (Integrabile / Regolare)")
    print(f" Riferimento Teorico - GOE         : {r_goe:.4f} (Caotico / Ergodico)")
    print("-" * 70)
    
    if abs(mean_r - r_poisson) < 0.08:
        print(" [ESITO SPETTRALE: ACCORDO CON REGIME POISSONIANO/CRITICO]")
        print(" Filtrando i gap di banda, lo spettro si allinea correttamente")
        print(" ai modelli aperiodici standard a transizione di mobilità.")
    elif mean_r > 0.45:
        print(" [ESITO SPETTRALE: REGIME CAOTICO / GOE]")
        print(" Lo spettro mostra forte repulsione dei livelli tipica dei sistemi ergodici.")
    else:
        print(" [ESITO SPETTRALE: STATISTICA INTERMEDIA / MULTIFRATTALE]")
        print(" Il valore riflette la natura multifrattale tipica del punto critico di Aubry-André.")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_rmt_analysis_corrected()
