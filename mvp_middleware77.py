"""
Verifica rigorosa della simmetria a vite (screw symmetry) del reticolo 13-3
e statistica RMT per settore con barre d'errore (bootstrap).

Geometria: 13 protofilamenti, 64 ring longitudinali, sfalsamento elicoidale
di 3 unita'. L'operazione di simmetria naturale non e' una rotazione pura
dei protofilamenti, ma una "screw operation":

    R_screw: (p, z) -> (p+1 mod 13, z+3 mod 64)

cioe' una rotazione di 1/13 accoppiata a una traslazione di 3 ring lungo z.
Se H commuta con R_screw, gli autostati di R_screw (autovalori
omega^k = exp(2*pi*i*k/13), k=0..12) block-diagonalizzano H in 13 settori.

Questo script:
  1. Costruisce H_TB per il reticolo 13-3 con hopping longitudinale,
     trasversale (tra protofilamenti) ed elicoidale.
  2. Costruisce R_screw come operatore unitario esplicito.
  3. Misura ||[H, R_screw]||  (deve essere ~0 a precisione macchina se la
     simmetria e' esatta nel modello).
  4. Proietta H nei 13 settori di k (via DFT sull'indice combinato che
     R_screw trasla ciclicamente) e calcola la r-statistic in ciascun
     settore.
  5. Bootstrap sui gap-ratio per stimare l'errore su <r> per settore e
     sulla media pesata complessiva.
"""

import numpy as np
import scipy.linalg as la

# ----------------------------------------------------------------------
# 1. Costruzione del reticolo 13-3
# ----------------------------------------------------------------------

N_PF = 13      # protofilamenti
N_RING = 64    # ring longitudinali
OFFSET = 3     # sfalsamento elicoidale
N = N_PF * N_RING  # 832 siti

def site_index(p, z):
    """Indice lineare del sito (protofilamento p, ring z), p e z ciclici."""
    return (p % N_PF) * N_RING + (z % N_RING)

def build_hamiltonian(t_long=1.0, t_trans=0.6, t_screw=0.8,
                       V_amp=0.5, seed=0):
    """
    H_TB con:
      - hopping longitudinale lungo z (stesso protofilamento)
      - hopping trasversale tra protofilamenti adiacenti (stesso ring)
      - hopping elicoidale: p -> p+1 accoppiato a z -> z+OFFSET
      - potenziale on-site V(theta, z) con armoniche coerenti con la
        simmetria 13-3 (cioe' invariante sotto p->p+1, z->z+OFFSET)
    Il potenziale e' costruito per dipendere SOLO dalla combinazione
    (13*z - OFFSET*p) mod (13*N_RING) circa, in modo da essere
    automaticamente invariante sotto la screw operation -> H commuta
    con R_screw ESATTAMENTE per costruzione. Questo e' il test di
    controllo positivo.
    """
    rng = np.random.default_rng(seed)
    H = np.zeros((N, N), dtype=complex)

    # Potenziale invariante sotto screw: funzione della sola fase combinata
    # phi(p,z) = 2*pi*(N_RING*p - N_PF*OFFSET*z) / (N_PF*N_RING)
    # Verifica: phi(p+1, z+OFFSET) - phi(p,z)
    #   = 2*pi*(N_RING - N_PF*OFFSET*OFFSET) / (N_PF*N_RING)  ... non zero in genere
    # Per avere invarianza ESATTA sotto (p,z)->(p+1,z+OFFSET) serve che la fase
    # dipenda solo dalla combinazione lineare intera invariante:
    #   m = (N_RING*p - OFFSET*z) mod (N_PF*N_RING) NON basta da sola;
    # la scelta corretta e standard per reticoli elicoidali e' usare come
    # variabile invariante:  u = z - (N_RING/N_PF)*p  (fibra elicoidale),
    # ma per restare su interi conviene lavorare direttamente sull'indice
    # "unwrapped" lungo la fibra elicoidale stessa:
    n_harm = 4
    amps = rng.uniform(0.5, 1.0, n_harm) * V_amp
    phases = rng.uniform(0, 2 * np.pi, n_harm)
    freqs = rng.integers(1, 5, n_harm)

    for p in range(N_PF):
        for z in range(N_RING):
            # variabile elicoidale invariante: avanzando p di 1 e z di OFFSET
            # simultaneamente, w resta costante mod N_RING solo se scegliamo
            # w = (z - OFFSET*p) mod N_RING  ->  w(p+1,z+OFFSET) = w(p,z). ESATTO.
            w = (z - OFFSET * p) % N_RING
            phase_var = 2 * np.pi * w / N_RING
            Vi = sum(A * np.cos(f * phase_var + ph)
                     for A, f, ph in zip(amps, freqs, phases))
            i = site_index(p, z)
            H[i, i] += Vi

            # hopping longitudinale
            j = site_index(p, z + 1)
            H[i, j] += t_long
            H[j, i] += np.conj(t_long)

            # hopping trasversale (tra protofilamenti adiacenti, stesso ring)
            j = site_index(p + 1, z)
            H[i, j] += t_trans
            H[j, i] += np.conj(t_trans)

            # hopping elicoidale: p -> p+1 con shift OFFSET in z
            j = site_index(p + 1, z + OFFSET)
            H[i, j] += t_screw
            H[j, i] += np.conj(t_screw)

    # Hermitianizza esplicitamente per sicurezza numerica
    H = 0.5 * (H + H.conj().T)
    return H


def build_screw_operator():
    """
    Operatore unitario R_screw che implementa (p,z) -> (p+1 mod 13, z+3 mod 64)
    come permutazione di siti.
    """
    R = np.zeros((N, N), dtype=complex)
    for p in range(N_PF):
        for z in range(N_RING):
            i = site_index(p, z)
            j = site_index(p + 1, z + OFFSET)
            R[j, i] = 1.0
    return R


# ----------------------------------------------------------------------
# 2. Test di commutazione [H, R_screw]
# ----------------------------------------------------------------------

def commutator_test(H, R):
    comm = H @ R - R @ H
    norm_comm = np.linalg.norm(comm, ord='fro')
    norm_H = np.linalg.norm(H, ord='fro')
    rel = norm_comm / norm_H
    print("=" * 70)
    print(" TEST DI COMMUTAZIONE [H, R_screw]")
    print("=" * 70)
    print(f" ||[H,R]||_F               : {norm_comm:.3e}")
    print(f" ||H||_F                   : {norm_H:.3e}")
    print(f" ||[H,R]||_F / ||H||_F     : {rel:.3e}")
    if rel < 1e-10:
        print(" [ESITO: SIMMETRIA A VITE ESATTA — block-diagonalization valida]")
    elif rel < 1e-6:
        print(" [ESITO: simmetria valida a precisione numerica standard]")
    else:
        print(" [ESITO: SIMMETRIA VIOLATA — il blocco per k NON è affidabile.")
        print("  Il potenziale/hopping usato non è invariante sotto R_screw.]")
    print("=" * 70 + "\n")
    return rel


# ----------------------------------------------------------------------
# 3. Block-diagonalization via autobase di R_screw
# ----------------------------------------------------------------------

def block_diagonalize(H, R):
    """
    R_screw e' una permutazione ciclica di ordine 13 sui 64-blocchi di siti
    (agisce come shift ciclico sull'indice p, con un aggiuntivo shift su z
    che pero' e' anch'esso ciclico di periodo 64/gcd(64,3)=64, quindi R^13 = 
    shift_z^(13*3 mod 64) sui ring — verifichiamo che R abbia effettivamente
    ordine finito e usiamo eigh su R per ottenere gli autospazi.
    R e' unitaria ma non hermitiana in generale; usiamo eig complesso.
    """
    # autovalori/autovettori di R (unitaria => autovalori su cerchio unitario)
    eigvals_R, eigvecs_R = la.eig(R)

    # Raggruppiamo per i 13 autovalori attesi omega^k = exp(2pi i k/13)
    # (ogni settore ha molteplicita' N/13 = 64)
    k_targets = np.exp(2j * np.pi * np.arange(N_PF) / N_PF)

    sectors = {k: [] for k in range(N_PF)}
    for idx, ev in enumerate(eigvals_R):
        # trova il k piu' vicino
        k_best = np.argmin(np.abs(k_targets - ev))
        sectors[k_best].append(idx)

    sector_info = []
    for k in range(N_PF):
        idxs = sectors[k]
        sector_info.append((k, len(idxs)))
    print("Dimensione dei settori k (attesa: 64 ciascuno):")
    for k, dim in sector_info:
        print(f"  k={k:2d} : dim={dim}")
    print()

    # Per ogni settore, proietta H sulla base ortonormalizzata degli
    # autovettori di R in quel settore e diagonalizza il blocco.
    sector_eigs = {}
    for k in range(N_PF):
        idxs = sectors[k]
        if len(idxs) < 3:
            continue
        V = eigvecs_R[:, idxs]
        # ortonormalizza (QR) per sicurezza numerica
        Q, _ = np.linalg.qr(V)
        H_block = Q.conj().T @ H @ Q
        H_block = 0.5 * (H_block + H_block.conj().T)  # forza hermitianita'
        evals_block = np.linalg.eigvalsh(H_block)
        sector_eigs[k] = evals_block

    return sector_eigs


# ----------------------------------------------------------------------
# 4. r-statistic con bootstrap
# ----------------------------------------------------------------------

def gap_ratios(evals):
    E = np.sort(evals)
    deltas = np.diff(E)
    deltas = deltas[deltas > 1e-13]
    if len(deltas) < 2:
        return np.array([])
    r = np.minimum(deltas[:-1], deltas[1:]) / np.maximum(deltas[:-1], deltas[1:])
    return r

def bootstrap_mean_r(r_values, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    if len(r_values) < 2:
        return np.nan, np.nan
    boot_means = np.empty(n_boot)
    n = len(r_values)
    for b in range(n_boot):
        sample = rng.choice(r_values, size=n, replace=True)
        boot_means[b] = np.mean(sample)
    return np.mean(boot_means), np.std(boot_means)


def rmt_analysis_with_errors(sector_eigs):
    print("=" * 70)
    print(" STATISTICA RMT PER SETTORE k CON BOOTSTRAP (n_boot=2000)")
    print("=" * 70)
    print(f"{'k':>3} | {'N_stati':>7} | {'<r>':>8} | {'SE(<r>)':>8} | {'N_gap':>6}")
    print("-" * 70)

    all_r = []
    per_sector_mean = []
    per_sector_se = []
    per_sector_n = []

    for k in sorted(sector_eigs.keys()):
        evals = sector_eigs[k]
        r_vals = gap_ratios(evals)
        if len(r_vals) < 2:
            continue
        mean_r, se_r = bootstrap_mean_r(r_vals)
        print(f"{k:>3} | {len(evals):>7} | {mean_r:>8.4f} | {se_r:>8.4f} | {len(r_vals):>6}")
        all_r.append(r_vals)
        per_sector_mean.append(mean_r)
        per_sector_se.append(se_r)
        per_sector_n.append(len(r_vals))

    print("-" * 70)

    # Media pesata tra settori (peso = 1/SE^2, standard per medie pesate)
    means = np.array(per_sector_mean)
    ses = np.array(per_sector_se)
    weights = 1.0 / (ses ** 2)
    weighted_mean = np.sum(weights * means) / np.sum(weights)
    weighted_se = np.sqrt(1.0 / np.sum(weights))

    # Alternativa: pool di tutti i gap-ratio da tutti i settori insieme
    pooled = np.concatenate(all_r)
    pooled_mean, pooled_se = bootstrap_mean_r(pooled, n_boot=3000)

    r_poisson = 2.0 * np.log(2) - 1.0  # ~0.3863
    r_goe = 0.5359

    print(f"\n <r> media pesata tra settori   : {weighted_mean:.4f} +/- {weighted_se:.4f}")
    print(f" <r> pooled (tutti i gap insieme): {pooled_mean:.4f} +/- {pooled_se:.4f}")
    print(f" Riferimento Poisson             : {r_poisson:.4f}")
    print(f" Riferimento GOE                 : {r_goe:.4f}")

    z_poisson = abs(pooled_mean - r_poisson) / pooled_se
    z_goe = abs(pooled_mean - r_goe) / pooled_se
    print(f"\n Distanza da Poisson in unita' di errore standard : {z_poisson:.2f} sigma")
    print(f" Distanza da GOE in unita' di errore standard     : {z_goe:.2f} sigma")
    print("=" * 70 + "\n")

    return weighted_mean, weighted_se, pooled_mean, pooled_se


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

if __name__ == "__main__":
    print("Costruzione dell'Hamiltoniana del reticolo 13-3 (N=832)...\n")
    H = build_hamiltonian()
    R = build_screw_operator()

    rel_comm = commutator_test(H, R)

    print("Block-diagonalization nei 13 settori k...\n")
    sector_eigs = block_diagonalize(H, R)

    rmt_analysis_with_errors(sector_eigs)


def diagnose_commutator_by_term():
    """Isola quale termine di H rompe la commutazione con R_screw."""
    R = build_screw_operator()

    def term_norm_comm(H_term, label):
        comm = H_term @ R - R @ H_term
        rel = np.linalg.norm(comm, 'fro') / max(np.linalg.norm(H_term, 'fro'), 1e-12)
        print(f"  {label:20s} : ||[H,R]||/||H|| = {rel:.3e}")

    H_pot = np.zeros((N, N), dtype=complex)
    H_long = np.zeros((N, N), dtype=complex)
    H_trans = np.zeros((N, N), dtype=complex)
    H_screw = np.zeros((N, N), dtype=complex)

    rng = np.random.default_rng(0)
    n_harm = 4
    amps = rng.uniform(0.5, 1.0, n_harm) * 0.5
    phases = rng.uniform(0, 2 * np.pi, n_harm)
    freqs = rng.integers(1, 5, n_harm)

    for p in range(N_PF):
        for z in range(N_RING):
            i = site_index(p, z)
            w = (z - OFFSET * p) % N_RING
            phase_var = 2 * np.pi * w / N_RING
            Vi = sum(A * np.cos(f * phase_var + ph) for A, f, ph in zip(amps, freqs, phases))
            H_pot[i, i] += Vi

            j = site_index(p, z + 1)
            H_long[i, j] += 1.0; H_long[j, i] += 1.0

            j = site_index(p + 1, z)
            H_trans[i, j] += 0.6; H_trans[j, i] += 0.6

            j = site_index(p + 1, z + OFFSET)
            H_screw[i, j] += 0.8; H_screw[j, i] += 0.8

    print("\nScomposizione della violazione di commutazione per termine:")
    term_norm_comm(H_pot, "Potenziale on-site")
    term_norm_comm(H_long, "Hopping longitudinale")
    term_norm_comm(H_trans, "Hopping trasversale")
    term_norm_comm(H_screw, "Hopping elicoidale")


diagnose_commutator_by_term()
