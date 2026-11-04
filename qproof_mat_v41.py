#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MORPHODYNAMIC RG — STARLINGS / ROMA
===================================

Analisi:
  1. caricamento dati_uccelli.txt
  2. costruzione dei blocchi kNN
  3. tensore di forma S
  4. autovalori lambda1 >= lambda2 >= lambda3
  5. R12 = lambda1/lambda2
  6. R23 = lambda2/lambda3
  7. anisotropia Ab
  8. null model spaziale
  9. Z-score
 10. bootstrap
 11. ricerca del minimo rispetto a phi
 12. correzione look-elsewhere
 13. stima di un eventuale plateau / fixed-point morfologico

IMPORTANTE:
Il valore phi NON viene usato per costruire il modello.
Viene introdotto soltanto nel test finale.
"""

import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import bootstrap


# ============================================================
# CONFIGURAZIONE
# ============================================================

FILE = "dati_uccelli.txt"

SEED = 42

# Scale di coarse-graining
K_VALUES = np.array([
    4, 6, 8, 10, 12, 16, 20, 24,
    32, 40, 48, 56, 64, 80, 96, 112, 128
])

# Numero randomizzazioni null model
N_NULL = 2000

# Bootstrap
N_BOOT = 2000

# Selezione bordo:
# un blocco viene considerato "bulk" se il suo centro
# non è troppo vicino alla superficie stimata.
USE_BULK_FILTER = True

# soglia relativa alla scala locale
EDGE_FACTOR = 1.0

# se True, il null model preserva la nube globale
# attraverso una permutazione delle coordinate
NULL_MODE = "shuffle"

# Golden ratio — usato SOLO nel test finale
PHI = (1.0 + np.sqrt(5.0)) / 2.0


# ============================================================
# RANDOM GENERATOR
# ============================================================

rng = np.random.default_rng(SEED)


# ============================================================
# CARICAMENTO DATI
# ============================================================

def load_data(filename):
    """
    Cerca automaticamente colonne x,y,z.

    Formati supportati:

        x y z

        id x y z

        id frame x y z

        frame id x y z

    Le ultime tre colonne numeriche vengono utilizzate
    come coordinate se non è possibile riconoscere
    esplicitamente le colonne.
    """

    if not os.path.exists(filename):
        raise FileNotFoundError(
            f"File non trovato: {filename}\n"
            "Metti dati_uccelli.txt nella stessa directory dello script."
        )

    rows = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            parts = line.replace(",", " ").split()

            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue

            if len(vals) < 3:
                continue

            rows.append(vals)

    if len(rows) == 0:
        raise ValueError("Nessun dato numerico trovato.")

    data = np.asarray(rows, dtype=float)

    print("=" * 72)
    print("CARICAMENTO DATI")
    print("=" * 72)
    print("File:", filename)
    print("Shape:", data.shape)

    # ultime 3 colonne
    xyz = data[:, -3:]

    print("Coordinate utilizzate:")
    print("x =", xyz[:, 0].min(), "→", xyz[:, 0].max())
    print("y =", xyz[:, 1].min(), "→", xyz[:, 1].max())
    print("z =", xyz[:, 2].min(), "→", xyz[:, 2].max())

    return xyz


# ============================================================
# TENSORE DI FORMA
# ============================================================

def shape_tensor(points):
    """
    Tensore di forma / gyration tensor:

        S_ab = 1/N sum_i
               (r_i,a - R_a)(r_i,b - R_b)

    Restituisce:
        centroid
        eigenvalues
        eigenvectors
    """

    center = np.mean(points, axis=0)

    X = points - center

    S = (X.T @ X) / len(points)

    eigvals, eigvecs = np.linalg.eigh(S)

    # ordine decrescente
    idx = np.argsort(eigvals)[::-1]

    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # protezione numerica
    eigvals = np.maximum(eigvals, 0.0)

    return center, eigvals, eigvecs


# ============================================================
# OSSERVABILI MORFOLOGICHE
# ============================================================

def morphology_from_eigenvalues(lam):
    """
    lambda1 >= lambda2 >= lambda3

    R12 = lambda1/lambda2
    R23 = lambda2/lambda3

    Ab = anisotropia normalizzata
    """

    l1, l2, l3 = lam

    eps = 1e-15

    R12 = l1 / max(l2, eps)
    R23 = l2 / max(l3, eps)

    denom = 2.0 * (l1 + l2 + l3)**2

    if denom > eps:
        Ab = (
            (l1 - l2)**2
            + (l2 - l3)**2
            + (l3 - l1)**2
        ) / denom
    else:
        Ab = 0.0

    return R12, R23, Ab


# ============================================================
# kNN BLOCS
# ============================================================

def knn_blocks(points, k):
    """
    Per ogni punto costruisce un blocco kNN.

    Restituisce gli indici dei blocchi.
    """

    tree = cKDTree(points)

    _, neighbors = tree.query(points, k=k)

    return neighbors


# ============================================================
# SCALE FISICHE
# ============================================================

def estimate_scale(points, k):
    """
    Stima l(k) come distanza media dal k-esimo vicino.
    """

    tree = cKDTree(points)

    distances, _ = tree.query(points, k=k)

    # ultima colonna = distanza k-esimo vicino
    d_k = distances[:, -1]

    return np.mean(d_k)


# ============================================================
# OSSERVABILI SU TUTTI I BLOCCHI
# ============================================================

def analyze_k(points, k):
    """
    Calcola la distribuzione dei parametri morfologici
    sui blocchi kNN.
    """

    neighbors = knn_blocks(points, k)

    R12_values = []
    R23_values = []
    Ab_values = []

    centers = []

    for idx in neighbors:

        block = points[idx]

        center, eigvals, eigvecs = shape_tensor(block)

        R12, R23, Ab = morphology_from_eigenvalues(eigvals)

        R12_values.append(R12)
        R23_values.append(R23)
        Ab_values.append(Ab)

        centers.append(center)

    R12_values = np.asarray(R12_values)
    R23_values = np.asarray(R23_values)
    Ab_values = np.asarray(Ab_values)

    return {
        "R12": R12_values,
        "R23": R23_values,
        "Ab": Ab_values,
        "centers": np.asarray(centers)
    }


# ============================================================
# NULL MODEL
# ============================================================

def null_randomize(points):
    """
    Null model spaziale.

    Mantiene esattamente la distribuzione marginale
    delle coordinate ma distrugge la struttura
    spaziale congiunta.

    Ogni asse viene permutato indipendentemente.
    """

    null_points = points.copy()

    for axis in range(3):
        null_points[:, axis] = rng.permutation(
            null_points[:, axis]
        )

    return null_points


# ============================================================
# Z SCORE
# ============================================================

def z_score(real, null_mean, null_std):

    if null_std <= 1e-15:
        return 0.0

    return (real - null_mean) / null_std


# ============================================================
# BOOTSTRAP
# ============================================================

def bootstrap_mean(values, n_boot=N_BOOT):

    values = np.asarray(values)

    n = len(values)

    if n < 2:
        return np.mean(values), np.nan, np.nan

    boot_means = np.empty(n_boot)

    for i in range(n_boot):

        sample = rng.choice(
            values,
            size=n,
            replace=True
        )

        boot_means[i] = np.mean(sample)

    lo = np.percentile(boot_means, 2.5)
    hi = np.percentile(boot_means, 97.5)

    return np.mean(values), lo, hi


# ============================================================
# ANALISI PRINCIPALE
# ============================================================

def run_analysis(points):

    results = []

    print()
    print("=" * 72)
    print("MORPHODYNAMIC COARSE-GRAINING")
    print("=" * 72)

    for k in K_VALUES:

        if k > len(points):
            continue

        print(f"\nAnalisi k = {k}")

        scale = estimate_scale(points, k)

        real = analyze_k(points, k)

        R12_real = np.mean(real["R12"])
        R23_real = np.mean(real["R23"])
        Ab_real = np.mean(real["Ab"])

        # ----------------------------------------------------
        # NULL ENSEMBLE
        # ----------------------------------------------------

        null_R12 = []
        null_R23 = []
        null_Ab = []

        for n in range(N_NULL):

            null_points = null_randomize(points)

            null = analyze_k(null_points, k)

            null_R12.append(np.mean(null["R12"]))
            null_R23.append(np.mean(null["R23"]))
            null_Ab.append(np.mean(null["Ab"]))

        null_R12 = np.asarray(null_R12)
        null_R23 = np.asarray(null_R23)
        null_Ab = np.asarray(null_Ab)

        # ----------------------------------------------------
        # Z SCORES
        # ----------------------------------------------------

        Z12 = z_score(
            R12_real,
            np.mean(null_R12),
            np.std(null_R12, ddof=1)
        )

        Z23 = z_score(
            R23_real,
            np.mean(null_R23),
            np.std(null_R23, ddof=1)
        )

        ZAb = z_score(
            Ab_real,
            np.mean(null_Ab),
            np.std(null_Ab, ddof=1)
        )

        # ----------------------------------------------------
        # BOOTSTRAP REAL
        # ----------------------------------------------------

        _, R12_lo, R12_hi = bootstrap_mean(
            real["R12"]
        )

        _, R23_lo, R23_hi = bootstrap_mean(
            real["R23"]
        )

        _, Ab_lo, Ab_hi = bootstrap_mean(
            real["Ab"]
        )

        results.append({

            "k": k,
            "l": scale,

            "R12": R12_real,
            "R12_null": np.mean(null_R12),
            "Z12": Z12,

            "R23": R23_real,
            "R23_null": np.mean(null_R23),
            "Z23": Z23,

            "Ab": Ab_real,
            "Ab_null": np.mean(null_Ab),
            "ZAb": ZAb,

            "R12_lo": R12_lo,
            "R12_hi": R12_hi,

            "R23_lo": R23_lo,
            "R23_hi": R23_hi,

            "Ab_lo": Ab_lo,
            "Ab_hi": Ab_hi,

            "null_R12": null_R12,
            "null_R23": null_R23,
            "null_Ab": null_Ab
        })

        print(
            f"k={k:3d} | "
            f"l={scale:.3f} | "
            f"R12={R12_real:.4f}/{np.mean(null_R12):.4f} | "
            f"Z12={Z12:7.2f} | "
            f"R23={R23_real:.4f}/{np.mean(null_R23):.4f} | "
            f"Z23={Z23:7.2f} | "
            f"Ab={Ab_real:.4f}/{np.mean(null_Ab):.4f} | "
            f"ZAb={ZAb:7.2f}"
        )

    return results


# ============================================================
# GOLDEN RATIO TEST
# ============================================================

def golden_ratio_test(results):

    print()
    print("=" * 72)
    print("GOLDEN RATIO TEST")
    print("=" * 72)

    ks = np.array([r["k"] for r in results])
    ls = np.array([r["l"] for r in results])

    R12 = np.array([r["R12"] for r in results])

    # distanza assoluta
    distances = np.abs(R12 - PHI)

    idx = np.argmin(distances)

    best_k = ks[idx]
    best_l = ls[idx]
    best_R = R12[idx]

    print()
    print(f"phi = {PHI:.12f}")
    print()
    print("Migliore coincidenza osservata:")
    print(f"k*      = {best_k}")
    print(f"l*      = {best_l:.6f}")
    print(f"R12*    = {best_R:.8f}")
    print(f"|R-phi| = {distances[idx]:.8f}")
    print(
        f"errore relativo = "
        f"{100.0 * distances[idx] / PHI:.4f}%"
    )

    # --------------------------------------------------------
    # LOOK ELSEWHERE TEST
    # --------------------------------------------------------

    null_min_distances = []

    for r in results:

        null_R = r["null_R12"]

        # Per ogni null realization calcoliamo la distanza
        # minima da phi sulla scala corrispondente.
        null_min_distances.append(
            np.abs(null_R - PHI)
        )

    # ATTENZIONE:
    # qui abbiamo bisogno della distribuzione globale
    # del minimo tra tutte le scale.
    #
    # Ricostruiamo correttamente il test usando la stessa
    # realizzazione nulla per tutte le scale sarebbe ideale.
    #
    # Con l'ensemble già calcolato, facciamo una stima
    # conservativa aggregando i minimi per scala.

    null_min = np.min(
        np.asarray(null_min_distances),
        axis=0
    )

    real_min = np.min(distances)

    p_lee = np.mean(
        null_min <= real_min
    )

    print()
    print("Look-Elsewhere test:")
    print(f"D_real = {real_min:.8f}")
    print(f"p_LookElsewhere ≈ {p_lee:.6f}")

    if p_lee < 0.01:
        print(
            ">>> Evidenza forte: la vicinanza a phi "
            "è rara nei null model."
        )
    elif p_lee < 0.05:
        print(
            ">>> Evidenza moderata: possibile struttura "
            "compatibile con phi."
        )
    else:
        print(
            ">>> Nessuna evidenza sufficiente per "
            "attribuire significato speciale a phi."
        )

    return {
        "best_k": best_k,
        "best_l": best_l,
        "best_R": best_R,
        "distance": real_min,
        "p_lee": p_lee
    }


# ============================================================
# IDENTIFICAZIONE SCALA MORFOLOGICA
# ============================================================

def find_morphological_scale(results):

    print()
    print("=" * 72)
    print("MORPHOLOGICAL SCALE")
    print("=" * 72)

    R12 = np.array([r["R12"] for r in results])
    Ab = np.array([r["Ab"] for r in results])

    ks = np.array([r["k"] for r in results])
    ls = np.array([r["l"] for r in results])

    i_R = np.argmin(R12)
    i_A = np.argmin(Ab)

    print(
        f"min R12 : k={ks[i_R]}, "
        f"l={ls[i_R]:.4f}, "
        f"R12={R12[i_R]:.6f}"
    )

    print(
        f"min Ab  : k={ks[i_A]}, "
        f"l={ls[i_A]:.4f}, "
        f"Ab={Ab[i_A]:.6f}"
    )

    if ks[i_R] == ks[i_A]:
        print()
        print(
            ">>> Le due osservabili identificano "
            "la stessa scala morfologica."
        )
    else:
        print()
        print(
            ">>> Le due osservabili non hanno "
            "lo stesso minimo."
        )


# ============================================================
# STIMA PLATEAU
# ============================================================

def plateau_analysis(results):

    print()
    print("=" * 72)
    print("PLATEAU ANALYSIS")
    print("=" * 72)

    ks = np.array([r["k"] for r in results])
    ls = np.array([r["l"] for r in results])
    R = np.array([r["R12"] for r in results])

    # Regione candidata intorno a R ~ phi
    mask = np.abs(R - PHI) < 0.05

    if np.sum(mask) < 2:

        print(
            "Nessun plateau sufficientemente vicino "
            "a phi con soglia |R12-phi|<0.05."
        )

        return

    print()
    print("Punti compatibili con |R12 - phi| < 0.05:")

    for k, l, r in zip(
        ks[mask],
        ls[mask],
        R[mask]
    ):
        print(
            f"k={k:3d} | l={l:.4f} | R12={r:.6f}"
        )

    R_sel = R[mask]

    print()
    print(
        f"Media plateau = {np.mean(R_sel):.6f}"
    )

    print(
        f"Std plateau   = {np.std(R_sel, ddof=1):.6f}"
    )

    print(
        f"Distanza da phi = "
        f"{abs(np.mean(R_sel)-PHI):.6f}"
    )


# ============================================================
# SALVATAGGIO RISULTATI
# ============================================================

def save_results(results, filename="morphodynamics_results.csv"):

    header = [
        "k",
        "l",
        "R12",
        "R12_null",
        "Z12",
        "R23",
        "R23_null",
        "Z23",
        "Ab",
        "Ab_null",
        "ZAb",
        "R12_lo",
        "R12_hi",
        "R23_lo",
        "R23_hi",
        "Ab_lo",
        "Ab_hi"
    ]

    with open(filename, "w", encoding="utf-8") as f:

        f.write(",".join(header) + "\n")

        for r in results:

            row = [
                r["k"],
                r["l"],
                r["R12"],
                r["R12_null"],
                r["Z12"],
                r["R23"],
                r["R23_null"],
                r["Z23"],
                r["Ab"],
                r["Ab_null"],
                r["ZAb"],
                r["R12_lo"],
                r["R12_hi"],
                r["R23_lo"],
                r["R23_hi"],
                r["Ab_lo"],
                r["Ab_hi"]
            ]

            f.write(
                ",".join(
                    f"{x:.10g}" if isinstance(x, float)
                    else str(x)
                    for x in row
                )
                + "\n"
            )

    print()
    print("Risultati salvati in:", filename)


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 72)
    print("MORPHODYNAMIC RG — ROMA STARLINGS")
    print("=" * 72)
    print()
    print(f"Golden ratio phi = {PHI:.12f}")
    print(f"Seed             = {SEED}")
    print(f"Null realizations= {N_NULL}")
    print(f"Bootstrap        = {N_BOOT}")

    points = load_data(FILE)

    if len(points) < max(K_VALUES):

        raise ValueError(
            f"Numero di punti insufficiente: "
            f"{len(points)} < {max(K_VALUES)}"
        )

    results = run_analysis(points)

    find_morphological_scale(results)

    plateau_analysis(results)

    phi_result = golden_ratio_test(results)

    save_results(results)

    print()
    print("=" * 72)
    print("VERDETTO")
    print("=" * 72)

    print()

    print(
        "Il test phi è stato effettuato EX-POST e "
        "non è stato utilizzato nel coarse-graining."
    )

    print(
        f"Migliore R12 = {phi_result['best_R']:.8f} "
        f"a k={phi_result['best_k']}, "
        f"l={phi_result['best_l']:.4f}"
    )

    print(
        f"Distanza da phi = "
        f"{phi_result['distance']:.8f}"
    )

    print(
        f"p look-elsewhere ≈ "
        f"{phi_result['p_lee']:.6f}"
    )

    print()
    print(
        "Interpretazione: il valore di phi può essere "
        "considerato significativo solo se il test "
        "look-elsewhere sui null model risulta significativo."
    )


if __name__ == "__main__":
    main()
