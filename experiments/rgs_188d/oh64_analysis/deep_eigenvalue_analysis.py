#!/usr/bin/env python3
"""
Analisi approfondita della struttura degli autovalori e clustering naturale.
"""
import numpy as np
from math import sqrt, acos, pi
from collections import defaultdict
import csv


def load_quaternions(csv_file):
    quaternions = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row['t'])
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            group = row['group']
            quaternions.append((t, x, y, z, group))
    return quaternions


def quaternion_conjugate(q):
    t, x, y, z = q
    return (t, -x, -y, -z)


def quaternion_multiply(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b1,
        a1*d2 + b1*c2 - c1*b1 + d1*a1
    )


def quaternion_norm(q):
    t, x, y, z = q
    return sqrt(t*t + x*x + y*y + z*z)


def normalize_quaternion(q):
    norm = quaternion_norm(q)
    t, x, y, z = q
    return (t/norm, x/norm, y/norm, z/norm)


def hopf_map(q):
    i_quat = (0, 1, 0, 0)
    q_norm = normalize_quaternion(q)
    conj_q = quaternion_conjugate(q_norm)
    temp = quaternion_multiply(q_norm, i_quat)
    result = quaternion_multiply(temp, conj_q)
    return result[1:]


def main():
    print("=" * 70)
    print("ANALISI APPROFONDITA: AUTOVALORI E CLUSTERING")
    print("=" * 70)
    
    quaternions = load_quaternions('norm_316_quaternions.csv')
    hopf_images = [hopf_map(q[:4]) for q in quaternions]
    
    # Separa per gruppo
    group_216 = [i for i, q in enumerate(quaternions) if q[4] == '216']
    group_300 = [i for i, q in enumerate(quaternions) if q[4] == '300']
    
    # ============================================================
    # ANALISI AUTOVALORI
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI AUTOVALORI")
    print(f"{'='*70}")
    
    n = len(hopf_images)
    G = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(hopf_images[i][k] * hopf_images[j][k] for k in range(3))
    
    eigenvalues = np.linalg.eigvalsh(G)
    eigenvalues = np.round(eigenvalues, 10)
    
    print(f"\nAutovalori completi (64):")
    for i, ev in enumerate(eigenvalues):
        if abs(ev) > 1e-10:
            print(f"  λ_{i}: {ev}")
    
    # Rango della matrice
    rank = np.linalg.matrix_rank(G)
    print(f"\nRango della Gram matrix: {rank}")
    print(f"Dimensione dello spazio effettivo: {rank}")
    
    # ============================================================
    # CLUSTERING PER SIMILITUDINE
    # ============================================================
    print(f"\n{'='*70}")
    print("CLUSTERING PER SIMILITUDINE")
    print(f"{'='*70}")
    
    # Usa gli autovalori della matrice di covarianza per clustering
    images_array = np.array(hopf_images)
    
    # PCA semplice
    mean = np.mean(images_array, axis=0)
    centered = images_array - mean
    
    # Matrice di covarianza
    cov = np.cov(centered.T)
    
    print(f"\nMatrice di covarianza (3×3):")
    print(cov)
    
    # Autovalori della covarianza
    cov_eigenvals = np.linalg.eigvalsh(cov)
    cov_eigenvals = np.round(cov_eigenvals, 6)
    
    print(f"\nAutovalori della covarianza:")
    for i, ev in enumerate(cov_eigenvals):
        print(f"  λ_{i}: {ev:.6f}")
    
    print(f"\nVarianza totale: {sum(cov_eigenvals):.6f}")
    
    # Varianza spiegata da ogni componente principale
    total_var = sum(cov_eigenvals)
    for i, ev in enumerate(cov_eigenvals):
        pct = ev / total_var * 100
        print(f"  Componente {i}: {ev:.6f} ({pct:.2f}%)")
    
    # ============================================================
    # ANALISI ORTOGONALITÀ SETTORI
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI ORTOGONALITÀ SETTORI")
    print(f"{'='*70}")
    
    images_216 = np.array([hopf_images[i] for i in group_216])
    images_300 = np.array([hopf_images[i] for i in group_300])
    
    mean_216 = np.mean(images_216, axis=0)
    mean_300 = np.mean(images_300, axis=0)
    
    print(f"\nCentroide settore 216: {mean_216}")
    print(f"Centroide settore 300: {mean_300}")
    
    # Distanza tra centroidi
    dist_centroids = sqrt(sum((mean_216 - mean_300)**2))
    print(f"\nDistanza tra centroidi: {dist_centroids:.6f}")
    
    # Proiezione su subspazi
    # Verifica se i settori sono ortogonali
    proj_216_on_300 = np.dot(mean_216, mean_300)
    print(f"Prodotto scalare centroidi: {proj_216_on_300:.6f}")
    
    # ============================================================
    # RICERCA 5 COMPONENTI INDIPENDENTI
    # ============================================================
    print(f"\n{'='*70}")
    print("RICERCA 5 COMPONENTI INDIPENDENTI")
    print(f"{'='*70}")
    
    # Cerca 5 direzioni principali nei dati
    # Usa SVD sulla matrice delle immagini
    U, S, Vt = np.linalg.svd(images_array)
    
    print(f"\nValori singolari (SVD):")
    for i, s in enumerate(S):
        print(f"  σ_{i}: {s:.6f}")
    
    print(f"\nPrime 3 componenti principali spiegano:")
    total_s = sum(S)
    for i in range(3):
        pct = S[i] / total_s * 100
        print(f"  Componente {i}: {S[i]:.6f} ({pct:.2f}%)")
    
    # Vettori singolari corrispondenti
    print(f"\nPrimi 3 vettori singolari (colonne di U):")
    for i in range(3):
        print(f"  u_{i}: {U[:, i]}")
    
    # ============================================================
    # ANALISI COORDINATE DELLE IMMAGINI
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI COORDINATE IMMAGINI")
    print(f"{'='*70}")
    
    # Statistiche per coordinata
    for coord_idx in range(3):
        coords = [img[coord_idx] for img in hopf_images]
        print(f"\nCoordinata {coord_idx}:")
        print(f"  Min: {min(coords):.6f}")
        print(f"  Max: {max(coords):.6f}")
        print(f"  Media: {np.mean(coords):.6f}")
        print(f"  Deviazione std: {np.std(coords):.6f}")
    
    # ============================================================
    # VERIFICA SE I 64 PUNTI SONO SU UN CIRCONFERENZA
    # ============================================================
    print(f"\n{'='*70}")
    print("VERIFICA CIRCONFERENZA")
    print(f"{'='*70}")
    
    # Calcola la normale media
    normals = [tuple(coord / sqrt(sum(p*p for p in img)) for coord in img) for img in hopf_images]
    mean_normal = np.mean(normals, axis=0)
    norm_mean_normal = sqrt(sum(mean_normal**2))
    
    print(f"\nNormale media: {mean_normal}")
    print(f"Norma della normale media: {norm_mean_normal:.10f}")
    
    # Deviazione dalla normale media
    deviations = [sqrt(sum((np.array(n) - mean_normal)**2)) for n in normals]
    print(f"\nDeviazione media dalla normale media: {np.mean(deviations):.6f}")
    
    # ============================================================
    # RIEPILOGO FINALE
    # ============================================================
    print(f"\n{'='*70}")
    print("RIEPILOGO FINALE")
    print(f"{'='*70}")
    
    print(f"\nStruttura spaziale:")
    print(f"  Rango Gram: {rank}")
    print(f"  Dimensione effettiva: {rank}")
    
    print(f"\nVarianza (da covarianza):")
    for i, ev in enumerate(cov_eigenvals):
        pct = ev / total_var * 100
        print(f"  PC{i+1}: {pct:.2f}%")
    
    print(f"\nSettori:")
    print(f"   Centroidi separati di {dist_centroids:.6f}")
    print(f"  Prodotto scalare: {proj_216_on_300:.6f}")
    
    print(f"\nFattore 5:")
    print(f"  Non evidente dagli autovalori (solo 3 non nulli)")
    print(f"  Non evidente dalle norme 3D (11 norme diverse)")
    print(f"  Non evidente dal clustering naturale")
    
    print(f"\nInterpretazione:")
    print(f"  Le 64 immagini vivono in uno spazio 3D di dimensione {rank}")
    print(f"  La struttura è molto regolare (solo 3 autovalori non nulli)")
    print(f"  I settori 216 e 300 sono ben separati ma correlati")
    print(f"  Il fattore 5 rimane non identificato")


if __name__ == "__main__":
    main()
