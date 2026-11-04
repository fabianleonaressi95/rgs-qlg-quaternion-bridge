#!/usr/bin/env python3
"""
Analisi dettagliata della Gram matrix delle 64 immagini Hopf.
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
        a1*d2 + b1*c2 - c1*b1 + d1*a2
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


def round_point(p, decimals=6):
    return tuple(round(coord, decimals) for coord in p)


def main():
    print("=" * 70)
    print("ANALISI DETTAGLIATA GRAM 64 IMMAGINI HOPF")
    print("=" * 70)
    
    quaternions = load_quaternions('norm_316_quaternions.csv')
    
    # Calcola immagini Hopf
    hopf_images = [hopf_map(q[:4]) for q in quaternions]
    
    # Separa per gruppo
    group_216 = [i for i, q in enumerate(quaternions) if q[4] == '216']
    group_300 = [i for i, q in enumerate(quaternions) if q[4] == '300']
    
    print(f"\nGruppo 216: {len(group_216)} quaternioni")
    print(f"Gruppo 300: {len(group_300)} quaternioni")
    
    # ============================================================
    # GRAM MATRIX COMPLETA 64×64
    # ============================================================
    print(f"\n{'='*70}")
    print("GRAM MATRIX COMPLETA 64×64")
    print(f"{'='*70}")
    
    n = len(hopf_images)
    G = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(hopf_images[i][k] * hopf_images[j][k] for k in range(3))
    
    # Autovalori
    eigenvalues = np.linalg.eigvalsh(G)
    eigenvalues = np.round(eigenvalues, 6)
    
    print(f"\nAutovalori della Gram (64×64):")
    for i, ev in enumerate(eigenvalues):
        print(f"  λ_{i}: {ev:.6f}")
    
    # Distribuzione degli autovalori
    unique_eigenvals, counts = np.unique(eigenvalues, return_counts=True)
    print(f"\nAutovalori unici e molteplicità:")
    for ev, count in zip(unique_eigenvals, counts):
        print(f"  {ev:.6f}: molteplicità {count}")
    
    # Valori della Gram matrix
    gram_values = [G[i, j] for i in range(n) for j in range(n)]
    unique_gram = sorted(set(round(v, 6) for v in gram_values))
    
    print(f"\nValori distinti nella Gram: {len(unique_gram)}")
    print(f"  Intervallo: [{min(unique_gram):.6f}, {max(unique_gram):.6f}]")
    
    # ============================================================
    # ANALISI PER SETTORE
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI PER SETTORE")
    print(f"{'='*70}")
    
    # Gram settore 216
    images_216 = [hopf_images[i] for i in group_216]
    n_216 = len(images_216)
    G_216 = np.zeros((n_216, n_216))
    
    for i in range(n_216):
        for j in range(n_216):
            G_216[i, j] = sum(images_216[i][k] * images_216[j][k] for k in range(3))
    
    eigenvalues_216 = np.linalg.eigvalsh(G_216)
    eigenvalues_216 = np.round(eigenvalues_216, 6)
    
    print(f"\nSettore 216 (48×48):")
    print(f"  Autovalori: {eigenvalues_216}")
    
    # Gram settore 300
    images_300 = [hopf_images[i] for i in group_300]
    n_300 = len(images_300)
    G_300 = np.zeros((n_300, n_300))
    
    for i in range(n_300):
        for j in range(n_300):
            G_300[i, j] = sum(images_300[i][k] * images_300[j][k] for k in range(3))
    
    eigenvalues_300 = np.linalg.eigvalsh(G_300)
    eigenvalues_300 = np.round(eigenvalues_300, 6)
    
    print(f"\nSettore 300 (16×16):")
    print(f"  Autovalori: {eigenvalues_300}")
    
    # ============================================================
    # ANALISI BLOCCO 64×64
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI STRUTTURA A BLOCCHI")
    print(f"{'='*70}")
    
    # Dividi la Gram in blocchi 48×48, 48×16, 16×48, 16×16
    G_11 = G[:48, :48]  # 216-216
    G_12 = G[:48, 48:]  # 216-300
    G_21 = G[48:, :48]  # 300-216
    G_22 = G[48:, 48:]  # 300-300
    
    print(f"\nBlocco G_11 (216-216):")
    print(f"  Forma: {G_11.shape}")
    print(f"  Somma elementi: {np.sum(G_11):.6f}")
    print(f"  Media: {np.mean(G_11):.6f}")
    
    print(f"\nBlocco G_12 (216-300):")
    print(f"  Forma: {G_12.shape}")
    print(f"  Somma elementi: {np.sum(G_12):.6f}")
    print(f"  Media: {np.mean(G_12):.6f}")
    
    print(f"\nBlocco G_22 (300-300):")
    print(f"  Forma: {G_22.shape}")
    print(f"  Somma elementi: {np.sum(G_22):.6f}")
    print(f"  Media: {np.mean(G_22):.6f}")
    
    # Verifica se G_12 è trasposta di G_21
    if np.allclose(G_12, G_21.T):
        print(f"\n✓ G_12 = G_21^T (simmetrica tra settori)")
    else:
        print(f"\n✗ G_12 ≠ G_21^T")
    
    # ============================================================
    # RICERCA DEL FATTORE 5
    # ============================================================
    print(f"\n{'='*70}")
    print("RICERCA DEL FATTORE 5")
    print(f"{'='*70}")
    
    # Cerca 5 classi nelle immagini
    # Raggruppa per norma 3D delle immagini
    norms_3d = [sum(p*p for p in img) for img in hopf_images]
    unique_norms = sorted(set(round(n, 6) for n in norms_3d))
    
    print(f"\nNorme 3D delle immagini: {unique_norms}")
    print(f"Numero di norme diverse: {len(unique_norms)}")
    
    # Cerca clustering naturale usando k-means semplice
    if len(unique_norms) >= 5:
        print(f"\nTentativo clustering semplice per norma 3D:")
        
        # Raggruppa per norma più vicina
        cluster_assignment = []
        for img in hopf_images:
            norm = sum(p*p for p in img)
            # Trova la norma più vicina
            closest_idx = min(range(len(unique_norms)), key=lambda i: abs(norm - unique_norms[i]))
            cluster_assignment.append(closest_idx)
        
        cluster_sizes = defaultdict(int)
        for c in cluster_assignment:
            cluster_sizes[c] += 1
        
        print(f"  Cluster per norma 3D: {dict(cluster_sizes)}")
    else:
        print(f"\nSolo {len(unique_norms)} norme diverse - non sufficienti per 5 cluster")
    
    # Cerca 5 autovalori significativi
    nonzero_eigenvals = [ev for ev in eigenvalues if abs(ev) > 1e-6]
    print(f"\nAutovalori non nulli: {len(nonzero_eigenvals)}")
    
    if len(nonzero_eigenvals) >= 5:
        print(f"Autovalori significativi ( primi 5): {nonzero_eigenvals[:5]}")
    else:
        print(f"Tutti gli autovalori non nulli: {nonzero_eigenvals}")
    
    # Cerca pattern nei prodotti scalari
    print(f"\nAnalisi pattern nei prodotti scalari:")
    
    # Prodotto scalare medio intra-settore
    mean_216 = np.mean(G_11)
    mean_300 = np.mean(G_22)
    mean_cross = np.mean(G_12)
    
    print(f"  Prodotto medio intra-216: {mean_216:.6f}")
    print(f"  Prodotto medio intra-300: {mean_300:.6f}")
    print(f"  Prodotto medio tra settori: {mean_cross:.6f}")
    
    # Rapporti
    if mean_300 > 1e-6:
        ratio = mean_216 / mean_300
        print(f"  Rapporto 216/300: {ratio:.6f}")
    
    # ============================================================
    # DISTANZE ANGOLARI
    # ============================================================
    print(f"\n{'='*70}")
    print("ANALISI DISTANZE ANGOLARI")
    print(f"{'='*70}")
    
    def angular_distance(p1, p2):
        dot = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
        dot = max(-1, min(1, dot))
        return acos(dot)
    
    # Calcola tutte le distanze
    distances = []
    for i in range(n):
        for j in range(i+1, n):
            dist = angular_distance(hopf_images[i], hopf_images[j])
            distances.append(dist)
    
    # Statistiche
    print(f"\nStatistiche distanze angolari:")
    print(f"  Minima: {min(distances):.6f} rad ({min(distances)*180/pi:.2f}°)")
    print(f"  Massima: {max(distances):.6f} rad ({max(distances)*180/pi:.2f}°)")
    print(f"  Media: {np.mean(distances):.6f} rad ({np.mean(distances)*180/pi:.2f}°)")
    
    # Istogramma
    hist_bins = [0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]
    hist_counts = [0] * (len(hist_bins) - 1)
    
    for d in distances:
        for i in range(len(hist_bins) - 1):
            if hist_bins[i] <= d < hist_bins[i+1]:
                hist_counts[i] += 1
                break
    
    print(f"\nIstogramma distanze:")
    for i in range(len(hist_counts)):
        print(f"  [{hist_bins[i]*180/pi:.0f}°, {hist_bins[i+1]*180/pi:.0f}°): {hist_counts[i]}")
    
    # ============================================================
    # RIEPILOGO
    # ============================================================
    print(f"\n{'='*70}")
    print("RIEPILOGO")
    print(f"{'='*70}")
    
    print(f"\nStruttura 64 immagini:")
    print(f"  64 immagini distinte su S²")
    print(f"  Autovalori Gram: {len(unique_eigenvals)} valori unici")
    
    print(f"\nDecomposizione 48+16:")
    print(f"  Settore 216: 48 immagini")
    print(f"  Settore 300: 16 immagini")
    print(f"  Correlazione tra settori: {'forte' if abs(mean_cross) > 0.1 else 'debole'}")
    
    print(f"\nFattore 5:")
    if len(unique_norms) == 5:
        print(f"  ✓ Trovate 5 norme 3D diverse - possibile collegamento")
    else:
        print(f"  Norme 3D trovate: {len(unique_norms)}")
        print(f"  Fattore 5 non evidente dalle norme 3D")


if __name__ == "__main__":
    main()
