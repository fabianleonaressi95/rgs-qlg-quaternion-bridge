#!/usr/bin/env python3
"""
Audit sistematico della norma 316 nell'ordine indicato:
1. Fibre Hopf: distribuzione esatta delle 16 fibre
2. Gram delle 16 immagini
3. Gram separata del settore 48→16
4. Gram del settore 16→8
5. Verifica se gli 8 punti sono un ottaedro
6. Azione di O_h sui 64 punti
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
    """h(q) = q*i*conj(q) / |q|²"""
    i_quat = (0, 1, 0, 0)
    q_norm = normalize_quaternion(q)
    conj_q = quaternion_conjugate(q_norm)
    temp = quaternion_multiply(q_norm, i_quat)
    result = quaternion_multiply(temp, conj_q)
    return result[1:]  # Parte immaginaria


def round_point(p, decimals=6):
    return tuple(round(coord, decimals) for coord in p)


def main():
    print("=" * 70)
    print("AUDIT SISTEMATICO NORMA 316")
    print("=" * 70)
    
    quaternions = load_quaternions('norm_316_quaternions.csv')
    
    # Separa per gruppo
    group_216 = [i for i, q in enumerate(quaternions) if q[4] == '216']
    group_300 = [i for i, q in enumerate(quaternions) if q[4] == '300']
    
    print(f"\nGruppo 216: {len(group_216)} quaternioni")
    print(f"Gruppo 300: {len(group_300)} quaternioni")
    
    # ============================================================
    # TEST 1: Fibre Hopf - distribuzione esatta
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 1: FIBRE HOPF - DISTRIBUZIONE ESATTA")
    print(f"{'='*70}")
    
    # Calcola immagini Hopf per tutti i 64 quaternioni
    hopf_images = []
    for q in quaternions:
        image = hopf_map(q[:4])
        hopf_images.append(round_point(image))
    
    # Calcola fibre: per ogni immagine, quali quaternioni mappano ad essa
    fibers = defaultdict(list)
    for i, image in enumerate(hopf_images):
        fibers[image].append(i)
    
    print(f"\nTotale immagini distinte: {len(fibers)}")
    print(f"Distribuzione delle fibre:")
    
    fiber_sizes = sorted([len(fibers[image]) for image in fibers])
    size_distribution = defaultdict(int)
    for size in fiber_sizes:
        size_distribution[size] += 1
    
    for size, count in sorted(size_distribution.items()):
        print(f"  Fibre di dimensione {size}: {count} immagini")
    
    # Verifica se è 16×4
    if len(fibers) == 16 and all(len(f) == 4 for f in fibers.values()):
        print(f"\n✓ Distribuzione esatta: 16 fibre × 4 punti = 64")
    else:
        print(f"\n✗ Distribuzione non uniforme")
    
    # ============================================================
    # TEST 2: Gram matrix delle 16 immagini
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 2: GRAM MATRIX DELLE 16 IMMAGINI")
    print(f"{'='*70}")
    
    # Estrai le 16 immagini distinte
    unique_images = list(fibers.keys())
    
    # Calcola Gram matrix
    n = len(unique_images)
    G = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(unique_images[i][k] * unique_images[j][k] for k in range(3))
    
    print(f"\nGram matrix (16×16):")
    print(G)
    
    # Trova valori distinti
    unique_values = sorted(set(round(G[i, j], 6) for i in range(n) for j in range(n)))
    print(f"\nValori distinti nella Gram: {unique_values}")
    
    # Calcola autovalori
    eigenvalues = np.linalg.eigvalsh(G)
    eigenvalues = np.round(eigenvalues, 6)
    print(f"\nAutovalori della Gram: {eigenvalues}")
    
    # ============================================================
    # TEST 3: Gram separata settore 48→16
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 3: GRAM SETTORE 48→16")
    print(f"{'='*70}")
    
    # Immagini dal gruppo 216
    images_216 = [hopf_images[i] for i in group_216]
    unique_216 = sorted(set(images_216))
    
    print(f"\nGruppo 216: {len(group_216)} quaternioni")
    print(f"Imagini distinte: {len(unique_216)}")
    
    # Gram delle immagini del gruppo 216
    n_216 = len(unique_216)
    G_216 = np.zeros((n_216, n_216))
    
    for i in range(n_216):
        for j in range(n_216):
            G_216[i, j] = sum(unique_216[i][k] * unique_216[j][k] for k in range(3))
    
    print(f"\nGram matrix settore 216 ({n_216}×{n_216}):")
    print(G_216)
    
    eigenvalues_216 = np.linalg.eigvalsh(G_216)
    eigenvalues_216 = np.round(eigenvalues_216, 6)
    print(f"\nAutovalori: {eigenvalues_216}")
    
    # ============================================================
    # TEST 4: Gram settore 16→8
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 4: GRAM SETTORE 16→8")
    print(f"{'='*70}")
    
    # Immagini dal gruppo 300
    images_300 = [hopf_images[i] for i in group_300]
    unique_300 = sorted(set(images_300))
    
    print(f"\nGruppo 300: {len(group_300)} quaternioni")
    print(f"Imagini distinte: {len(unique_300)}")
    
    # Gram delle immagini del gruppo 300
    n_300 = len(unique_300)
    G_300 = np.zeros((n_300, n_300))
    
    for i in range(n_300):
        for j in range(n_300):
            G_300[i, j] = sum(unique_300[i][k] * unique_300[j][k] for k in range(3))
    
    print(f"\nGram matrix settore 300 ({n_300}×{n_300}):")
    print(G_300)
    
    eigenvalues_300 = np.linalg.eigvalsh(G_300)
    eigenvalues_300 = np.round(eigenvalues_300, 6)
    print(f"\nAutovalori: {eigenvalues_300}")
    
    # ============================================================
    # TEST 5: Verifica se gli 8 punti sono un ottaedro
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 5: VERIFICA OTTAEDRO (SETTORE 300)")
    print(f"{'='*70}")
    
    print(f"\nCoordinate delle 8 immagini dal settore 300:")
    for i, p in enumerate(unique_300):
        print(f"  p_{i}: {p}")
    
    # Un ottaedro standard ha coordinate (±1, 0, 0) e permutazioni
    # Verifica se dopo normalizzazione e rotazione corrisponde
    
    # Calcola distanze tra punti
    print(f"\nDistanze tra punti:")
    distances = []
    for i in range(n_300):
        for j in range(i+1, n_300):
            dist = sqrt(sum((unique_300[i][k] - unique_300[j][k])**2 for k in range(3)))
            distances.append(dist)
    
    unique_distances = sorted(set(round(d, 6) for d in distances))
    print(f"  Distanze uniche: {unique_distances}")
    
    # In un ottaedro regolare, ci dovrebbero essere 2 distanze: spigolo e faccia
    if len(unique_distances) == 2:
        print(f"  ✓ Due distanze - compatibile con ottaedro")
        ratio = unique_distances[1] / unique_distances[0]
        print(f"  Rapporto: {ratio:.6f}")
        print(f"  In un ottaedro: √2 ≈ {sqrt(2):.6f}")
    else:
        print(f"  ✗ Non due distanze - non un ottaedro regolare")
    
    # ============================================================
    # TEST 6: Azione di O_h sui 64 punti
    # ============================================================
    print(f"\n{'='*70}")
    print("TEST 6: AZIONE DI O_h SUI 64 PUNTI")
    print(f"{'='*70}")
    
    # Verifica invarianza sotto permutazioni di coordinate e cambi di segno
    # O_h include 48 elementi: permutazioni (24) × cambi di segno (2)
    
    # Genera alcune trasformazioni O_h
    # Permutazioni di coordinate
    from itertools import permutations
    
    coords_4d = [q[:4] for q in quaternions]
    
    # Testa invarianza sotto cambio di segno
    sign_invariant = True
    for q in coords_4d:
        q_neg = (-q[0], -q[1], -q[2], -q[3])
        if round_point(q_neg) not in [round_point(q2) for q2 in coords_4d]:
            sign_invariant = False
            break
    
    print(f"\nInvarianza sotto cambio di segno globale: {'✓' if sign_invariant else '✗'}")
    
    # Testa permutazioni delle coordinate x,y,z (non t)
    coord_set = set(round_point(q) for q in coords_4d)
    perm_invariant = True
    
    for perm in permutations([1, 2, 3]):
        permuted_set = set()
        for q in coords_4d:
            permuted = (q[0], q[perm[0]], q[perm[1]], q[perm[2]])
            permuted_set.add(round_point(permuted))
        
        if permuted_set != coord_set:
            perm_invariant = False
            break
    
    print(f"Invarianza sotto permutazioni xyz: {'✓' if perm_invariant else '✗'}")
    
    # Verifica se i settori 216 e 300 sono invarianti separatamente
    coords_216 = [quaternions[i][:4] for i in group_216]
    coords_300 = [quaternions[i][:4] for i in group_300]
    
    set_216 = set(round_point(q) for q in coords_216)
    set_300 = set(round_point(q) for q in coords_300)
    
    # Testa invarianza settore 216
    perm_216_invariant = True
    for perm in permutations([1, 2, 3]):
        permuted_216 = set()
        for q in coords_216:
            permuted = (q[0], q[perm[0]], q[perm[1]], q[perm[2]])
            permuted_216.add(round_point(permuted))
        if permuted_216 != set_216:
            perm_216_invariant = False
            break
    
    print(f"Invarianza settore 216 sotto perm xyz: {'✓' if perm_216_invariant else '✗'}")
    
    # Testa invarianza settore 300
    perm_300_invariant = True
    for perm in permutations([1, 2, 3]):
        permuted_300 = set()
        for q in coords_300:
            permuted = (q[0], q[perm[0]], q[perm[1]], q[perm[2]])
            permuted_300.add(round_point(permuted))
        if permuted_300 != set_300:
            perm_300_invariant = False
            break
    
    print(f"Invarianza settore 300 sotto perm xyz: {'✓' if perm_300_invariant else '✗'}")
    
    # ============================================================
    # RIEPILOGO
    # ============================================================
    print(f"\n{'='*70}")
    print("RIEPILOGO AUDIT")
    print(f"{'='*70}")
    
    print(f"\n1. Fibre Hopf:")
    print(f"   64 → {len(fibers)} immagini")
    print(f"   Distribuzione: {dict(size_distribution)}")
    
    print(f"\n2. Gram 16 immagini:")
    print(f"   Autovalori: {eigenvalues}")
    
    print(f"\n3. Gram settore 216:")
    print(f"   {len(group_216)} → {len(unique_216)} immagini")
    print(f"   Autovalori: {eigenvalues_216}")
    
    print(f"\n4. Gram settore 300:")
    print(f"   {len(group_300)} → {len(unique_300)} immagini")
    print(f"   Autovalori: {eigenvalues_300}")
    
    print(f"\n5. Ottaedro settore 300:")
    print(f"   Distanze: {unique_distances}")
    
    print(f"\n6. Azione O_h:")
    print(f"   Invarianza settore 216: {'✓' if perm_216_invariant else '✗'}")
    print(f"   Invarianza settore 300: {'✓' if perm_300_invariant else '✗'}")


if __name__ == "__main__":
    main()
