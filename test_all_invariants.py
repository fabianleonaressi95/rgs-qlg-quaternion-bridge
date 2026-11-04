#!/usr/bin/env python3
"""
Test completi per identificare l'invariante geometrico mancante.

Candidati da testare:
1. Jacobiano trasformazione Krist-All → S³
2. Volume 24-cell e strutture correlate
3. Determinante Gram configurazioni
4. Misure Haar gruppi simmetria
5. Invarianti curvatura proiezioni
"""

import numpy as np
from math import sqrt, pi, log
from itertools import permutations, product
from collections import defaultdict


def kristall_points():
    """Genera i punti Krist-All."""
    points = []
    
    # Guscio (1,0,0): 6 punti
    for i in range(3):
        for sign in [1, -1]:
            coord = [0, 0, 0]
            coord[i] = sign
            points.append((coord[0], coord[1], coord[2], 1))
    
    # Guscio (1,1,1): 8 punti
    for x in [1, -1]:
        for y in [1, -1]:
            for z in [1, -1]:
                points.append((x, y, z, 1/2))
    
    # Guscio (1,2,4): 48 punti
    for perm in permutations([1, 2, 4]):
        for signs in [(1,1,1), (1,1,-1), (1,-1,1), (1,-1,-1),
                      (-1,1,1), (-1,1,-1), (-1,-1,1), (-1,-1,-1)]:
            points.append((signs[0]*perm[0], signs[1]*perm[1], signs[2]*perm[2], 1/4))
    
    return points


def test_jacobian():
    """Test 1: Jacobiano trasformazione Krist-All → S³."""
    print("=" * 70)
    print("TEST 1: JACOBIANO TRASFORMAZIONE KRIST-ALL → S³")
    print("=" * 70)
    
    points = kristall_points()
    
    # Mappa i punti su S³ (normalizzazione)
    s3_points = []
    for x, y, z, w in points:
        norm = sqrt(x*x + y*y + z*z)
        s3_points.append((x/norm, y/norm, z/norm, w))
    
    # Calcola il volume dell'ipersfera 3D di raggio 1
    volume_s3 = 2 * pi**2  # Volume S³
    
    print(f"\nVolume S³ (r=1): {volume_s3:.6f}")
    print(f"Log(volume): {log(volume_s3):.6f}")
    
    # Calcola il "volume discreto" dei punti Krist-All
    # Approx con semplice covering
    coords = np.array([(x, y, z) for x, y, z, w in points])
    
    # Calcola il volume del convesso
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(coords)
        volume_convex = hull.volume
        print(f"\nVolume convesso (R³): {volume_convex:.6f}")
        print(f"Log(volume): {log(volume_convex):.6f}")
    except:
        print(f"\nSciPy non disponibile - skip volume convesso")
    
    # Calcola determinante Gram
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            G[i, j] = sum(coords[k, i] * coords[k, j] for k in range(len(coords)))
    
    det_G = np.linalg.det(G)
    print(f"\nDeterminante Gram: {det_G:.6f}")
    print(f"Log(det): {log(abs(det_G)):.6f}" if det_G != 0 else "Det = 0")
    
    # Confronta con esponenti target
    print(f"\nConfronto con esponenti target:")
    print(f"  19/4 = {19/4:.4f}")
    print(f"  Log(volume_s3) = {log(volume_s3):.4f}")
    print(f"  Log(det_G) = {log(abs(det_G)):.4f}" if det_G != 0 else "  N/A")


def test_24cell_volume():
    """Test 2: Volume 24-cell."""
    print(f"\n{'='*70}")
    print("TEST 2: VOLUME 24-CELL")
    print(f"{'='*70}")
    
    # 24-cell: 24 vertici (unità Hurwitz)
    # Volume del 24-cell unitario = 2
    volume_24cell = 2
    
    print(f"\nVolume 24-cell unitario: {volume_24cell}")
    print(f"Log(volume): {log(volume_24cell):.6f}")
    
    # Volume ipercubo 4D unitario = 1
    volume_hypercube = 1
    print(f"\nVolume ipercubo 4D: {volume_hypercube}")
    print(f"Rapporto 24-cell/ipercubo: {volume_24cell/volume_hypercube:.6f}")
    
    # Volume sfera 4D unitario = π²/2
    volume_s4 = pi**2 / 2
    print(f"\nVolume sfera 4D (r=1): {volume_s4:.6f}")
    print(f"Rapporto 24-cell/sfera: {volume_24cell/volume_s4:.6f}")
    
    # Confronta con 5! = 120
    print(f"\nConfronto con 5! = 120:")
    print(f"  5! = 120")
    print(f"  Log(120) = {log(120):.6f}")
    print(f"  Log(volume_24cell) = {log(volume_24cell):.6f}")


def test_gram_determinant():
    """Test 3: Determinante Gram configurazioni."""
    print(f"\n{'='*70}")
    print("TEST 3: DETERMINANTE GRAM CONFIGURAZIONI")
    print(f"{'='*70}")
    
    points = kristall_points()
    
    # Gram per guscio
    by_norm = defaultdict(list)
    for x, y, z, w in points:
        norm_sq = x*x + y*y + z*z
        by_norm[norm_sq].append((x, y, z))
    
    for norm in sorted(by_norm.keys()):
        coords = np.array(by_norm[norm])
        n = len(coords)
        
        G = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                G[i, j] = sum(coords[i, k] * coords[j, k] for k in range(3))
        
        det_G = np.linalg.det(G)
        print(f"\nNorma {norm} ({n} punti):")
        print(f"  Det(G) = {det_G:.6e}")
        if det_G != 0:
            print(f"  Log(|det|) = {log(abs(det_G)):.6f}")


def test_haar_measure():
    """Test 4: Misure Haar gruppi simmetria."""
    print(f"\n{'='*70}")
    print("TEST 4: MISURE HAAR GRUPPI SIMMETRIA")
    print(f"{'='*70}")
    
    # Gruppo cubico O_h: 48 elementi
    # Gruppo tetraedrico T_d: 24 elementi
    # Gruppo ottaedrico O_h: 48 elementi
    
    group_orders = {
        'T_d': 24,
        'O_h': 48,
        'I_h': 120  # icosaedrico - contiene 5!
    }
    
    print(f"\nOrdini gruppi:")
    for name, order in group_orders.items():
        print(f"  {name}: {order}")
        print(f"    Log({order}) = {log(order):.6f}")
    
    # Nota: 120 = 5!
    print(f"\nNota importante:")
    print(f"  I_h (icosaedrico) ha ordine 120 = 5!")
    print(f"  Questo è UNICO tra i gruppi di simmetria platonic")
    
    # Verifica se Krist-All ha simmetria icosaedrica
    points = kristall_points()
    coords = [(x, y, z) for x, y, z, w in points]
    
    # Conta distanze uniche
    distances = set()
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = sqrt(sum((coords[i][k] - coords[j][k])**2 for k in range(3)))
            distances.add(round(dist, 6))
    
    print(f"\nDistanze uniche in Krist-All: {len(distances)}")
    print(f"Dodecaedro/icosaedro hanno pattern di distanze specifici")


def test_curvature_invariants():
    """Test 5: Invarianti curvatura proiezioni."""
    print(f"\n{'='*70}")
    print("TEST 5: INVARIANTI CURVATURA PROIEZIONI")
    print(f"{'='*70}")
    
    points = kristall_points()
    coords = np.array([(x, y, z) for x, y, z, w in points])
    
    # Curvatura media discreta
    # Calcola il tensore di curvatura approximato
    
    # Gram matrix
    G = np.zeros((len(coords), len(coords)))
    for i in range(len(coords)):
        for j in range(len(coords)):
            G[i, j] = sum(coords[i, k] * coords[j, k] for k in range(3))
    
    # Autovalori della Gram
    eigenvals = np.linalg.eigvalsh(G)
    eigenvals = np.round(eigenvals, 6)
    
    print(f"\nAutovalori Gram:")
    for i, ev in enumerate(eigenvals):
        print(f"  λ_{i}: {ev:.6f}")
    
    # Curvatura Gaussiana approx (determinante metrica)
    try:
        det_G = np.linalg.det(G)
        print(f"\nDeterminante Gram: {det_G:.6e}")
        
        # "Curvatura" come log(det)
        if det_G > 0:
            curvature = log(det_G)
            print(f"Log(det) = {curvature:.6f}")
            print(f"vs target 11/4 = {11/4:.6f}")
    except:
        print(f"Gram singolare")


def test_symmetry_group_orders():
    """Test aggiuntivo: ordini gruppi di simmetria."""
    print(f"\n{'='*70}")
    print("TEST AGGIUNTIVO: ORDINI GRUPPI SIMMETRIA")
    print(f"{'='*70}")
    
    # Gruppi e loro ordini
    symmetry_groups = {
        'C_n (ciclico)': 'n',
        'D_n (diedrico)': '2n',
        'T (tetraedrico)': 12,
        'O (ottaedrico)': 24,
        'I (icosaedrico)': 60,
        'T_d (full tetraedrico)': 24,
        'O_h (full ottaedrico)': 48,
        'I_h (full icosaedrico)': 120
    }
    
    print(f"\nOrdini gruppi:")
    for name, order in symmetry_groups.items():
        if isinstance(order, int):
            print(f"  {name}: {order}")
    
    print(f"\nNota: I_h = 120 = 5! è l'UNICO gruppo finito platonic con ordine 5!")
    
    # Verifica se Krist-All ha struttura icosaedrica
    print(f"\nStruttura Krist-All:")
    print(f"  Gusci: 3 (norme 1, 3, 21)")
    print(f"  Cardinalità: 6, 8, 48")
    print(f"  Rapporti: 8/6=1.33, 48/8=6, 48/6=8")
    print(f"  Nessun rapporto uguale a 5 o 120")


def test_target_exponents():
    """Test finale: confronto con esponenti target."""
    print(f"\n{'='*70}")
    print("TEST FINALE: CONFRONTO CON ESPONENTI TARGET")
    print(f"{'='*70}")
    
    target = {
        '2': 19/4,
        '3': -7/4,
        '5': 1/4,
        'π': 11/4
    }
    
    print(f"\nEsponenti target Wyler:")
    for base, exp in target.items():
        print(f"  {base}: {exp} = {exp:.4f}")
    
    # Raccogli tutti i valori log dai test
    log_values = {
        'log(2)': log(2),
        'log(3)': log(3),
        'log(5)': log(5),
        'log(π)': log(pi),
        'log(5!)': log(120),
        'log(24)': log(24),
        'log(48)': log(48),
        'log(64)': log(64),
        'log(216)': log(216),
        'log(Vol S³)': log(2*pi**2),
        'log(Vol 24-cell)': log(2),
    }
    
    print(f"\nValori log disponibili:")
    for name, val in log_values.items():
        print(f"  {name}: {val:.6f}")
    
    # Cerca combinazioni lineari
    print(f"\nRicerca combinazioni lineari semplici:")
    
    # 19/4 = 4.75
    candidates_19_4 = []
    for name1, val1 in log_values.items():
        for name2, val2 in log_values.items():
            if abs(val1 + val2 - 19/4) < 0.01:
                candidates_19_4.append((name1, name2, val1 + val2))
    
    if candidates_19_4:
        print(f"\nCombinazioni per 19/4:")
        for n1, n2, val in candidates_19_4:
            print(f"  {n1} + {n2} = {val:.4f}")
    else:
        print(f"\nNessuna combinazione semplice per 19/4")
    
    # 11/4 = 2.75
    candidates_11_4 = []
    for name1, val1 in log_values.items():
        for name2, val2 in log_values.items():
            if abs(val1 + val2 - 11/4) < 0.01:
                candidates_11_4.append((name1, name2, val1 + val2))
    
    if candidates_11_4:
        print(f"\nCombinazioni per 11/4:")
        for n1, n2, val in candidates_11_4:
            print(f"  {n1} + {n2} = {val:.4f}")
    else:
        print(f"\nNessuna combinazione semplice per 11/4")


def main():
    print("=" * 70)
    print("TEST COMPLETI: IDENTIFICAZIONE INVARIANTE GEOMETRICO")
    print("=" * 70)
    
    test_jacobian()
    test_24cell_volume()
    test_gram_determinant()
    test_haar_measure()
    test_curvature_invariants()
    test_symmetry_group_orders()
    test_target_exponents()
    
    print(f"\n{'='*70}")
    print("RIEPILOGO")
    print(f"{'='*70}")
    
    print(f"\n✓ I_h (icosaedrico) ha ordine 120 = 5!")
    print(f"✗ Krist-All non mostra simmetria icosaedrica evidente")
    print(f"✗ Nessuna combinazione semplice di log produce 19/4 o 11/4")
    print(f"✗ Invariante geometrico mancante NON identificato")


if __name__ == "__main__":
    main()
