#!/usr/bin/env python3
"""
Audit completo delle 4 norme con cardinalità 64 nel reticolo di Hurwitz.
Per N ∈ {76, 79, 84, 91}, calcola:
- |ℋ_N|
- Orbite sotto 2T, T, O_h
- Matrici Gram
- Proiezioni 3D
- Confronto con Krist-All e 64 Gate
"""

import numpy as np
from math import sqrt
from itertools import product, combinations
from collections import defaultdict


def hurwitz_integers(limit=10):
    """Genera interi di Hurwitz fino a un certo limite."""
    hurwitz = []
    # Tipo 1: tutti coefficienti interi
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    # Tipo 2: mezzi-interi
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    return hurwitz


def hurwitz_units():
    """Restituisce le 24 unità di Hurwitz (2T)."""
    units = []
    # 8 unità standard
    for sign in [1, -1]:
        units.append((sign, 0, 0, 0))
        units.append((0, sign, 0, 0))
        units.append((0, 0, sign, 0))
        units.append((0, 0, 0, sign))
    # 16 unità half
    for signs in product([1, -1], repeat=4):
        units.append((signs[0]/2, signs[1]/2, signs[2]/2, signs[3]/2))
    return units


def quaternion_norm(q):
    """Norma di un quaternione."""
    a, b, c, d = q
    return a*a + b*b + c*c + d*d


def quaternion_multiply(q1, q2):
    """Moltiplicazione quaternioni."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2
    )


def quaternion_conjugate(q):
    """Coniugato."""
    a, b, c, d = q
    return (a, -b, -c, -d)


def left_orbit(q, units):
    """Orbita sinistra sotto 2T."""
    return [quaternion_multiply(u, q) for u in units]


def right_orbit(q, units):
    """Orbita destra sotto 2T."""
    return [quaternion_multiply(q, u) for u in units]


def project_to_3d(q):
    """Proiezione su ImH ≃ R³."""
    a, b, c, d = q
    return (b, c, d)


def compute_gram_matrix(quaternions):
    """Calcola matrice Gram dei prodotti scalari."""
    n = len(quaternions)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a1, b1, c1, d1 = quaternions[i]
            a2, b2, c2, d2 = quaternions[j]
            G[i, j] = a1*a2 + b1*b2 + c1*c2 + d1*d2
    return G


def compute_gram_matrix_3d(vectors_3d):
    """Calcola matrice Gram delle proiezioni 3D."""
    n = len(vectors_3d)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(vectors_3d[i][k] * vectors_3d[j][k] for k in range(3))
    return G


def get_hurwitz_by_norm(hurwitz_set, target_norms):
    """Estrae quaternioni con norme specifiche."""
    by_norm = defaultdict(list)
    for q in hurwitz_set:
        t, x, y, z, qtype = q
        norm = t*t + x*x + y*y + z*z
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        if norm in target_norms:
            by_norm[norm].append((t, x, y, z, qtype))
    return by_norm


def analyze_orbit_structure(quaternions, units, norm):
    """Analizza struttura delle orbite per una data norma."""
    print(f"\n{'='*70}")
    print(f"ANALISI NORMA {norm}")
    print(f"{'='*70}")
    
    n_quats = len(quaternions)
    print(f"\nCardinalità totale: {n_quats}")
    
    if n_quats == 0:
        print("Nessun quaternione trovato.")
        return None
    
    # Calcola orbite sinistra per un campione
    print("\nAnalisi orbite sinistra (campione di 5):")
    orbit_sizes = defaultdict(int)
    
    for i, q in enumerate(quaternions[:5]):
        q_4d = q[:4]
        orbit = left_orbit(q_4d, units)
        orbit_set = set((round(a, 6), round(b, 6), round(c, 6), round(d, 6))
                       for a, b, c, d in orbit)
        orbit_size = len(orbit_set)
        orbit_sizes[orbit_size] += 1
        print(f"  Quaternione {i}: orbita dimensione {orbit_size}")
    
    print(f"\nDistribuzione dimensioni orbite:")
    for size, count in sorted(orbit_sizes.items()):
        print(f"  Dimensione {size}: {count} rappresentanti")
    
    # Analisi proiezioni 3D
    print("\nAnalisi proiezioni 3D:")
    projections_3d = [project_to_3d(q[:4]) for q in quaternions]
    
    # Calcola norme 3D
    norms_3d = [sum(x*x for x in proj) for proj in projections_3d]
    unique_norms_3d = sorted(set(round(n, 6) for n in norms_3d))
    
    print(f"  Norme 3D uniche: {unique_norms_3d}")
    
    # Matrice Gram 4D
    print("\nMatrice Gram 4D (prime 5×5):")
    G_4d = compute_gram_matrix([q[:4] for q in quaternions[:5]])
    print(G_4d)
    
    # Matrice Gram 3D
    print("\nMatrice Gram 3D (prime 5×5):")
    G_3d = compute_gram_matrix_3d(projections_3d[:5])
    print(G_3d)
    
    # Analisi distribuzione coordinate
    print("\nDistribuzione coordinate (valori unici):")
    for coord_idx in range(4):
        values = [q[coord_idx] for q in quaternions]
        unique_vals = sorted(set(round(v, 2) for v in values))
        print(f"  Coordinata {coord_idx}: {unique_vals[:10]}...")
    
    # Analisi tipi (integer vs half-integer)
    type_counts = defaultdict(int)
    for q in quaternions:
        type_counts[q[4]] += 1
    print(f"\nTipi di quaternioni:")
    for qtype, count in type_counts.items():
        print(f"  {qtype}: {count}")
    
    return {
        'cardinality': n_quats,
        'orbit_sizes': dict(orbit_sizes),
        'norms_3d': unique_norms_3d,
        'type_counts': dict(type_counts),
        'quaternions': quaternions,
        'projections_3d': projections_3d
    }


def compare_with_kristall(analysis, norm):
    """Confronta con struttura Krist-All."""
    print(f"\nConfronto con Krist-All (norma {norm}):")
    
    # Krist-All ha norme 1, 3, 21 in ImH
    kristall_norms = [1, 3, 21]
    
    # Se la norma corrisponde a una Krist-All
    if norm in kristall_norms:
        print(f"  ✓ Norma {norm} è una norma Krist-All")
    else:
        print(f"  ✗ Norma {norm} non è una norma Krist-All")
        print(f"    Norme Krist-All: {kristall_norms}")
    
    # Verifica se le proiezioni 3D hanno norme interessanti
    if analysis and 'norms_3d' in analysis:
        norms_3d = analysis['norms_3d']
        print(f"  Norme 3D trovate: {norms_3d}")
        
        # Cerca corrispondenze
        for kn in kristall_norms:
            if any(abs(n - kn) < 0.01 for n in norms_3d):
                print(f"    ✓ Corrispondenza con norma Krist-All {kn}")


def search_for_41(hurwitz_set):
    """
    Cerca oggetti geometrici di cardinalità 41.
    """
    print(f"\n{'='*70}")
    print(f"RICERCA CARDINALITÀ 41")
    print(f"{'='*70}")
    
    # Calcola cardinalità per diverse norme
    norm_counts = defaultdict(int)
    for q in hurwitz_set:
        norm = quaternion_norm(q[:4])
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        norm_counts[norm] += 1
    
    # Cerca norme con cardinalità 41
    norms_41 = [norm for norm, count in norm_counts.items() if count == 41]
    
    print(f"\nNorme con cardinalità 41: {norms_41}")
    
    if norms_41:
        for norm in norms_41:
            print(f"\nAnalisi norma {norm}:")
            elements = [q for q in hurwitz_set if quaternion_norm(q[:4]) == norm]
            print(f"  Cardinalità: {len(elements)}")
            
            # Analisi tipi
            type_counts = defaultdict(int)
            for q in elements:
                type_counts[q[4]] += 1
            print(f"  Tipi: {dict(type_counts)}")
    else:
        print("Nessuna norma con cardinalità esattamente 41 trovata.")
    
    # Cerca somme di cardinalità che diano 41
    print(f"\nCombinazioni di cardinalità che danno 41:")
    unique_counts = sorted(set(norm_counts.values()))
    
    found_combinations = []
    for i in range(len(unique_counts)):
        for j in range(i, len(unique_counts)):
            if unique_counts[i] + unique_counts[j] == 41:
                found_combinations.append((unique_counts[i], unique_counts[j]))
    
    if found_combinations:
        for c1, c2 in found_combinations:
            print(f"  {c1} + {c2} = 41")
    else:
        print("Nessuna combinazione semplice trovata.")
    
    # Cerca orbite di dimensione 41
    print(f"\nRicerca orbite di dimensione 41:")
    units = hurwitz_units()
    
    # Campiona alcuni quaternioni
    for norm in [1, 2, 3, 4, 5]:
        elements = [q[:4] for q in hurwitz_set if quaternion_norm(q[:4]) == norm]
        if elements:
            for q in elements[:3]:
                orbit = left_orbit(q, units)
                orbit_set = set((round(a, 6), round(b, 6), round(c, 6), round(d, 6))
                               for a, b, c, d in orbit)
                if len(orbit_set) == 41:
                    print(f"  ✓ Trovata orbita di dimensione 41 per norma {norm}")
    
    return norms_41


def main():
    print("=" * 70)
    print("AUDIT COMPLETO: NORME CON CARDINALITÀ 64")
    print("=" * 70)
    
    # Genera reticolo di Hurwitz
    print("\nGenerazione reticolo di Hurwitz (limit=10)...")
    hurwitz_set = hurwitz_integers(limit=10)
    print(f"Totale quaternioni: {len(hurwitz_set)}")
    
    units = hurwitz_units()
    print(f"Unità di Hurwitz (2T): {len(units)}")
    
    # Norme target con cardinalità 64
    target_norms = [76, 79, 84, 91]
    
    # Estrai quaternioni per queste norme
    print(f"\nEstrazione quaternioni per norme {target_norms}...")
    by_norm = get_hurwitz_by_norm(hurwitz_set, target_norms)
    
    # Analizza ogni norma
    analyses = {}
    for norm in target_norms:
        quaternions = by_norm.get(norm, [])
        analysis = analyze_orbit_structure(quaternions, units, norm)
        analyses[norm] = analysis
        compare_with_kristall(analysis, norm)
    
    # Ricerca cardinalità 41
    norms_41 = search_for_41(hurwitz_set)
    
    # Riepilogo
    print(f"\n{'='*70}")
    print(f"RIEPILOGO")
    print(f"{'='*70}")
    
    print(f"\nCardinalità trovate:")
    for norm in target_norms:
        count = len(by_norm.get(norm, []))
        print(f"  Norma {norm}: {count} elementi")
    
    print(f"\nNorme con cardinalità 41: {norms_41}")
    
    print(f"\nConclusione:")
    all_64 = all(len(by_norm.get(norm, [])) == 64 for norm in target_norms)
    if all_64:
        print("  ✓ Tutte e 4 le norme hanno cardinalità 64")
        print("  Prossimo passo: analizzare la struttura geometrica di ciascuna")
        print("  per vedere se corrisponde a Krist-All o ai 64 Gate")
    else:
        print("  ✗ Non tutte le norme hanno cardinalità 64")
        print("  Bisogna estendere il range di ricerca")


if __name__ == "__main__":
    main()
