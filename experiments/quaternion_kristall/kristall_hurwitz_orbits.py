#!/usr/bin/env python3
"""
Analisi approfondita delle orbite nel reticolo di Hurwitz e relazione con Krist-All.
Verifica della catena: Krist-All → ImH → H → Hurwitz lattice → S³
"""

import numpy as np
from math import isqrt, sqrt
from itertools import product, permutations
from collections import defaultdict


def hurwitz_integers(limit=5):
    """Genera gli interi di Hurwitz in un range limitato."""
    hurwitz = []
    
    # Tipo 1: tutti coefficienti interi
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    
    # Tipo 2: mezzi-interi
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    
    return hurwitz


def hurwitz_units():
    """Restituisce le 24 unità del reticolo di Hurwitz."""
    units = []
    
    # 8 unità standard
    for sign in [1, -1]:
        units.append((sign, 0, 0, 0))
        units.append((0, sign, 0, 0))
        units.append((0, 0, sign, 0))
        units.append((0, 0, 0, sign))
    
    # 16 unità di tipo (±1 ± i ± j ± k)/2
    for signs in product([1, -1], repeat=4):
        units.append((signs[0]/2, signs[1]/2, signs[2]/2, signs[3]/2))
    
    return units


def quaternion_multiply(q1, q2):
    """Moltiplicazione di quaternioni."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2
    )


def quaternion_conjugate(q):
    """Coniugato di un quaternione."""
    a, b, c, d = q
    return (a, -b, -c, -d)


def quaternion_norm(q):
    """Norma di un quaternione."""
    a, b, c, d = q
    return a*a + b*b + c*c + d*d


def quaternion_inverse(q):
    """Inverso di un quaternione non-zero."""
    norm = quaternion_norm(q)
    if norm == 0:
        return None
    conj = quaternion_conjugate(q)
    return (conj[0]/norm, conj[1]/norm, conj[2]/norm, conj[3]/norm)


def left_multiply(q, hurwitz_set):
    """Orbita sinistra: {q * h | h in hurwitz_set}."""
    return [quaternion_multiply(q, h) for h in hurwitz_set]


def right_multiply(hurwitz_set, q):
    """Orbita destra: {h * q | h in hurwitz_set}."""
    return [quaternion_multiply(h, q) for h in hurwitz_set]


def kristall_points():
    """
    Genera i 62 punti Krist-All.
    Basato sulla struttura: 64×6×6×6×5 ma con 2 punti mancanti.
    """
    # I tre gusci fondamentali
    shell1 = [(1, 0, 0)]
    shell2 = [(1, 1, 1)]
    shell3 = [(1, 2, 4)]
    
    # Genera permutazioni e varianti con segni
    points = set()
    
    # Guscio 1: permutazioni di (1,0,0) con segni
    for perm in set(permutations([1, 0, 0])):
        for signs in product([1, -1], repeat=3):
            if signs[0] * signs[1] * signs[2] == 1:  # Det = +1
                point = tuple(signs[i] * perm[i] for i in range(3))
                points.add(point)
    
    # Guscio 2: permutazioni di (1,1,1) con segni
    for signs in product([1, -1], repeat=3):
        if signs[0] * signs[1] * signs[2] == 1:
            points.add(tuple(signs[i] * 1 for i in range(3)))
    
    # Guscio 3: permutazioni di (1,2,4) con segni
    for perm in set(permutations([1, 2, 4])):
        for signs in product([1, -1], repeat=3):
            if signs[0] * signs[1] * signs[2] == 1:
                point = tuple(signs[i] * perm[i] for i in range(3))
                points.add(point)
    
    return list(points)


def kristall_to_quaternions(points):
    """Converte punti Krist-All (in ImH) a quaternioni puri."""
    return [(0, x, y, z) for (x, y, z) in points]


def find_orbit_representatives(hurwitz_set, target_norms, max_coord=5):
    """
    Trova rappresentanti di orbite con norme specifiche nel reticolo di Hurwitz.
    """
    orbit_reps = defaultdict(list)
    
    for q in hurwitz_set:
        t, x, y, z, qtype = q
        norm = t*t + x*x + y*y + z*z
        # Converti a intero se è .0
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)

        if norm in target_norms and max(abs(t), abs(x), abs(y), abs(z)) <= max_coord:
            orbit_reps[norm].append((t, x, y, z, qtype))
    
    return orbit_reps


def analyze_orbit_structure(representatives, units):
    """
    Analizza la struttura delle orbite: calcola dimensioni delle orbite
    tramite azione delle unità di Hurwitz.
    """
    orbit_analysis = {}
    
    for norm, reps in representatives.items():
        orbit_analysis[norm] = {
            'total_count': len(reps),
            'orbit_sizes': defaultdict(int),
            'sample_orbits': []
        }
        
        # Per ogni rappresentante, calcola l'orbita sinistra e destra
        for rep in reps[:10]:  # Campiona i primi 10
            q = rep[:4]
            
            # Orbita sinistra
            left_orbit = left_multiply(q, units)
            left_orbit_set = set((round(a, 6), round(b, 6), round(c, 6), round(d, 6))
                                 for a, b, c, d in left_orbit)
            left_size = len(left_orbit_set)
            
            # Orbita destra
            right_orbit = right_multiply(units, q)
            right_orbit_set = set((round(a, 6), round(b, 6), round(c, 6), round(d, 6))
                                  for a, b, c, d in right_orbit)
            right_size = len(right_orbit_set)
            
            orbit_analysis[norm]['orbit_sizes'][left_size] += 1
            orbit_analysis[norm]['sample_orbits'].append({
                'representative': q,
                'left_orbit_size': left_size,
                'right_orbit_size': right_size
            })
    
    return orbit_analysis


def check_kristall_in_hurwitz(kristall_points, hurwitz_set):
    """
    Verifica quali punti Krist-All appartengono al reticolo di Hurwitz.
    """
    kristall_in_hurwitz = []
    kristall_not_in_hurwitz = []
    
    hurwitz_set_tuples = set((round(t, 6), round(x, 6), round(y, 6), round(z, 6))
                           for t, x, y, z, _ in hurwitz_set)
    
    for point in kristall_points:
        # Converti a quaternione puro (t=0)
        q_pure = (0, point[0], point[1], point[2])
        
        if q_pure in hurwitz_set_tuples:
            kristall_in_hurwitz.append(point)
        else:
            kristall_not_in_hurwitz.append(point)
    
    return kristall_in_hurwitz, kristall_not_in_hurwitz


def analyze_cardinalities():
    """
    Analizza le cardinalità delle orbite nel reticolo di Hurwitz.
    Cerca orbite con cardinalità 4, 8, 16, 24, 48, 64.
    """
    hurwitz = hurwitz_integers(limit=5)
    units = hurwitz_units()
    
    norm_counts = defaultdict(int)
    orbit_cardinalities = defaultdict(set)
    
    for q in hurwitz:
        t, x, y, z, qtype = q
        norm = t*t + x*x + y*y + z*z
        # Converti a intero se è .0 (da mezzi-interi)
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        norm_counts[norm] += 1
    
    # Calcola cardinalità attese basate sulla struttura
    target_cardinalities = [4, 8, 16, 24, 48, 64]
    
    cardinality_analysis = {}
    for target in target_cardinalities:
        cardinality_analysis[target] = []
        for norm, count in norm_counts.items():
            if count == target:
                cardinality_analysis[target].append(norm)
    
    return norm_counts, cardinality_analysis


def main():
    print("=" * 70)
    print("ANALISI ORBITE KRIST-ALL NEL RETICOLO DI HURWITZ")
    print("=" * 70)
    
    # 1. Genera strutture
    print("\n1. GENERAZIONE STRUTTURE")
    print("-" * 70)
    
    kristall_pts = kristall_points()
    print(f"Punti Krist-All generati: {len(kristall_pts)}")
    
    kristall_quats = kristall_to_quaternions(kristall_pts)
    print(f"Quaternioni puri Krist-All: {len(kristall_quats)}")
    
    hurwitz = hurwitz_integers(limit=5)
    print(f"Interi di Hurwitz (limit=5): {len(hurwitz)}")
    
    units = hurwitz_units()
    print(f"Unità di Hurwitz: {len(units)}")
    
    # 2. Verifica appartenenza Krist-All nel reticolo di Hurwitz
    print("\n\n2. APPARTENENZA KRIST-ALL AL RETICOLO DI HURWITZ")
    print("-" * 70)
    
    in_hurwitz, not_in_hurwitz = check_kristall_in_hurwitz(kristall_pts, hurwitz)
    
    print(f"Punti Krist-All nel reticolo di Hurwitz: {len(in_hurwitz)}")
    print(f"Punti Krist-All NON nel reticolo: {len(not_in_hurwitz)}")
    
    if in_hurwitz:
        print("\nEsempi di punti Krist-All nel reticolo:")
        for p in in_hurwitz[:10]:
            print(f"  {p}")
    
    if not_in_hurwitz:
        print("\nEsempi di punti Krist-All NON nel reticolo:")
        for p in not_in_hurwitz[:10]:
            print(f"  {p}")
    
    # 3. Analisi norme di interesse
    print("\n\n3. ANALISI NORME DI INTERESSE (1, 3, 21)")
    print("-" * 70)
    
    target_norms = [1, 3, 21]
    orbit_reps = find_orbit_representatives(hurwitz, target_norms, max_coord=5)
    
    for norm in target_norms:
        print(f"\nNorma {norm}:")
        print(f"  Rappresentanti trovati: {len(orbit_reps[norm])}")
        if orbit_reps[norm]:
            print(f"  Esempi:")
            for rep in orbit_reps[norm][:5]:
                t, x, y, z, qtype = rep
                print(f"    ({t}, {x}, {y}, {z}) [{qtype}]")
    
    # 4. Analisi struttura orbite
    print("\n\n4. STRUTTURA DELLE ORBITE")
    print("-" * 70)
    
    orbit_analysis = analyze_orbit_structure(orbit_reps, units)
    
    for norm in target_norms:
        print(f"\nNorma {norm}:")
        print(f"  Cardinalità totale: {orbit_analysis[norm]['total_count']}")
        print(f"  Dimensioni orbite sinistra:")
        for size, count in sorted(orbit_analysis[norm]['orbit_sizes'].items()):
            print(f"    {size}: {count} rappresentanti")
        
        if orbit_analysis[norm]['sample_orbits']:
            print(f"  Esempio orbita:")
            sample = orbit_analysis[norm]['sample_orbits'][0]
            print(f"    Rappresentante: {sample['representative']}")
            print(f"    Dimensione orbita sinistra: {sample['left_orbit_size']}")
            print(f"    Dimensione orbita destra: {sample['right_orbit_size']}")
    
    # 5. Analisi cardinalità
    print("\n\n5. ANALISI CARDINALITÀ")
    print("-" * 70)
    
    norm_counts, cardinality_analysis = analyze_cardinalities()
    
    print("Cardinalità target e norme corrispondenti:")
    for target in [4, 8, 16, 24, 48, 64]:
        norms = cardinality_analysis[target]
        print(f"  Cardinalità {target:2d}: {len(norms)} norme")
        if norms:
            norms_str = [str(int(n) if isinstance(n, float) and n.is_integer() else n) for n in norms]
            print(f"    Norme: {norms_str[:10]}")
    
    print("\nNorme con cardinalità > 100:")
    for norm, count in sorted(norm_counts.items()):
        if count > 100:
            print(f"  Norma {norm:3d}: {count:4d} elementi")
    
    # 6. Verifica struttura 64
    print("\n\n6. VERIFICA STRUTTURA 64")
    print("-" * 70)
    
    # Cerca norme con cardinalità 64
    norms_64 = cardinality_analysis[64]
    print(f"Norme con cardinalità 64: {len(norms_64)}")
    if norms_64:
        print(f"  {norms_64}")
    else:
        print("  Nessuna norma con cardinalità esattamente 64 trovata.")
        print("  Cardinalità più vicine a 64:")
        closest = sorted([(abs(count - 64), norm, count) for norm, count in norm_counts.items()])[:10]
        for diff, norm, count in closest:
            norm_str = int(norm) if isinstance(norm, float) and norm.is_integer() else norm
            print(f"    Norma {norm_str}: cardinalità {count} (diff: {diff})")
    
    # 7. Riepilogo
    print("\n\n" + "=" * 70)
    print("RIEPILOGO")
    print("=" * 70)
    
    print(f"\nPunti Krist-All totali: {len(kristall_pts)}")
    print(f"Punti nel reticolo di Hurwitz: {len(in_hurwitz)}")
    print(f"Punti fuori dal reticolo: {len(not_in_hurwitz)}")
    
    print(f"\nCardinalità nel reticolo di Hurwitz:")
    for norm in [1, 3, 21]:
        count = norm_counts.get(norm, 0)
        print(f"  Norma {norm}: {count} elementi")
    
    if norms_64:
        print(f"\n✓ Trovata struttura con cardinalità 64")
    else:
        print(f"\n✗ Nessuna cardinalità 64 trovata - la struttura 64×6×6×6×5")
        print(f"  potrebbe richiedere un'interpretazione diversa (es. prodotto di orbite)")


if __name__ == "__main__":
    main()
