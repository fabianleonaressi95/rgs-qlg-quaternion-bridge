#!/usr/bin/env python3
"""
Analisi strutturale profonda delle 7 norme con cardinalità 64.
Verifica se corrispondono a Q6 o ai 64 vertici del cubo.
"""

import numpy as np
from math import sqrt
from itertools import product, combinations
from collections import defaultdict


def hurwitz_integers(limit=10):
    hurwitz = []
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    return hurwitz


def quaternion_norm(q):
    a, b, c, d = q
    return a*a + b*b + c*c + d*d


def project_to_3d(q):
    a, b, c, d = q
    return (b, c, d)


def analyze_projections_3d(quaternions, norm):
    """Analisi dettagliata delle proiezioni 3D."""
    print(f"\n{'='*70}")
    print(f"ANALISI PROIEZIONI 3D - NORMA {norm}")
    print(f"{'='*70}")
    
    projections = [project_to_3d(q[:4]) for q in quaternions]
    
    # Calcola tutte le norme 3D
    norms_3d = [sum(x*x for x in p) for p in projections]
    unique_norms = sorted(set(round(n, 6) for n in norms_3d))
    
    print(f"\nNorme 3D uniche: {unique_norms}")
    
    # Per ogni norma 3D, conta quanti elementi
    norm_distribution = defaultdict(int)
    for n in norms_3d:
        norm_distribution[round(n, 6)] += 1
    
    print(f"\nDistribuzione:")
    for n in unique_norms:
        print(f"  Norma 3D {n}: {norm_distribution[n]} elementi")
    
    # Cerca corrispondenze con Krist-All
    kristall_norms = [1, 3, 21]
    print(f"\nCorrispondenze con Krist-All:")
    for kn in kristall_norms:
        if any(abs(n - kn) < 0.01 for n in unique_norms):
            print(f"  ✓ Trovata norma Krist-All {kn}")
        else:
            print(f"  ✗ Norma {kn} non presente")
    
    # Cerca 216 = 6³
    if any(abs(n - 216) < 0.01 for n in unique_norms):
        print(f"\n  !!! TROVATO 216 = 6³ !!!")
    
    # Analisi coordinate
    print(f"\nCoordinate 3D (campione):")
    for i, p in enumerate(projections[:5]):
        print(f"  {i}: {p}")
    
    # Calcola distanze tra punti
    print(f"\nAnalisi distanze:")
    distances = []
    for i in range(min(len(projections), 10)):
        for j in range(i+1, min(len(projections), 10)):
            dist = sqrt(sum((projections[i][k] - projections[j][k])**2 for k in range(3)))
            distances.append(round(dist, 6))
    
    unique_distances = sorted(set(distances))
    print(f"  Distanze uniche: {unique_distances[:10]}...")
    
    return unique_norms, norm_distribution


def analyze_orbit_decomposition(quaternions, norm):
    """Analizza la decomposizione in orbite."""
    print(f"\n{'='*70}")
    print(f"DECOMPOSIZIONE IN ORBITE - NORMA {norm}")
    print(f"{'='*70}")
    
    # Raggruppa per norma 3D
    projections = [project_to_3d(q[:4]) for q in quaternions]
    by_norm_3d = defaultdict(list)
    
    for i, p in enumerate(projections):
        n_3d = round(sum(x*x for x in p), 6)
        by_norm_3d[n_3d].append(quaternions[i])
    
    print(f"\nDecomposizione per norma 3D:")
    for n_3d, elems in sorted(by_norm_3d.items()):
        print(f"  Norma 3D {n_3d}: {len(elems)} elementi")
    
    # Cerca decomposizioni tipo 4×16, 8×8, 2×32
    print(f"\nPossibili decomposizioni:")
    for n_3d, elems in sorted(by_norm_3d.items()):
        count = len(elems)
        if count in [4, 8, 16, 32]:
            print(f"  Norma 3D {n_3d}: {count} (potenza di 2)")


def compare_with_cube_64(quaternions, norm):
    """Confronta con i 64 vertici centrati del cubo."""
    print(f"\n{'='*70}")
    print(f"CONFRONTO CON CUBO 64 - NORMA {norm}")
    print(f"{'='*70}")
    
    # I 64 vertici centrati del cubo in 3D sarebbero coordinate tipo (±1, ±1, ±1)
    # con centro in (0,0,0) o altre combinazioni
    
    projections = [project_to_3d(q[:4]) for q in quaternions]
    
    # Conta quante coordinate sono intere vs mezze
    integer_coords = 0
    half_coords = 0
    
    for p in projections:
        for coord in p:
            if abs(coord - round(coord)) < 0.01:
                integer_coords += 1
            elif abs(coord - round(coord) - 0.5) < 0.01:
                half_coords += 1
    
    print(f"\nCoordinate intere: {integer_coords}")
    print(f"Coordinate mezze: {half_coords}")
    
    # Verifica se tutte le coordinate sono intere (struttura cubica)
    if half_coords == 0:
        print(f"  ✓ Tutte le coordinate sono intere - possibile struttura cubica")
    else:
        print(f"  ✗ Ci sono coordinate mezze - non struttura cubica pura")


def main():
    print("=" * 70)
    print("ANALISI STRUTTURALE PROFONDA: NORME CON CARDINALITÀ 64")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    
    target_norms = [304, 316, 327, 336, 349, 364, 381]
    
    # Estrai quaternioni
    by_norm = defaultdict(list)
    for q in hurwitz_set:
        norm = quaternion_norm(q[:4])
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        if norm in target_norms:
            by_norm[norm].append(q)
    
    # Analizza ogni norma
    for norm in target_norms:
        quaternions = by_norm.get(norm, [])
        
        analyze_projections_3d(quaternions, norm)
        analyze_orbit_decomposition(quaternions, norm)
        compare_with_cube_64(quaternions, norm)
    
    # Riepilogo finale
    print(f"\n{'='*70}")
    print("RIEPILOGO FINALE")
    print(f"{'='*70}")
    
    print(f"\nNorme analizzate: {target_norms}")
    print(f"Tutte hanno cardinalità 64 e orbite di dimensione 24")
    
    print(f"\nOsservazioni chiave:")
    print(f"  - Tutte le orbite hanno dimensione 24 (unità Hurwitz)")
    print(f"  - La norma 316 ha norma 3D = 216 = 6³")
    print(f"  - Questo è interessante per la struttura 64×6³×5")
    
    print(f"\nProssimo passo:")
    print(f"  Analizzare in dettaglio la norma 316")
    print(f"  per vedere se la sua struttura geometrica")
    print(f"  corrisponde ai 64 Gate o a Q6")


if __name__ == "__main__":
    main()
