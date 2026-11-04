#!/usr/bin/env python3
"""
Analisi focalizzata sulla norma 316.
Questa norma ha 216 = 6³ come norma 3D, chiave per 64×6³×5.
"""

import numpy as np
from math import sqrt, pi
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


def main():
    print("=" * 70)
    print("FOCUS: NORMA 316 (CONTIENE 216 = 6³)")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    
    # Estrai quaternioni con norma 316
    quaternions_316 = [q for q in hurwitz_set if quaternion_norm(q[:4]) == 316]
    
    print(f"\nCardinalità: {len(quaternions_316)}")
    
    # Separa per norma 3D
    by_norm_3d = defaultdict(list)
    for q in quaternions_316:
        proj = project_to_3d(q[:4])
        n_3d = sum(x*x for x in proj)
        by_norm_3d[round(n_3d, 6)].append(q)
    
    print(f"\nDecomposizione per norma 3D:")
    for n_3d, elems in sorted(by_norm_3d.items()):
        print(f"  Norma 3D {n_3d}: {len(elems)} elementi")
    
    # Analisi dettagliata del gruppo con norma 3D = 216
    group_216 = by_norm_3d[216.0]
    print(f"\n{'='*70}")
    print(f"ANALISI GRUPPO CON NORMA 3D = 216 = 6³")
    print(f"{'='*70}")
    
    print(f"\nCardinalità: {len(group_216)}")
    
    # Mostra coordinate
    print(f"\nCoordinate 3D (prime 10):")
    for i, q in enumerate(group_216[:10]):
        proj = project_to_3d(q[:4])
        print(f"  {i}: {proj}")
    
    # Analisi coordinate
    coords_216 = [project_to_3d(q[:4]) for q in group_216]
    
    # Trova valori unici per ogni coordinata
    for coord_idx in range(3):
        values = sorted(set(round(c[coord_idx], 2) for c in coords_216))
        print(f"\nCoordinata {coord_idx}: {values}")
    
    # Analisi struttura
    print(f"\n{'='*70}")
    print(f"ANALISI STRUTTURA")
    print(f"{'='*70}")
    
    # Verifica se forma un cubo
    print(f"\n48 elementi con norma 3D = 216")
    print(f"Decomposizioni possibili:")
    print(f"  48 = 16 × 3")
    print(f"  48 = 8 × 6")
    print(f"  48 = 4 × 12")
    print(f"  48 = 2 × 24")
    
    # Calcola il raggio della sfera 3D
    radius_216 = sqrt(216)
    print(f"\nRaggio sfera 3D: √216 = {radius_216:.6f}")
    print(f"Questo è 6√6 = {6*sqrt(6):.6f}")
    
    # Analisi del gruppo con norma 3D = 300
    group_300 = by_norm_3d[300.0]
    print(f"\n{'='*70}")
    print(f"ANALISI GRUPPO CON NORMA 3D = 300")
    print(f"{'='*70}")
    
    print(f"\nCardinalità: {len(group_300)}")
    print(f"Coordinate 3D (prime 10):")
    for i, q in enumerate(group_300[:10]):
        proj = project_to_3d(q[:4])
        print(f"  {i}: {proj}")
    
    radius_300 = sqrt(300)
    print(f"\nRaggio sfera 3D: √300 = {radius_300:.6f}")
    print(f"Questo è 10√3 = {10*sqrt(3):.6f}")
    
    # Connessione con 64×6³×5
    print(f"\n{'='*70}")
    print(f"CONNESSIONE CON 64×6³×5 = 69120")
    print(f"{'='*70}")
    
    print(f"\nStruttura trovata:")
    print(f"  64 quaternioni totali")
    print(f"  48 con norma 3D = 216 = 6³")
    print(f"  16 con norma 3D = 300")
    
    print(f"\nInterpretazione:")
    print(f"  64 = 48 + 16")
    print(f"  48 = 8 × 6 (potrebbe collegarsi a 6³)")
    print(f"  16 = 2⁴")
    
    print(f"\nSe moltiplichiamo:")
    print(f"  64 × 6³ = 64 × 216 = 13824")
    print(f"  64 × 6³ × 5 = 13824 × 5 = 69120")
    
    print(f"\nLa struttura 64×6³×5 potrebbe emergere da:")
    print(f"  - 64: cardinalità della norma 316")
    print(f"  - 6³: norma 3D del sottogruppo di 48 elementi")
    print(f"  - 5: da investigare (forse fattore geometrico)")
    
    # Cerca il fattore 5
    print(f"\n{'='*70}")
    print(f"RICERCA FATTORE 5")
    print(f"{'='*70}")
    
    print(f"\nPossibili origini del fattore 5:")
    print(f"  - 5 dimensioni?")
    print(f"  - 5 gruppi di simmetria?")
    print(f"  - 5 = (216 + 300) / √(...) ")
    print(f"  - 5 emerge da rapporto tra raggi?")
    
    ratio = radius_300 / radius_216
    print(f"\nRapporto raggi: √300 / √216 = {ratio:.6f}")
    print(f"≈ 1.1785 ≈ 3π/8 = {3*pi/8:.6f}")
    
    print(f"\nConclusione:")
    print(f"  La norma 316 è un candidato forte per")
    print(f"  il 64 nella struttura 64×6³×5.")
    print(f"  Il sottogruppo di 48 elementi ha norma 3D = 6³.")
    print(f"  Il fattore 5 rimane da investigare.")


if __name__ == "__main__":
    main()
