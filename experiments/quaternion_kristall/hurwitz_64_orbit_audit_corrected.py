#!/usr/bin/env python3
"""
Audit corretto delle norme con cardinalità 64 nel reticolo di Hurwitz.
Norme corrette: 304, 316, 327, 336, 349, 364, 381
"""

import numpy as np
from math import sqrt
from itertools import product
from collections import defaultdict


def hurwitz_integers(limit=10):
    """Genera interi di Hurwitz."""
    hurwitz = []
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    return hurwitz


def hurwitz_units():
    """Restituisce le 24 unità di Hurwitz."""
    units = []
    for sign in [1, -1]:
        units.append((sign, 0, 0, 0))
        units.append((0, sign, 0, 0))
        units.append((0, 0, sign, 0))
        units.append((0, 0, 0, sign))
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


def left_orbit(q, units):
    """Orbita sinistra."""
    return [quaternion_multiply(u, q) for u in units]


def project_to_3d(q):
    """Proiezione su ImH."""
    a, b, c, d = q
    return (b, c, d)


def analyze_norm(quaternions, units, norm):
    """Analizza una norma specifica."""
    print(f"\n{'='*70}")
    print(f"NORMA {norm}")
    print(f"{'='*70}")
    
    n = len(quaternions)
    print(f"Cardinalità: {n}")
    
    if n == 0:
        return None
    
    # Analisi orbite
    print("\nOrbite sinistra (campione):")
    orbit_sizes = defaultdict(int)
    for i, q in enumerate(quaternions[:min(5, n)]):
        q_4d = q[:4]
        orbit = left_orbit(q_4d, units)
        orbit_set = set((round(a, 6), round(b, 6), round(c, 6), round(d, 6))
                       for a, b, c, d in orbit)
        orbit_size = len(orbit_set)
        orbit_sizes[orbit_size] += 1
        print(f"  Quaternione {i}: orbita {orbit_size}")
    
    # Proiezioni 3D
    projections = [project_to_3d(q[:4]) for q in quaternions]
    norms_3d = sorted(set(round(sum(x*x for x in p), 2) for p in projections))
    print(f"\nNorme 3D: {norms_3d[:10]}...")
    
    # Tipi
    type_counts = defaultdict(int)
    for q in quaternions:
        type_counts[q[4]] += 1
    print(f"Tipi: {dict(type_counts)}")
    
    return {
        'cardinality': n,
        'orbit_sizes': dict(orbit_sizes),
        'norms_3d': norms_3d,
        'type_counts': dict(type_counts)
    }


def main():
    print("=" * 70)
    print("AUDIT CORRETTO: NORME CON CARDINALITÀ 64")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    units = hurwitz_units()
    
    # Norme corrette con cardinalità 64
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
        analyze_norm(quaternions, units, norm)
    
    # Riepilogo
    print(f"\n{'='*70}")
    print("RIEPILOGO")
    print(f"{'='*70}")
    
    for norm in target_norms:
        count = len(by_norm.get(norm, []))
        print(f"Norma {norm}: {count} elementi")
    
    all_64 = all(len(by_norm.get(norm, [])) == 64 for norm in target_norms)
    if all_64:
        print("\n✓ Tutte le norme hanno cardinalità 64")
    else:
        print("\n✗ Non tutte le norme hanno cardinalità 64")


if __name__ == "__main__":
    main()
