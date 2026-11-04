#!/usr/bin/env python3
"""
Ricalcolo corretto delle cardinalità nel reticolo di Hurwitz.
Trova quali norme hanno realmente cardinalità 64.
"""

from collections import defaultdict
from itertools import product


def hurwitz_integers(limit=10):
    """Genera interi di Hurwitz."""
    hurwitz = []
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    return hurwitz


def quaternion_norm(q):
    """Norma di un quaternione."""
    a, b, c, d = q
    return a*a + b*b + c*c + d*d


def main():
    print("=" * 70)
    print("RICALCOLO CARDINALITÀ NEL RETICOLO DI HURWITZ")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    print(f"\nTotale quaternioni: {len(hurwitz_set)}")
    
    # Calcola cardinalità per tutte le norme
    norm_counts = defaultdict(int)
    for q in hurwitz_set:
        norm = quaternion_norm(q[:4])
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        norm_counts[norm] += 1
    
    # Trova norme con cardinalità 64
    norms_64 = [norm for norm, count in norm_counts.items() if count == 64]
    
    print(f"\nNorme con cardinalità 64: {norms_64}")
    
    # Trova norme con cardinalità vicine a 64
    print(f"\nNorme con cardinalità vicine a 64:")
    for norm in sorted(norm_counts.keys()):
        count = norm_counts[norm]
        if 50 <= count <= 80:
            print(f"  Norma {norm}: {count} elementi")
    
    # Mostra distribuzione delle cardinalità
    print(f"\nDistribuzione delle cardinalità:")
    cardinality_dist = defaultdict(int)
    for count in norm_counts.values():
        cardinality_dist[count] += 1
    
    for count in sorted(cardinality_dist.keys())[:20]:
        print(f"  Cardinalità {count}: {cardinality_dist[count]} norme")
    
    # Verifica le norme Krist-All
    print(f"\nCardinalità delle norme Krist-All:")
    for norm in [1, 3, 21]:
        count = norm_counts.get(norm, 0)
        print(f"  Norma {norm}: {count} elementi")


if __name__ == "__main__":
    main()
