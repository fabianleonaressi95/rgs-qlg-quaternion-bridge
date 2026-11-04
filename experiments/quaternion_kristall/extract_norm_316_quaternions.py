#!/usr/bin/env python3
"""
Estrazione delle coordinate precise dei 64 quaternioni della norma 316.
Questi sono i quaternioni necessari per il test Hopf corretto.
"""

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


def normalize_quaternion(q):
    """Normalizza un quaternione a norma 1."""
    a, b, c, d = q
    norm = sqrt(a*a + b*b + c*c + d*d)
    return (a/norm, b/norm, c/norm, d/norm)


def main():
    from math import sqrt
    
    print("=" * 70)
    print("ESTRAZIONE QUATERNIONI NORMA 316 PER TEST HOPF")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    
    # Estrai quaternioni con norma 316
    quaternions_316 = [q for q in hurwitz_set if quaternion_norm(q[:4]) == 316]
    
    print(f"\nTotale quaternioni con norma 316: {len(quaternions_316)}")
    
    # Separa per norma 3D
    from collections import defaultdict
    by_norm_3d = defaultdict(list)
    for q in quaternions_316:
        proj = (q[1], q[2], q[3])
        n_3d = proj[0]*proj[0] + proj[1]*proj[1] + proj[2]*proj[2]
        by_norm_3d[round(n_3d, 6)].append(q)
    
    print(f"\nDecomposizione per norma 3D:")
    for n_3d, elems in sorted(by_norm_3d.items()):
        print(f"  Norma 3D {n_3d}: {len(elems)} elementi")
    
    # Esporta coordinate in formato utilizzabile
    print(f"\n{'='*70}")
    print(f"COORDINATE QUATERNIONI (FORMATO: t, x, y, z)")
    print(f"{'='*70}")
    
    # Gruppo 216 (48 elementi)
    print(f"\n# GRUPPO 216 (48 elementi, norma 3D = 216 = 6³)")
    group_216 = by_norm_3d[216.0]
    for i, q in enumerate(group_216):
        t, x, y, z, qtype = q
        print(f"q_{i}_216 = ({t}, {x}, {y}, {z})  # {qtype}")
    
    # Gruppo 300 (16 elementi)
    print(f"\n# GRUPPO 300 (16 elementi, norma 3D = 300, vertici cubo)")
    group_300 = by_norm_3d[300.0]
    for i, q in enumerate(group_300):
        t, x, y, z, qtype = q
        print(f"q_{i}_300 = ({t}, {x}, {y}, {z})  # {qtype}")
    
    # Normalizza per il test Hopf
    print(f"\n{'='*70}")
    print(f"QUATERNIONI NORMALIZZATI (norma = 1)")
    print(f"{'='*70}")
    
    print(f"\n# GRUPPO 216 NORMALIZZATO")
    for i, q in enumerate(group_216):
        t, x, y, z, qtype = q
        norm = sqrt(t*t + x*x + y*y + z*z)
        tn, xn, yn, zn = t/norm, x/norm, y/norm, z/norm
        print(f"q_{i}_216_norm = ({tn:.10f}, {xn:.10f}, {yn:.10f}, {zn:.10f})")
    
    print(f"\n# GRUPPO 300 NORMALIZZATO")
    for i, q in enumerate(group_300):
        t, x, y, z, qtype = q
        norm = sqrt(t*t + x*x + y*y + z*z)
        tn, xn, yn, zn = t/norm, x/norm, y/norm, z/norm
        print(f"q_{i}_300_norm = ({tn:.10f}, {xn:.10f}, {yn:.10f}, {zn:.10f})")
    
    # Esporta in formato CSV per facile import
    print(f"\n{'='*70}")
    print(f"FORMATO CSV PER IMPORT")
    print(f"{'='*70}")
    
    print(f"\nt,x,y,z,group")
    for q in group_216:
        t, x, y, z, _ = q
        print(f"{t},{x},{y},{z},216")
    for q in group_300:
        t, x, y, z, _ = q
        print(f"{t},{x},{y},{z},300")
    
    # Verifica distintezza
    print(f"\n{'='*70}")
    print(f"VERIFICA DISTINTEZZA")
    print(f"{'='*70}")
    
    all_quats = [(q[0], q[1], q[2], q[3]) for q in quaternions_316]
    unique_quats = set(all_quats)
    
    print(f"\nTotale quaternioni: {len(all_quats)}")
    print(f"Quaternioni distinti: {len(unique_quats)}")
    
    if len(all_quats) == len(unique_quats):
        print("✓ Tutti i quaternioni sono distinti")
    else:
        print("✗ Ci sono quaternioni duplicati")
    
    # Calcola distanze minime
    print(f"\nAnalisi distanze:")
    min_dist = float('inf')
    for i in range(len(all_quats)):
        for j in range(i+1, len(all_quats)):
            dist = sqrt(sum((all_quats[i][k] - all_quats[j][k])**2 for k in range(4)))
            if dist < min_dist:
                min_dist = dist
    
    print(f"  Distanza minima tra quaternioni: {min_dist:.6f}")
    
    if min_dist > 0.01:
        print("  ✓ I quaternioni sono ben separati")
    else:
        print("  ✗ Alcuni quaternioni sono molto vicini")


if __name__ == "__main__":
    main()
