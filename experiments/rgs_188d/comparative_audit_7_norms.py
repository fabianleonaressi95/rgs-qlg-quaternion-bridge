#!/usr/bin/env python3
"""
Audit comparativo delle 7 norme con cardinalità 64:
304, 316, 327, 336, 349, 364, 381

Pipeline per ogni norma:
- Cardinalità
- Decomposizione per norma 3D
- Gram matrix (rango, autovalori)
- Invarianza O_h completa
- Fibre Hopf
- Firma geometrica completa
"""

import numpy as np
from math import sqrt, acos, pi
from itertools import product, permutations
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
        a1*d2 + b1*c2 - c1*b1 + d1*a1
    )


def normalize_quaternion(q):
    norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
    return (q[0]/norm, q[1]/norm, q[2]/norm, q[3]/norm)


def hopf_map(q):
    i_quat = (0, 1, 0, 0)
    q_norm = normalize_quaternion(q)
    conj_q = quaternion_conjugate(q_norm)
    temp = quaternion_multiply(q_norm, i_quat)
    result = quaternion_multiply(temp, conj_q)
    return result[1:]


def test_oh_invariance(quaternions):
    """Testa invarianza sotto O_h completo (permutazioni + cambi di segno)."""
    coords = [q[:4] for q in quaternions]
    coord_set = set((round(c, 6) for c in q) for q in coords)
    
    # Testa permutazioni di x,y,z
    perm_invariant = True
    for perm in permutations([1, 2, 3]):
        permuted_set = set()
        for q in coords:
            permuted = (q[0], q[perm[0]], q[perm[1]], q[perm[2]])
            permuted_set.add(tuple(round(c, 6) for c in permuted))
        if permuted_set != coord_set:
            perm_invariant = False
            break
    
    # Testa cambi di segno di tutte le coordinate (incluso t)
    sign_invariant = True
    signs_list = list(product([1, -1], repeat=4))
    for signs in signs_list:
        signed_set = set()
        for q in coords:
            signed = (signs[0]*q[0], signs[1]*q[1], signs[2]*q[2], signs[3]*q[3])
            signed_set.add(tuple(round(c, 6) for c in signed))
        if signed_set != coord_set:
            sign_invariant = False
            break
    
    return perm_invariant, sign_invariant


def analyze_norm(hurwitz_set, target_norm):
    """Analisi completa di una norma."""
    print(f"\n{'='*70}")
    print(f"NORMA {target_norm}")
    print(f"{'='*70}")
    
    # Estrai quaternioni
    quaternions = [q for q in hurwitz_set if quaternion_norm(q[:4]) == target_norm]
    
    if len(quaternions) == 0:
        print(f"NESSUN QUATERNIONE TROVATO")
        return None
    
    print(f"\nCardinalità: {len(quaternions)}")
    
    # Decomposizione per norma 3D
    by_norm_3d = defaultdict(list)
    for q in quaternions:
        proj = (q[1], q[2], q[3])
        n_3d = proj[0]*proj[0] + proj[1]*proj[1] + proj[2]*proj[2]
        by_norm_3d[round(n_3d, 6)].append(q)
    
    print(f"\nDecomposizione per norma 3D:")
    norm_3d_dist = sorted([(n, len(elems)) for n, elems in by_norm_3d.items()])
    for n_3d, count in norm_3d_dist:
        print(f"  r² = {n_3d}: {count} elementi")
    
    # Applica Hopf
    hopf_images = [hopf_map(q[:4]) for q in quaternions]
    
    # Verifica distintezza
    unique_images = set((round(p[0], 6), round(p[1], 6), round(p[2], 6)) for p in hopf_images)
    print(f"\nHopf: {len(quaternions)} → {len(unique_images)} immagini")
    
    # Gram matrix
    n = len(hopf_images)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(hopf_images[i][k] * hopf_images[j][k] for k in range(3))
    
    eigenvalues = np.linalg.eigvalsh(G)
    eigenvalues = np.round(eigenvalues, 6)
    rank = np.linalg.matrix_rank(G)
    
    nonzero_eigenvals = [ev for ev in eigenvalues if abs(ev) > 1e-6]
    
    print(f"\nGram matrix:")
    print(f"  Rango: {rank}")
    print(f"  Autovalori non nulli: {nonzero_eigenvals}")
    
    # Invarianza O_h
    perm_inv, sign_inv = test_oh_invariance(quaternions)
    print(f"\nInvarianza O_h:")
    print(f"  Permutazioni xyz: {'✓' if perm_inv else '✗'}")
    print(f"  Cambi di segno: {'✓' if sign_inv else '✗'}")
    
    # Firma geometrica
    signature = {
        'norm': target_norm,
        'cardinality': len(quaternions),
        'norm_3d_distribution': norm_3d_dist,
        'hopf_collapse': len(quaternions) - len(unique_images),
        'gram_rank': rank,
        'gram_eigenvalues': nonzero_eigenvals,
        'oh_perm_invariant': perm_inv,
        'oh_sign_invariant': sign_inv
    }
    
    return signature


def compare_signatures(signatures):
    """Confronta le firme delle diverse norme."""
    print(f"\n{'='*70}")
    print("CONFRONTO FIRME GEOMETRICHE")
    print(f"{'='*70}")
    
    print(f"\n{'Norma':<10} {'Card':<6} {'Rango':<6} {'Autovalori':<30} {'O_h perm':<10} {'O_h sign':<10}")
    print("-" * 90)
    
    for sig in signatures:
        if sig is None:
            continue
        
        ev_str = str(sig['gram_eigenvalues'])[:28]
        perm_str = '✓' if sig['oh_perm_invariant'] else '✗'
        sign_str = '✓' if sig['oh_sign_invariant'] else '✗'
        
        print(f"{sig['norm']:<10} {sig['cardinality']:<6} {sig['gram_rank']:<6} {ev_str:<30} {perm_str:<10} {sign_str:<10}")
    
    # Analisi decomposizioni
    print(f"\n{'='*70}")
    print("DECOMPOSIZIONI NORMA 3D")
    print(f"{'='*70}")
    
    for sig in signatures:
        if sig is None:
            continue
        
        print(f"\nNorma {sig['norm']}:")
        for n_3d, count in sig['norm_3d_distribution']:
            print(f"  r² = {n_3d}: {count} elementi")
    
    # Cerca pattern comuni
    print(f"\n{'='*70}")
    print("PATTERN COMUNI")
    print(f"{'='*70}")
    
    # Tutte hanno rango 3?
    all_rank_3 = all(sig['gram_rank'] == 3 for sig in signatures if sig is not None)
    print(f"\nTutte hanno rango Gram = 3: {'✓' if all_rank_3 else '✗'}")
    
    # Tutte hanno invarianza permutazioni?
    all_perm_inv = all(sig['oh_perm_invariant'] for sig in signatures if sig is not None)
    print(f"Tutte hanno invarianza permutazioni xyz: {'✓' if all_perm_inv else '✗'}")
    
    # Tutte hanno invarianza segno?
    all_sign_inv = all(sig['oh_sign_invariant'] for sig in signatures if sig is not None)
    print(f"Tutte hanno invarianza cambi di segno: {'✓' if all_sign_inv else '✗'}")
    
    # Cerca norme 3D comuni
    all_norm_3d = set()
    for sig in signatures:
        if sig is not None:
            for n_3d, _ in sig['norm_3d_distribution']:
                all_norm_3d.add(n_3d)
    
    print(f"\nNorme 3D trovate complessivamente: {sorted(all_norm_3d)}")
    
    # Cerca 216 = 6³
    has_216 = []
    for sig in signatures:
        if sig is not None:
            for n_3d, count in sig['norm_3d_distribution']:
                if abs(n_3d - 216) < 0.01:
                    has_216.append((sig['norm'], count))
    
    if has_216:
        print(f"\nNorme con r² = 216 = 6³:")
        for norm, count in has_216:
            print(f"  Norma {norm}: {count} elementi")
    else:
        print(f"\nNessuna norma ha r² = 216")


def main():
    print("=" * 70)
    print("AUDIT COMPARATIVO: 7 NORME CON CARDINALITÀ 64")
    print("=" * 70)
    
    hurwitz_set = hurwitz_integers(limit=10)
    target_norms = [304, 316, 327, 336, 349, 364, 381]
    
    signatures = []
    
    for norm in target_norms:
        sig = analyze_norm(hurwitz_set, norm)
        signatures.append(sig)
    
    # Confronta firme
    compare_signatures(signatures)
    
    # Riepilogo
    print(f"\n{'='*70}")
    print("RIEPILOGO")
    print(f"{'='*70}")
    
    print(f"\nNorme analizzate: {target_norms}")
    print(f"Cardinalità: tutte 64")
    
    # Identifica 316 come speciale
    sig_316 = next((sig for sig in signatures if sig and sig['norm'] == 316), None)
    
    if sig_316:
        print(f"\nNorma 316:")
        print(f"  ✓ Ha r² = 216 = 6³")
        print(f"  ✓ Rango Gram = 3")
        print(f"  ✓ Invarianza O_h")
        
        # Confronta con altre
        other_has_216 = any(sig for sig in signatures if sig and sig['norm'] != 316 
                           and any(abs(n - 216) < 0.01 for n, _ in sig['norm_3d_distribution']))
        
        if not other_has_216:
            print(f"\n  ✗ 316 è L'UNICA con r² = 6³")
        else:
            print(f"\n  → Altre norme hanno anche r² = 6³")


if __name__ == "__main__":
    main()
