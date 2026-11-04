#!/usr/bin/env python3
"""
Analisi matematica di Krist-All nel contesto quaternionico.
Investigazione delle equazioni diofantee e struttura del reticolo di Hurwitz.
"""

import numpy as np
from math import isqrt, sqrt
from itertools import product
from collections import defaultdict


def is_perfect_square(n):
    """Verifica se n è un quadrato perfetto."""
    if n < 0:
        return False
    root = isqrt(n)
    return root * root == n


def find_diophantine_solutions(max_R2=10000):
    """
    Cerca R² tale che R² - 1, R² - 3, R² - 21 siano tutti quadrati interi.
    
    Problema: trovare R² intero tale che:
    - R² - 1 = a²
    - R² - 3 = b²
    - R² - 21 = c²
    
    per alcuni interi a, b, c.
    """
    solutions = []
    
    for R2 in range(21, max_R2 + 1):
        diffs = [R2 - 1, R2 - 3, R2 - 21]
        
        if all(is_perfect_square(d) for d in diffs):
            a, b, c = [isqrt(d) for d in diffs]
            solutions.append({
                'R2': R2,
                'R': sqrt(R2),
                'a': a,
                'b': b,
                'c': c,
                'diffs': diffs
            })
    
    return solutions


def analyze_partial_solutions(max_R2=10000):
    """
    Analizza soluzioni parziali dove solo alcune delle differenze sono quadrati.
    """
    partial_stats = defaultdict(int)
    partial_examples = defaultdict(list)
    
    for R2 in range(21, max_R2 + 1):
        squares = []
        for idx, offset in enumerate([1, 3, 21]):
            if is_perfect_square(R2 - offset):
                squares.append(idx)
        
        if squares:
            key = tuple(squares)
            partial_stats[key] += 1
            if len(partial_examples[key]) < 5:
                partial_examples[key].append(R2)
    
    return partial_stats, partial_examples


def quaternion_norm(t, x, y, z):
    """Calcola la norma del quaternione q = t + xi + yj + zk."""
    return t*t + x*x + y*y + z*z


def hurwitz_integers(limit=5):
    """
    Genera gli interi di Hurwitz in un range limitato.
    Gli interi di Hurwitz sono quaternioni della forma:
    - a + bi + cj + dk con tutti coefficienti interi
    - (a + 1/2) + (b + 1/2)i + (c + 1/2)j + (d + 1/2)k
    """
    hurwitz = []
    
    # Tipo 1: tutti coefficienti interi
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    
    # Tipo 2: mezzi-interi (tutti componenti sono interi + 1/2)
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    
    return hurwitz


def hurwitz_units():
    """
    Restituisce le 24 unità del reticolo di Hurwitz.
    """
    units = []
    
    # 8 unità standard: ±1, ±i, ±j, ±k
    for sign in [1, -1]:
        units.append((sign, 0, 0, 0))  # ±1
        units.append((0, sign, 0, 0))  # ±i
        units.append((0, 0, sign, 0))  # ±j
        units.append((0, 0, 0, sign))  # ±k
    
    # 16 unità di tipo (±1 ± i ± j ± k)/2
    for signs in product([1, -1], repeat=4):
        units.append((signs[0]/2, signs[1]/2, signs[2]/2, signs[3]/2))
    
    return units


def kristall_shells():
    """
    Definisce i tre gusci Krist-All come vettori in ImH.
    Gusci originali: (1,0,0), (1,1,1), (1,2,4)
    """
    shells = [
        np.array([1, 0, 0]),    # norma = 1
        np.array([1, 1, 1]),    # norma = 3
        np.array([1, 2, 4])     # norma = 21
    ]
    return shells


def kristall_quaternion_completions(max_t=20):
    """
    Per ogni guscio Krist-All, trova i valori di t che portano su una sfera S³
    di raggio R.
    
    Per un vettore v = (x, y, z) con norma N = x² + y² + z²,
    cerchiamo t tale che t² + N = R².
    """
    shells = kristall_shells()
    completions = {}
    
    for idx, shell in enumerate(shells):
        N = np.sum(shell**2)
        completions[idx] = {
            'shell': shell,
            'norm': N,
            'solutions': []
        }
        
        for t in range(-max_t, max_t + 1):
            R2 = t*t + N
            if is_perfect_square(R2):
                R = isqrt(R2)
                completions[idx]['solutions'].append({
                    't': t,
                    'R': R,
                    'R2': R2
                })
    
    return completions


def analyze_hurwitz_orbits(max_norm=50):
    """
    Analizza le orbite nel reticolo di Hurwitz e verifica se compaiono
    le norme 1, 3, 21.
    """
    units = hurwitz_units()
    hurwitz = hurwitz_integers(limit=5)
    
    # Calcola le norme di tutti gli interi di Hurwitz
    norm_histogram = defaultdict(int)
    norm_examples = defaultdict(list)
    
    for q in hurwitz:
        t, x, y, z, qtype = q
        norm = quaternion_norm(t, x, y, z)
        
        if norm <= max_norm:
            norm_histogram[norm] += 1
            if len(norm_examples[norm]) < 3:
                norm_examples[norm].append((t, x, y, z, qtype))
    
    # Verifica norme di interesse
    target_norms = [1, 3, 21]
    target_info = {}
    
    for norm in target_norms:
        target_info[norm] = {
            'count': norm_histogram.get(norm, 0),
            'examples': norm_examples.get(norm, [])
        }
    
    return norm_histogram, target_info


def main():
    print("=" * 70)
    print("ANALISI KRIST-ALL E QUATERNIONI")
    print("=" * 70)
    
    # 1. Equazioni diofantee
    print("\n1. EQUAZIONI DIOFANTEE")
    print("-" * 70)
    print("Cerchiamo R² tale che R² - 1, R² - 3, R² - 21 siano tutti quadrati interi.")
    
    solutions = find_diophantine_solutions(max_R2=10000)
    
    if solutions:
        print(f"\nTrovate {len(solutions)} soluzioni:")
        for sol in solutions:
            print(f"  R² = {sol['R2']}, R = {sol['R']:.4f}")
            print(f"    R² - 1 = {sol['diffs'][0]} = {sol['a']}²")
            print(f"    R² - 3 = {sol['diffs'][1]} = {sol['b']}²")
            print(f"    R² - 21 = {sol['diffs'][2]} = {sol['c']}²")
    else:
        print("\nNessuna soluzione trovata fino a R² = 10000.")
        print("Questo suggerisce che 1, 3, 21 non possono essere completati")
        print("simultaneamente su una stessa sfera quaternionica con coordinate intere.")
    
    # Analisi parziale
    print("\n\nAnalisi soluzioni parziali:")
    partial_stats, partial_examples = analyze_partial_solutions(max_R2=10000)
    
    for key, count in sorted(partial_stats.items()):
        offsets = ["1", "3", "21"]
        active = [offsets[i] for i in key]
        print(f"  R² - {', '.join(active)} è quadrato: {count} soluzioni")
        print(f"    Esempi: {partial_examples[key][:5]}")
    
    # 2. Completamento quaternionico
    print("\n\n2. COMPLETAMENTO QUATERNIONICO")
    print("-" * 70)
    print("Per ogni guscio Krist-All, cerchiamo t tale che t² + N = R².")
    
    completions = kristall_quaternion_completions(max_t=50)
    
    for idx, comp in completions.items():
        print(f"\nGuscio {idx}: {comp['shell']}, norma N = {comp['norm']}")
        print(f"  Soluzioni (t, R):")
        for sol in comp['solutions'][:10]:  # Mostra solo le prime 10
            print(f"    t = {sol['t']:3d}, R = {sol['R']:3d}, R² = {sol['R2']:4d}")
        if len(comp['solutions']) > 10:
            print(f"    ... ({len(comp['solutions']) - 10} altre soluzioni)")
    
    # 3. Analisi reticolo di Hurwitz
    print("\n\n3. RETICOLO DI HURWITZ")
    print("-" * 70)
    
    norm_histogram, target_info = analyze_hurwitz_orbits(max_norm=50)
    
    print("\nUnità di Hurwitz (24 totali):")
    units = hurwitz_units()
    print(f"  Totale: {len(units)} unità")
    
    print("\nNorme di interesse nel reticolo di Hurwitz:")
    for norm in [1, 3, 21]:
        info = target_info[norm]
        print(f"\n  Norma {norm}:")
        print(f"    Cardinalità: {info['count']}")
        if info['examples']:
            print(f"    Esempi:")
            for ex in info['examples'][:5]:
                t, x, y, z, qtype = ex
                print(f"      ({t}, {x}, {y}, {z}) [{qtype}]")
        else:
            print(f"    Nessun esempio trovato")
    
    print("\n\nDistribuzione delle norme nel reticolo di Hurwitz (≤ 50):")
    for norm in sorted(norm_histogram.keys())[:20]:
        print(f"  Norma {norm:2d}: {norm_histogram[norm]:3d} elementi")
    
    # 4. Riepilogo
    print("\n\n" + "=" * 70)
    print("RIEPILOGO")
    print("=" * 70)
    
    if not solutions:
        print("\n✗ Le norme 1, 3, 21 non possono essere completate simultaneamente")
        print("  su una stessa sfera quaternionica S³ con coordinate intere.")
        print("\nQuesto è un risultato significativo: indica che la struttura Krist-All")
        print("non emerge semplicemente da una sezione di S³, ma richiede una")
        print("struttura più complessa.")
    else:
        print(f"\n✓ Trovate {len(solutions)} soluzioni simultanee.")
    
    print("\nIl prossimo passo è investigare se:")
    print("  - I punti Krist-All appartengono a orbite del reticolo di Hurwitz")
    print("  - La struttura 64×6×6×6×5 può essere interpretata come prodotto")
    print("    di orbite quaternioniche × settori di simmetria")


if __name__ == "__main__":
    main()
