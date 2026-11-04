#!/usr/bin/env python3
"""
Test combinatorio-geometrico: derivare gli esponenti di Wyler
dagli invarianti canonici di Krist-All.

Target: (19,-7,1,11) per (2,3,5,π)
Equivalente: (8, 9/4, 7/4, 1/4, 11/4) nella forma con 2/3

Invarianti Krist-All da testare:
- Norme: 1, 3, 21
- Pesi: 1/2, 1/4 (se presente)
- M₂ = 26I
- c₄
- Rapporto 2/3
- Cardinalità: 62, 34, ecc.
- Struttura geometrica (tetraedro, ottaedro, cubo)
"""

from math import sqrt, pi
from collections import defaultdict
from itertools import permutations


def kristall_points():
    """Genera i punti Krist-All."""
    points = []
    
    # Guscio (1,0,0): 6 punti (tutte le permutazioni con segni)
    for i in range(3):
        for sign in [1, -1]:
            coord = [0, 0, 0]
            coord[i] = sign
            points.append((coord[0], coord[1], coord[2], 1))
    
    # Guscio (1,1,1): 8 punti (tutti i segni)
    for x in [1, -1]:
        for y in [1, -1]:
            for z in [1, -1]:
                points.append((x, y, z, 1/2))
    
    # Guscio (1,2,4): 24 punti (tutte le permutazioni con segni)
    for perm in permutations([1, 2, 4]):
        for signs in [(1,1,1), (1,1,-1), (1,-1,1), (1,-1,-1),
                      (-1,1,1), (-1,1,-1), (-1,-1,1), (-1,-1,-1)]:
            points.append((signs[0]*perm[0], signs[1]*perm[1], signs[2]*perm[2], 1/4))
    
    return points


def analyze_norms(points):
    """Analizza le norme Krist-All."""
    print("=" * 70)
    print("ANALISI NORME KRIST-ALL")
    print("=" * 70)
    
    norms = {}
    for x, y, z, w in points:
        norm_sq = x*x + y*y + z*z
        if norm_sq not in norms:
            norms[norm_sq] = []
        norms[norm_sq].append((x, y, z, w))
    
    print(f"\nNorme trovate:")
    for norm in sorted(norms.keys()):
        count = len(norms[norm])
        print(f"  r² = {norm}: {count} punti")
        
        # Fattorizzazione
        factors = []
        n = norm
        for p in [2, 3, 5, 7]:
            exp = 0
            while n % p == 0 and n > 1:
                exp += 1
                n //= p
            if exp > 0:
                factors.append(f"{p}^{exp}")
        if n > 1:
            factors.append(str(n))
        
        print(f"    Fattorizzazione: {' × '.join(factors) if factors else str(norm)}")
    
    return norms


def analyze_exponents_from_geometry(points):
    """Tenta di derivare esponenti da invarianti geometrici."""
    print(f"\n{'='*70}")
    print("DERIVAZIONE ESPONENTI DA GEOMETRIA")
    print(f"{'='*70}")
    
    # ============================================================
    # 1. Analisi delle norme
    # ============================================================
    print(f"\n1. Analisi norme (1, 3, 21):")
    
    # Norme: 1 = 1, 3 = 3¹, 21 = 3¹·7¹
    # Questi suggeriscono esponenti per 3: 0, 1, 1
    # Ma non c'è uno schema chiaro per 2, 5, π
    
    # ============================================================
    # 2. Analisi dei pesi
    # ============================================================
    print(f"\n2. Analisi pesi (1, 1/2, 1/4):")
    
    weights = defaultdict(int)
    for x, y, z, w in points:
        weights[w] += 1
    
    for w in sorted(weights.keys()):
        print(f"  Peso {w}: {weights[w]} punti")
    
    # I pesi 1, 1/2, 1/4 sono potenze di 2
    # Suggeriscono esponenti per 2: 0, -1, -2
    # Ma non corrispondono a 19/4
    
    # ============================================================
    # 3. Rapporto 2/3
    # ============================================================
    print(f"\n3. Rapporto 2/3 (doppio tetraedro/cubo):")
    print(f"  2/3 = 2¹ × 3⁻¹")
    print(f"  Questo suggerisce esponenti (1, -1) per (2, 3)")
    print(f"  Wyler richiede (19/4, -7/4) = (4.75, -1.75)")
    print(f"  Ratio mancano: 19/4 - 1 = 15/4 = 3.75 per 2")
    print(f"               -7/4 - (-1) = -3/4 = -0.75 per 3")
    
    # ============================================================
    # 4. Cardinalità
    # ============================================================
    print(f"\n4. Cardinalità (34 punti pesati, 62 totali):")
    
    total_points = len(points)
    weighted_sum = sum(w for _, _, _, w in points)
    
    print(f"  Totale punti: {total_points}")
    print(f"  Somma pesi: {weighted_sum}")
    
    # Fattorizzazione
    for n in [total_points, weighted_sum]:
        factors = []
        temp = n
        for p in [2, 3, 5, 7]:
            exp = 0
            while temp % p == 0 and temp > 1:
                exp += 1
                temp //= p
            if exp > 0:
                factors.append(f"{p}^{exp}")
        if temp > 1:
            factors.append(str(temp))
        print(f"  {n} = {' × '.join(factors) if factors else str(n)}")
    
    # ============================================================
    # 5. 21 = 3·7
    # ============================================================
    print(f"\n5. Analisi 21 = 3·7:")
    print(f"  Krist-All introduce 7")
    print(f"  Wyler non contiene 7 (esponente 0)")
    print(f"  Questo suggerisce che 21 deve CANCELLARSI o appartenere a altro settore")
    
    # Se 21 dovesse generare l'esponente 7/4:
    # 21/12 = 7/4, quindi 7 = 21 × 4/12
    # Ma 12 non emerge chiaramente
    
    # ============================================================
    # 6. Struttura dimensionale 5?
    # ============================================================
    print(f"\n6. Ricerca struttura dimensionale 5:")
    print(f"  5! = 120 = 2³·3·5")
    print(f"  Krist-All alfabeto: {2, 3, 5}")
    print(f"  ✓ I primi sono presenti")
    print(f"  ✗ Manca la struttura geometrica di dimensione 5")
    
    # Cerca se c'è un invariante che produce 5
    # Cardinalità gusci: 4, 8, 24
    # Rapporti: 24/4 = 6, 24/8 = 3, 8/4 = 2
    # Nessuno è 5
    
    print(f"  Cardinalità gusci: 4, 8, 24")
    print(f"  Rapporti: 24/4=6, 24/8=3, 8/4=2")
    print(f"  Nessun rapporto uguale a 5")
    
    # ============================================================
    # 7. Verifica esponenti target
    # ============================================================
    print(f"\n{'='*70}")
    print("VERIFICA ESPONENTI TARGET")
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
    
    # Verifica se qualche combinazione di invarianti produce questi
    print(f"\nTentativo combinazioni invarianti:")
    
    # Combinazioni possibili:
    # - Cardinalità: 34, 62
    # - Norme: 1, 3, 21
    # - Rapporti: 2/3
    # - Pesi: 1/2, 1/4
    
    print(f"\nLog(2) esponenti da cardinalità:")
    for n in [34, 62]:
        log2_n = np.log2(n)
        print(f"  log₂({n}) = {log2_n:.4f}")
        print(f"  vs target 19/4 = {19/4:.4f}")
        print(f"  differenza: {log2_n - 19/4:.4f}")
    
    print(f"\nLog(3) esponenti da norme:")
    for n in [1, 3, 21]:
        if n > 0:
            log3_n = np.log(n) / np.log(3)
            print(f"  log₃({n}) = {log3_n:.4f}")
            print(f"  vs target -7/4 = {-7/4:.4f}")
            print(f"  differenza: {log3_n - (-7/4):.4f}")
    
    print(f"\nEsponenti da rapporti:")
    print(f"  2/3 → esponenti (1, -1)")
    print(f"  Target: (19/4, -7/4) = (4.75, -1.75)")
    print(f"  Differenza: (3.75, -0.75)")
    
    # ============================================================
    # 8. Conclusione
    # ============================================================
    print(f"\n{'='*70}")
    print("CONCLUSIONE")
    print(f"{'='*70}")
    
    print(f"\n✓ Alfabeto Krist-All ⊂ ingredienti Wyler")
    print(f"✗ Esponenti Wyler NON emergono da invarianti geometrici ovvi")
    print(f"✗ 7 compare in Krist-All ma non in Wyler")
    print(f"✗ Struttura dimensionale 5 non identificata")
    
    print(f"\nEsponenti derivabili:")
    print(f"  Da 2/3: (1, -1) per (2, 3)")
    print(f"  Da pesi: (0, -1, -2) per 2")
    print(f"  Da norme: (0, 1, 1) per 3")
    
    print(f"\nEsponenti target:")
    print(f"  (19/4, -7/4, 1/4, 11/4)")
    
    print(f"\nGap:")
    print(f"  Per 2: mancano 15/4 = 3.75")
    print(f"  Per 3: mancano -3/4 = -0.75")
    print(f"  Per 5: mancano 1/4 = 0.25")
    print(f"  Per π: mancano 11/4 = 2.75")


import numpy as np


def main():
    print("=" * 70)
    print("TEST COMBINATORIO-GEOMETRICO KRIST-ALL → WYLER")
    print("=" * 70)
    
    points = kristall_points()
    
    analyze_norms(points)
    analyze_exponents_from_geometry(points)
    
    print(f"\n{'='*70}")
    print("VERDETTO")
    print(f"{'='*70}")
    
    print(f"\nL'alfabeto Krist-All è sufficiente per RAPPRESENTARE Wyler,")
    print(f"ma gli invarianti geometrici attuali NON GENERANO gli esponenti.")
    print(f"\nServe un invariante più sofisticato che produca:")
    print(f"  - 7/4 per l'esponente di 2/3")
    print(f"  - 1/4 per l'esponente di 5")
    print(f"  - 11/4 per l'esponente di π")
    print(f"  - Cancellazione del 7 da 21")


if __name__ == "__main__":
    main()
