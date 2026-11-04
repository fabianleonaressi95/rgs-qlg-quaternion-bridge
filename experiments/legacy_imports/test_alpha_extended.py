#!/usr/bin/env python3
"""
Test esteso per relazioni geometriche con 1/(2α).
Usa tolleranze più ampie per cercare coincidenze interessanti.
"""

import numpy as np
from math import sqrt, pi, e
from itertools import product, combinations


def compute_target():
    """Calcola il valore target: 1/(2α)."""
    alpha = 1/137.035999084
    return 1 / (2 * alpha)


def test_approximate_relations(target, tolerance_percent=1.0):
    """
    Testa relazioni approssimate con tolleranza percentuale.
    """
    constants = {
        '2/3': 2/3,
        'sqrt(3)/2': sqrt(3)/2,
        'sqrt(2)': sqrt(2),
        'sqrt(5)': sqrt(5),
        'pi': pi,
        'e': e,
        'pi*sqrt(3)/2': pi * sqrt(3) / 2,
    }
    
    names = list(constants.keys())
    values = list(constants.values())
    
    close_matches = []
    
    # Testa prodotti semplici
    for i in range(len(names)):
        for j in range(i, len(names)):
            for power_a in [1, 2, 3, 1/2, 3/2]:
                for power_b in [1, 2, 3, 1/2, 3/2]:
                    result = values[i] ** power_a * values[j] ** power_b
                    
                    relative_error = abs(result - target) / target
                    if relative_error < tolerance_percent / 100:
                        close_matches.append({
                            'formula': f"{names[i]}^{power_a} * {names[j]}^{power_b}",
                            'result': result,
                            'relative_error': relative_error * 100
                        })
    
    # Testa rapporti
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j:
                for power_a in [1, 2]:
                    for power_b in [1, 2]:
                        result = (values[i] ** power_a) / (values[j] ** power_b)
                        
                        relative_error = abs(result - target) / target
                        if relative_error < tolerance_percent / 100:
                            close_matches.append({
                                'formula': f"{names[i]}^{power_a} / {names[j]}^{power_b}",
                                'result': result,
                                'relative_error': relative_error * 100
                            })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:15]


def test_kristall_ratios(target, tolerance_percent=5.0):
    """
    Testa rapporti che coinvolgono le norme Krist-All (1, 3, 21).
    """
    norms = [1, 3, 21]
    close_matches = []
    
    geometric = {
        '2/3': 2/3,
        'sqrt(3)/2': sqrt(3)/2,
        'pi': pi,
    }
    
    # Testa combinazioni del tipo: (norm_a / norm_b) * geometric
    for i in range(len(norms)):
        for j in range(len(norms)):
            if i != j:
                ratio = norms[i] / norms[j]
                sqrt_ratio = sqrt(ratio)
                
                for geo_name, geo_value in geometric.items():
                    result = sqrt_ratio * geo_value
                    
                    relative_error = abs(result - target) / target
                    if relative_error < tolerance_percent / 100:
                        close_matches.append({
                            'formula': f"sqrt({norms[i]}/{norms[j]}) * {geo_name}",
                            'result': result,
                            'relative_error': relative_error * 100
                        })
                    
                    result = sqrt_ratio / geo_value
                    relative_error = abs(result - target) / target
                    if relative_error < tolerance_percent / 100:
                        close_matches.append({
                            'formula': f"sqrt({norms[i]}/{norms[j]}) / {geo_name}",
                            'result': result,
                            'relative_error': relative_error * 100
                        })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_hurwitz_combinations(target, tolerance_percent=10.0):
    """
    Testa combinazioni con cardinalità Hurwitz (24, 96, 768).
    """
    card_values = [24, 96, 768]
    close_matches = []
    
    # Testa combinazioni lineari con potenze
    for i in range(len(card_values)):
        for power in [1/2, 1, 3/2, 2]:
            result = card_values[i] ** power
            
            relative_error = abs(result - target) / target
            if relative_error < tolerance_percent / 100:
                close_matches.append({
                    'formula': f"{card_values[i]}^{power}",
                    'result': result,
                    'relative_error': relative_error * 100
                })
    
    # Testa rapporti
    for i in range(len(card_values)):
        for j in range(len(card_values)):
            if i != j:
                result = card_values[i] / card_values[j]
                
                relative_error = abs(result - target) / target
                if relative_error < tolerance_percent / 100:
                    close_matches.append({
                        'formula': f"{card_values[i]} / {card_values[j]}",
                        'result': result,
                        'relative_error': relative_error * 100
                    })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_combined_geometric(target):
    """
    Testa combinazioni miste: geometriche * cardinalità.
    """
    close_matches = []
    
    geometric = [sqrt(3)/2, 2/3, pi, e]
    cardinalities = [24, 96, 768]
    
    for geo in geometric:
        for card in cardinalities:
            for power_geo in [1, 2]:
                for power_card in [1/2, 1, 3/2]:
                    result = (geo ** power_geo) * (card ** power_card)
                    
                    relative_error = abs(result - target) / target
                    if relative_error < 0.05:  # 5%
                        close_matches.append({
                            'formula': f"geo^{power_geo} * card^{power_card}",
                            'geo': geo,
                            'card': card,
                            'result': result,
                            'relative_error': relative_error * 100
                        })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def main():
    print("=" * 70)
    print("TEST ESTESO: RELAZIONI GEOMETRICHE PER 1/(2α)")
    print("=" * 70)
    
    target = compute_target()
    print(f"\nTarget: 1/(2α) = {target:.10f}")
    print(f"Valore α: {1/(2*target):.10f}")
    
    # Test 1: Relazioni approssimate
    print("\n\n1. RELAZIONI APPROSSIMATE (tolleranza 1%)")
    print("-" * 70)
    
    approx_matches = test_approximate_relations(target, tolerance_percent=1.0)
    
    if approx_matches:
        print(f"\nTrovate {len(approx_matches)} corrispondenze:")
        for match in approx_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.6f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata con tolleranza 1%.")
    
    # Test 2: Relazioni Krist-All
    print("\n\n2. RELAZIONI KRIST-ALL (tolleranza 5%)")
    print("-" * 70)
    
    kristall_matches = test_kristall_ratios(target, tolerance_percent=5.0)
    
    if kristall_matches:
        print(f"\nTrovate {len(kristall_matches)} corrispondenze:")
        for match in kristall_matches:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.6f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 3: Relazioni Hurwitz
    print("\n\n3. RELAZIONI HURWITZ (tolleranza 10%)")
    print("-" * 70)
    
    hurwitz_matches = test_hurwitz_combinations(target, tolerance_percent=10.0)
    
    if hurwitz_matches:
        print(f"\nTrovate {len(hurwitz_matches)} corrispondenze:")
        for match in hurwitz_matches:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.6f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 4: Combinazioni miste
    print("\n\n4. COMBINAZIONI MISTE GEOMETRICHE × CARDINALITÀ")
    print("-" * 70)
    
    combined_matches = test_combined_geometric(target)
    
    if combined_matches:
        print(f"\nTrovate {len(combined_matches)} corrispondenze:")
        for match in combined_matches:
            print(f"  geo={match['geo']:.4f}, card={match['card']}")
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.6f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Conclusione
    print("\n\n" + "=" * 70)
    print("CONCLUSIONI")
    print("=" * 70)
    
    total_matches = len(approx_matches) + len(kristall_matches) + len(hurwitz_matches) + len(combined_matches)
    
    if total_matches == 0:
        print("\n✗ Nessuna relazione geometrica trovata nemmeno con tolleranze ampie.")
        print("\nQuesto conferma che:")
        print("  1. 1/(2α) non emerge dalle costanti geometriche pure")
        print("  2. La scala elettromagnetica e²/(ε₀c) è necessaria")
        print("  3. La struttura Krist-All/Hurwitz fornisce un 'alfabeto'")
        print("     geometrico ma non genera direttamente α")
    else:
        print(f"\nTrovate {total_matches} corrispondenze approssimate.")
        print("Tuttavia, queste richiedono tolleranze significative")
        print("e potrebbero essere coincidenze numeriche.")


if __name__ == "__main__":
    main()
