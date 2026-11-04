#!/usr/bin/env python3
"""
Investigazione se α (o 1/α) può emergere dalla geometria
del reticolo di Hurwitz o del 24-cell come numero puro.
"""

import numpy as np
from math import sqrt, pi, e
from itertools import product, combinations
from fractions import Fraction


def alpha_target():
    """Restituisce il valore target per 1/α."""
    alpha = 1/137.035999084
    return 1 / alpha


def hurwitz_geometric_constants():
    """Restituisce costanti geometriche del reticolo di Hurwitz e 24-cell."""
    return {
        # Cardinalità
        '24': 24,
        '96': 96,
        '768': 768,
        '4': 4,
        '32': 32,
        # Rapporti cardinalità
        '96/24': 96/24,
        '768/24': 768/24,
        '768/96': 768/96,
        # Geometria 24-cell
        '24_cell_volume': 2,
        '24_cell_inradius': sqrt(2)/2,
        'packing_density': pi**2/16,
        # Costanti fondamentali
        'pi': pi,
        'pi^2': pi**2,
        'sqrt2': sqrt(2),
        'sqrt3': sqrt(3),
        'sqrt5': sqrt(5),
        'e': e,
        # Rapporti geometrici
        'sqrt3/2': sqrt(3)/2,
        '2/3': 2/3,
        'pi*sqrt3/2': pi*sqrt(3)/2,
    }


def test_linear_combinations_1_over_alpha(target, constants, tolerance=0.001):
    """
    Testa combinazioni lineari delle costanti per approssimare 1/α.
    Cerca: a*const1 + b*const2 ≈ target
    """
    names = list(constants.keys())
    values = list(constants.values())
    
    close_matches = []
    
    # Coefficienti interi da -100 a 100
    for a in range(-100, 101):
        for b in range(-100, 101):
            if a == 0 and b == 0:
                continue
            
            for i in range(len(names)):
                for j in range(len(names)):
                    result = a * values[i] + b * values[j]
                    
                    relative_error = abs(result - target) / target
                    if relative_error < tolerance:
                        close_matches.append({
                            'formula': f"{a}*{names[i]} + {b}*{names[j]}",
                            'result': result,
                            'error': abs(result - target),
                            'relative_error': relative_error * 100
                        })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_products_1_over_alpha(target, constants, tolerance=0.001):
    """
    Testa prodotti di costanti per approssimare 1/α.
    Cerca: const1^a * const2^b ≈ target
    """
    names = list(constants.keys())
    values = list(constants.values())
    
    close_matches = []
    
    # Esponenti razionali semplici
    exponents = [Fraction(n, d) for n in range(-3, 4) for d in [1, 2, 3, 4]
                 if not (d > 1 and n == 0)]
    
    for i in range(len(names)):
        for j in range(len(names)):
            for exp_a in exponents:
                for exp_b in exponents:
                    result = values[i] ** float(exp_a) * values[j] ** float(exp_b)
                    
                    if result > 0 and result < 1000:  # Range ragionevole
                        relative_error = abs(result - target) / target
                        if relative_error < tolerance:
                            close_matches.append({
                                'formula': f"{names[i]}^{exp_a} * {names[j]}^{exp_b}",
                                'result': result,
                                'error': abs(result - target),
                                'relative_error': relative_error * 100
                            })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_ratios_1_over_alpha(target, constants, tolerance=0.001):
    """
    Testa rapporti di costanti per approssimare 1/α.
    Cerca: const1 / const2 ≈ target
    """
    names = list(constants.keys())
    values = list(constants.values())
    
    close_matches = []
    
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j and values[j] != 0:
                result = values[i] / values[j]
                
                if result > 0 and result < 1000:
                    relative_error = abs(result - target) / target
                    if relative_error < tolerance:
                        close_matches.append({
                            'formula': f"{names[i]} / {names[j]}",
                            'result': result,
                            'error': abs(result - target),
                            'relative_error': relative_error * 100
                        })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_hurwitz_specific_combinations(target):
    """
    Testa combinazioni specifiche legate al reticolo di Hurwitz.
    """
    close_matches = []
    
    # Cardinalità
    card = [24, 96, 768]
    
    # Testa combinazioni del tipo: a*card1 + b*card2 + c*pi + d
    for a in range(-5, 6):
        for b in range(-5, 6):
            for c in range(-5, 6):
                for d in range(-10, 11):
                    result = a*card[0] + b*card[1] + c*pi + d
                    
                    relative_error = abs(result - target) / target
                    if relative_error < 0.01:  # 1%
                        close_matches.append({
                            'formula': f"{a}*24 + {b}*96 + {c}*π + {d}",
                            'result': result,
                            'relative_error': relative_error * 100
                        })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_pi_based_combinations(target):
    """
    Testa combinazioni basate su π, dato che π appare spesso in fisica.
    """
    close_matches = []
    
    # Testa: a*pi + b*pi^2 + c
    for a in range(-50, 51):
        for b in range(-10, 11):
            for c in range(-100, 101):
                result = a*pi + b*pi**2 + c
                
                relative_error = abs(result - target) / target
                if relative_error < 0.001:  # 0.1%
                    close_matches.append({
                        'formula': f"{a}*π + {b}*π² + {c}",
                        'result': result,
                        'relative_error': relative_error * 100
                    })
    
    return sorted(close_matches, key=lambda x: x['relative_error'])[:10]


def test_famous_approximations(target):
    """
    Testa approssimazioni famose di 1/α.
    """
    approximations = {
        '137': 137,
        '137.036': 137.036,
        '4π^3 + π^2 + π': 4*pi**3 + pi**2 + pi,
        'π^3 + π^2': pi**3 + pi**2,
        '2π^3': 2*pi**3,
        'π^4': pi**4,
        'π^3 + 29': pi**3 + 29,
        '29*π': 29*pi,
        '44*π': 44*pi,
    }
    
    print("Approssimazioni famose di 1/α:")
    for name, value in approximations.items():
        error = abs(value - target)
        relative_error = error / target * 100
        print(f"  {name:20s} = {value:.10f}")
        print(f"    Errore: {error:.10f}, Errore relativo: {relative_error:.6f}%")
    
    return approximations


def investigate_137_geometry():
    """
    Investiga se 137 (o 137.036) ha significato geometrico.
    """
    print("\n" + "=" * 70)
    print("INVESTIGAZIONE GEOMETRICA DI 137")
    print("=" * 70)
    
    # 137 è un numero primo
    print(f"\n137 è un numero primo: ✓")
    
    # Fattorizzazione di 137.036
    print(f"\nFattorizzazione di interi vicini:")
    for n in [136, 137, 138]:
        factors = []
        temp = n
        for i in range(2, int(sqrt(n)) + 1):
            while temp % i == 0:
                factors.append(i)
                temp //= i
        if temp > 1:
            factors.append(temp)
        print(f"  {n} = {' × '.join(map(str, factors))}")
    
    # Connessioni geometriche
    print(f"\nConnessioni geometriche potenziali:")
    print(f"  24 × 5 + 37 = {24*5 + 37}")
    print(f"  96 + 41 = {96 + 41}")
    print(f"  768 / 5.6 ≈ {768/5.6:.2f}")
    
    # Relazione con π
    print(f"\nRelazioni con π:")
    print(f"  137 / π = {137/pi:.10f}")
    print(f"  π × 43.6 ≈ {pi*43.6:.2f}")
    print(f"  π^2 × 13.9 ≈ {pi**2*13.9:.2f}")
    
    # Relazione con e
    print(f"\nRelazioni con e:")
    print(f"  137 / e = {137/e:.10f}")
    print(f"  e^5 ≈ {e**5:.2f}")


def main():
    print("=" * 70)
    print("INVESTIGAZIONE: α DALLA GEOMETRIA")
    print("=" * 70)
    
    target = alpha_target()
    print(f"\nTarget: 1/α = {target:.10f}")
    print(f"Valore α: {1/target:.10f}")
    
    constants = hurwitz_geometric_constants()
    
    # Test 1: Approssimazioni famose
    print("\n\n1. APPROSSIMAZIONI FAMOSE")
    print("-" * 70)
    test_famous_approximations(target)
    
    # Test 2: Combinazioni lineari
    print("\n\n2. COMBINAZIONI LINEARI")
    print("-" * 70)
    linear_matches = test_linear_combinations_1_over_alpha(target, constants, tolerance=0.01)
    
    if linear_matches:
        print(f"\nTrovate {len(linear_matches)} corrispondenze:")
        for match in linear_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.10f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 3: Prodotti
    print("\n\n3. PRODOTTI DI COSTANTI")
    print("-" * 70)
    product_matches = test_products_1_over_alpha(target, constants, tolerance=0.01)
    
    if product_matches:
        print(f"\nTrovate {len(product_matches)} corrispondenze:")
        for match in product_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.10f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 4: Rapporti
    print("\n\n4. RAPPORTI DI COSTANTI")
    print("-" * 70)
    ratio_matches = test_ratios_1_over_alpha(target, constants, tolerance=0.01)
    
    if ratio_matches:
        print(f"\nTrovate {len(ratio_matches)} corrispondenze:")
        for match in ratio_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.10f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 5: Combinazioni Hurwitz specifiche
    print("\n\n5. COMBINAZIONI HURWITZ SPECIFICHE")
    print("-" * 70)
    hurwitz_matches = test_hurwitz_specific_combinations(target)
    
    if hurwitz_matches:
        print(f"\nTrovate {len(hurwitz_matches)} corrispondenze:")
        for match in hurwitz_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.10f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Test 6: Combinazioni basate su π
    print("\n\n6. COMBINAZIONI BASATE SU π")
    print("-" * 70)
    pi_matches = test_pi_based_combinations(target)
    
    if pi_matches:
        print(f"\nTrovate {len(pi_matches)} corrispondenze:")
        for match in pi_matches[:5]:
            print(f"  {match['formula']}")
            print(f"    Risultato: {match['result']:.10f}")
            print(f"    Errore relativo: {match['relative_error']:.4f}%")
    else:
        print("\nNessuna corrispondenza trovata.")
    
    # Investigazione geometrica di 137
    investigate_137_geometry()
    
    # Conclusione
    print("\n\n" + "=" * 70)
    print("CONCLUSIONI")
    print("=" * 70)
    
    total_matches = len(linear_matches) + len(product_matches) + len(ratio_matches) + len(hurwitz_matches) + len(pi_matches)
    
    if total_matches == 0:
        print("\n✗ Nessuna relazione geometrica semplice trovata per 1/α.")
        print("\nQuesto suggerisce che:")
        print("  1. α = 1/137.036... non emerge da combinazioni semplici")
        print("     delle costanti geometriche del reticolo di Hurwitz")
        print("  2. Potrebbe richiedere una struttura più complessa")
        print("  3. O potrebbe essere una costante fondamentale non derivabile")
    else:
        print(f"\nTrovate {total_matches} corrispondenze approssimate.")
        print("Queste richiedono ulteriore analisi per determinare")
        print("se sono significative o coincidenze numeriche.")


if __name__ == "__main__":
    main()
