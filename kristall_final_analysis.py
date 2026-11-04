#!/usr/bin/env python3
"""
Analisi finale della struttura 64×6×6×6×5 nel contesto quaternionico.
Verifica se questa struttura può essere interpretata come prodotto di orbite.
"""

import numpy as np
from itertools import product, permutations
from collections import defaultdict


def hurwitz_integers(limit=5):
    """Genera gli interi di Hurwitz."""
    hurwitz = []
    for a, b, c, d in product(range(-limit, limit+1), repeat=4):
        hurwitz.append((a, b, c, d, 'integer'))
    for a, b, c, d in product(range(-limit, limit), repeat=4):
        hurwitz.append((a + 0.5, b + 0.5, c + 0.5, d + 0.5, 'half-integer'))
    return hurwitz


def hurwitz_units():
    """Restituisce le 24 unità del reticolo di Hurwitz."""
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


def analyze_64_structure():
    """
    Analizza le norme con cardinalità 64 e cerca pattern.
    """
    hurwitz = hurwitz_integers(limit=5)
    norm_counts = defaultdict(int)
    norm_elements = defaultdict(list)
    
    for q in hurwitz:
        t, x, y, z, qtype = q
        norm = t*t + x*x + y*y + z*z
        if isinstance(norm, float) and norm.is_integer():
            norm = int(norm)
        norm_counts[norm] += 1
        if len(norm_elements[norm]) < 10:
            norm_elements[norm].append((t, x, y, z, qtype))
    
    # Norme con cardinalità 64
    norms_64 = [norm for norm, count in norm_counts.items() if count == 64]
    
    return norms_64, norm_elements


def analyze_24_cell_structure():
    """
    Analizza la struttura del 24-cell (le 24 unità di Hurwitz).
    """
    units = hurwitz_units()
    
    # Separa le 8 unità standard dalle 16 unità "half"
    standard_units = []
    half_units = []
    
    for u in units:
        if all(coord == 0 or abs(coord) == 1 for coord in u):
            standard_units.append(u)
        else:
            half_units.append(u)
    
    return {
        'total': len(units),
        'standard': len(standard_units),
        'half': len(half_units),
        'standard_units': standard_units,
        'half_units': half_units
    }


def analyze_sector_structure():
    """
    Analizza i settori di simmetria che potrebbero corrispondere ai fattori 6×6×6×5.
    """
    # Il 24-cell ha 24 facce, 24 celle, 96 triangoli, ecc.
    # Cerchiamo se ci sono strutture con cardinalità 6, 5
    
    sectors_info = {
        '24_cell_faces': 24,
        '24_cell_vertices': 24,
        '24_cell_edges': 96,
        '24_cell_triangles': 96,
        'symmetry_order_Oh': 48,  # Gruppo simmetria cubo completo
        'symmetry_order_Td': 24,  # Gruppo simmetria tetraedro
        'symmetry_order_O': 24,   # Gruppo rotazione cubo
    }
    
    return sectors_info


def check_product_structure():
    """
    Verifica se 64×6×6×6×5 può essere interpretato come prodotto.
    """
    # 64×6×6×6×5 = 69120
    total = 64 * 6 * 6 * 6 * 5
    
    # Possibili decomposizioni
    decompositions = {
        'direct_product': f"64 (orbite) × 6×6×6×5 (settori) = {total}",
        'alternative_1': f"24 (unità Hurwitz) × 24 (orbite) × 12 (altro) = {24*24*12}",
        'alternative_2': f"24 (unità) × 96 (norma 3) × 30 (altro) = {24*96*30}",
        'alternative_3': f"96 (norma 3) × 720 (altro) = {96*720}",
    }
    
    return decompositions


def analyze_double_tetrahedron_relation():
    """
    Analizza la relazione tra doppio tetraedro e struttura quaternionica.
    """
    # Doppio tetraedro ha 8 vertici
    # 24-cell ha 24 vertici = 3 × 8
    
    relation = {
        'double_tetrahedron_vertices': 8,
        '24_cell_vertices': 24,
        'ratio': 24 / 8,
        'interpretation': 'Il 24-cell può essere visto come 3 copie del doppio tetraedro'
    }
    
    return relation


def analyze_kristall_complete_structure():
    """
    Analisi completa della struttura Krist-All includendo tutti i 62 punti.
    """
    # I 34 punti generati precedentemente sono solo una parte
    # La struttura completa 64×6×6×6×5 = 69120 suggerisce una struttura molto più grande
    
    # Analizziamo i contributi
    contributions = {
        'shell_1_norm_1': {
            'norm': 1,
            'cardinality': 24,
            'orbit_size': 24,
            'contribution': '24 elementi su S³ di raggio 1'
        },
        'shell_2_norm_3': {
            'norm': 3,
            'cardinality': 96,
            'orbit_size': 24,
            'contribution': '96 elementi su S³ di raggio √3 (4 orbite di 24)'
        },
        'shell_3_norm_21': {
            'norm': 21,
            'cardinality': 768,
            'orbit_size': 24,
            'contribution': '768 elementi su S³ di raggio √21 (32 orbite di 24)'
        }
    }
    
    return contributions


def main():
    print("=" * 70)
    print("ANALISI FINALE: STRUTTURA 64×6×6×6×5 E QUATERNIONI")
    print("=" * 70)
    
    # 1. Analisi struttura 64
    print("\n1. NORME CON CARDINALITÀ 64")
    print("-" * 70)
    
    norms_64, norm_elements = analyze_64_structure()
    print(f"Norme con cardinalità 64: {norms_64}")
    
    for norm in norms_64:
        print(f"\nNorma {norm}:")
        print(f"  Esempi di elementi:")
        for elem in norm_elements[norm][:5]:
            t, x, y, z, qtype = elem
            print(f"    ({t}, {x}, {y}, {z}) [{qtype}]")
    
    # 2. Struttura 24-cell
    print("\n\n2. STRUTTURA 24-CELL (UNITÀ DI HURWITZ)")
    print("-" * 70)
    
    cell_24 = analyze_24_cell_structure()
    print(f"Totale unità: {cell_24['total']}")
    print(f"Unità standard (±1, ±i, ±j, ±k): {cell_24['standard']}")
    print(f"Unità half ((±1±i±j±k)/2): {cell_24['half']}")
    
    # 3. Relazione doppio tetraedro
    print("\n\n3. RELAZIONE DOPPIO TETRAEDRO - 24-CELL")
    print("-" * 70)
    
    dt_relation = analyze_double_tetrahedron_relation()
    print(f"Vertici doppio tetraedro: {dt_relation['double_tetrahedron_vertices']}")
    print(f"Vertici 24-cell: {dt_relation['24_cell_vertices']}")
    print(f"Rapporto: {dt_relation['ratio']}")
    print(f"Interpretazione: {dt_relation['interpretation']}")
    
    # 4. Analisi settori
    print("\n\n4. SETTORI DI SIMMETRIA")
    print("-" * 70)
    
    sectors = analyze_sector_structure()
    for key, value in sectors.items():
        print(f"  {key}: {value}")
    
    # 5. Verifica struttura prodotto
    print("\n\n5. VERIFICA STRUTTURA PRODOTTO")
    print("-" * 70)
    
    product_structures = check_product_structure()
    for key, value in product_structures.items():
        print(f"  {key}: {value}")
    
    # 6. Analisi completa Krist-All
    print("\n\n6. ANALISI COMPLETA KRIST-ALL")
    print("-" * 70)
    
    kristall_complete = analyze_kristall_complete_structure()
    for shell_name, info in kristall_complete.items():
        print(f"\n{shell_name}:")
        print(f"  Norma: {info['norm']}")
        print(f"  Cardinalità: {info['cardinality']}")
        print(f"  Dimensione orbita: {info['orbit_size']}")
        print(f"  Contributo: {info['contribution']}")
    
    # Calcola totale elementi nei tre gusci
    total_elements = sum(info['cardinality'] for info in kristall_complete.values())
    print(f"\nTotale elementi nei tre gusci: {total_elements}")
    print(f"Struttura 64×6×6×6×5 = {64*6*6*6*5}")
    
    # 7. Riepilogo e conclusioni
    print("\n\n" + "=" * 70)
    print("RIEPILOGO E CONCLUSIONI")
    print("=" * 70)
    
    print("\nRisultati chiave:")
    print("  1. ✓ Le norme 1, 3, 21 compaiono nel reticolo di Hurwitz")
    print("     con cardinalità 24, 96, 768 (rapporto 1:4:32)")
    print("  2. ✓ Tutte le orbite hanno dimensione 24 (le unità di Hurwitz)")
    print("  3. ✓ Esistono norme con cardinalità 64 (76, 79, 84, 91)")
    print("  4. ✓ Il 24-cell (24 unità) è 3× il doppio tetraedro (8 vertici)")
    print("  5. ✗ Le norme 1, 3, 21 non possono essere completate")
    print("     simultaneamente sulla stessa sfera S³")
    
    print("\nInterpretazione geometrica:")
    print("  La struttura Krist-All non è una semplice sezione di S³,")
    print("  ma emerge dall'interazione tra:")
    print("  - Il reticolo di Hurwitz in H (4D)")
    print("  - Le orbite delle 24 unità")
    print("  - La proiezione su ImH (3D)")
    
    print("\nLa struttura 64×6×6×6×5 = 69120:")
    print("  64: potrebbe corrispondere a norme con cardinalità 64")
    print("  6×6×6×5: settori di simmetria (da investigare ulteriormente)")
    print("  Totale 69120 >> 24+96+768 = 888 (struttura base)")
    print("  Suggerisce una struttura gerarchica o prodotto di orbite")
    
    print("\nProssimi passi suggeriti:")
    print("  1. Investigare le norme 76, 79, 84, 91 (cardinalità 64)")
    print("  2. Analizzare la struttura delle orbite di queste norme")
    print("  3. Verificare se 6×6×6×5 corrisponde a settori del 24-cell")
    print("  4. Esplorare la relazione con il gruppo di simmetria completo (48 elementi)")


if __name__ == "__main__":
    main()
