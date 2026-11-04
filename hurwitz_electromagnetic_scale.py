#!/usr/bin/env python3
"""
Investigazione se la scala elettromagnetica e²/(ε₀c) possa emergere
dalle proprietà geometriche del reticolo di Hurwitz o del 24-cell.
"""

import numpy as np
from math import sqrt, pi
from itertools import product


def physical_constants():
    """Restituisce le costanti fisiche rilevanti."""
    return {
        'e': 1.602176634e-19,  # Carica elementare (C)
        'epsilon_0': 8.8541878128e-12,  # Permittività del vuoto (F/m)
        'c': 299792458,  # Velocità della luce (m/s)
        'h': 6.62607015e-34,  # Costante di Planck (J·s)
        'alpha': 1/137.035999084,  # Costante di struttura fine
    }


def electromagnetic_scale():
    """Calcola e²/(ε₀c)."""
    phys = physical_constants()
    return phys['e']**2 / (phys['epsilon_0'] * phys['c'])


def hurwitz_lattice_density():
    """
    Calcola la densità del reticolo di Hurwitz.
    Il reticolo di Hurwitz è il reticolo più denso in 4D.
    """
    # Il volume fondamentale del reticolo di Hurwitz è 1
    # La densità di impacchettamento massima in 4D è π²/16 ≈ 0.61685
    packing_density = pi**2 / 16
    
    return {
        'fundamental_volume': 1,
        'packing_density': packing_density,
        'packing_density_decimal': float(packing_density)
    }


def hurwitz_minimal_distance():
    """
    Calcola la distanza minima tra punti nel reticolo di Hurwitz.
    """
    # Nel reticolo di Hurwitz, la distanza minima è 1
    # (tra punti interi adiacenti o mezzi-interi adiacenti)
    return 1.0


def hurwitz_covolume():
    """
    Calcola il covolume (volume della cella fondamentale) del reticolo di Hurwitz.
    """
    # Per il reticolo di Hurwitz, il covolume è 1
    return 1.0


def cell_24_properties():
    """
    Calcola proprietà geometriche del 24-cell.
    """
    # Il 24-cell in 4D ha:
    # - 24 vertici
    # - 96 facce (triangoli)
    # - 24 celle (ottaedri)
    # - 96 spigoli
    
    # Coordinate dei vertici del 24-cell (versione standardizzata)
    vertices = []
    
    # 8 vertici tipo (±1, 0, 0, 0) e permutazioni
    for i in range(4):
        for sign in [1, -1]:
            v = [0, 0, 0, 0]
            v[i] = sign
            vertices.append(tuple(v))
    
    # 16 vertici tipo (±1/2, ±1/2, ±1/2, ±1/2)
    for signs in product([1, -1], repeat=4):
        v = [signs[i]/2 for i in range(4)]
        vertices.append(tuple(v))
    
    # Calcola volume del 24-cell
    # Il volume del 24-cell di raggio R è 2R^4
    # Per il 24-cell unitario (vertici sulla sfera unitaria), R = 1
    volume = 2
    
    # Calcola raggio circoscritto
    # Per il 24-cell standard, i vertici sono sulla sfera di raggio 1
    circumradius = 1
    
    # Calcola raggio inscritto
    # Per il 24-cell, il raggio inscritto è √2/2
    inradius = sqrt(2) / 2
    
    return {
        'n_vertices': 24,
        'n_edges': 96,
        'n_faces': 96,
        'n_cells': 24,
        'volume': volume,
        'circumradius': circumradius,
        'inradius': inradius,
        'vertices': vertices
    }


def hurwitz_lattice_constant():
    """
    Calcola una costante del reticolo basata sulle proprietà geometriche.
    """
    # Costante di impacchettamento
    packing = pi**2 / 16
    
    # Volume della sfera 4D unitaria
    sphere_4d_volume = pi**2 / 2
    
    # Rapporto tra volume sfera e volume cella
    ratio = sphere_4d_volume / 2  # 2 è il volume del 24-cell
    
    return {
        'packing_constant': packing,
        'sphere_4d_volume': sphere_4d_volume,
        'volume_ratio': ratio
    }


def test_geometric_scales(em_scale):
    """
    Testa se la scala elettromagnetica può essere espressa come funzione
    di scale geometriche del reticolo di Hurwitz.
    """
    properties = hurwitz_lattice_constant()
    cell_24 = cell_24_properties()
    
    geometric_scales = {
        'packing_density': properties['packing_constant'],
        'volume_ratio': properties['volume_ratio'],
        '24_cell_volume': cell_24['volume'],
        '24_cell_inradius': cell_24['inradius'],
        'pi_squared': pi**2,
        'pi_squared_over_2': pi**2 / 2,
        'sqrt2': sqrt(2),
    }
    
    print("Scale geometriche disponibili:")
    for name, value in geometric_scales.items():
        print(f"  {name:20s} = {value:.10f}")
    
    print(f"\nScala elettromagnetica e²/(ε₀c): {em_scale:.6e} J·s")
    
    # Cerca relazioni del tipo: geometric_scale^n ≈ em_scale
    # Nota: questo è problematico perché em_scale ha dimensioni fisiche
    # mentre le scale geometriche sono adimensionali
    
    print("\nNota: Le scale geometriche sono adimensionali.")
    print("Per ottenere una scala con dimensioni [J·s], serve una scala fisica.")
    
    return geometric_scales


def investigate_dimensional_bridge():
    """
    Investiga possibili "ponti dimensionali" tra geometria e fisica.
    """
    print("\n" + "=" * 70)
    print("INVESTIGAZIONE PONTE DIMENSIONALE")
    print("=" * 70)
    
    phys = physical_constants()
    em_scale = electromagnetic_scale()
    
    print(f"\nScala elettromagnetica: {em_scale:.6e} J·s")
    print(f"Costante di Planck: {phys['h']:.6e} J·s")
    print(f"Rapporto: {em_scale/phys['h']:.10f}")
    print(f"Questo è ≈ 2α = {2*phys['alpha']:.10f}")
    
    print("\nScala naturale:")
    print("  Per ottenere una scala [J·s] da quantità adimensionali,")
    print("  serve moltiplicare per una scala con dimensioni [J·s].")
    print("  Candidate naturali:")
    print("    - h (costante di Planck)")
    print("    - ħ = h/(2π)")
    print("    - Compton wavelength related scales")
    
    # Testa: (geometric factor) × h ≈ e²/(ε₀c)
    geometric_factors = {
        '2α': 2 * phys['alpha'],
        'α': phys['alpha'],
        '1/(2π)': 1/(2*pi),
        'π²/16': pi**2/16,
        '√2/2': sqrt(2)/2,
    }
    
    print("\nTest: (fattore geometrico) × h ≈ e²/(ε₀c)")
    for name, factor in geometric_factors.items():
        result = factor * phys['h']
        ratio = result / em_scale
        print(f"  {name:10s} × h = {result:.6e} J·s")
        print(f"    Rapporto con e²/(ε₀c): {ratio:.10f}")


def investigate_hurwitz_scales():
    """
    Investiga scale naturali nel reticolo di Hurwitz.
    """
    print("\n" + "=" * 70)
    print("SCALE NATURALI NEL RETICOLO DI HURWITZ")
    print("=" * 70)
    
    # Distanze nel reticolo
    distances = {
        'minimal_distance': 1,
        'second_shell': sqrt(2),  # Distanza al secondo guscio
        'third_shell': sqrt(3),  # Distanza al terzo guscio
    }
    
    print("\nDistanze nel reticolo:")
    for name, dist in distances.items():
        print(f"  {name:20s} = {dist:.10f}")
    
    # Cardinalità dei gusci
    shell_cardinalities = {
        'norm_1': 24,
        'norm_2': 24,
        'norm_3': 96,
        'norm_4': 24,
        'norm_5': 144,
    }
    
    print("\nCardinalità dei gusci:")
    for name, card in shell_cardinalities.items():
        print(f"  {name:10s}: {card:4d} elementi")
    
    # Rapporti interessanti
    print("\nRapporti interessanti:")
    print(f"  96/24 = {96/24}")
    print(f"  144/24 = {144/24}")
    print(f"  768/24 = {768/24}")
    
    # Costante di struttura fine
    phys = physical_constants()
    print(f"\nConfronto con α = {phys['alpha']:.10f}")
    print(f"  1/137 ≈ {1/137:.10f}")
    print(f"  24/3289 ≈ {24/3289:.10f}")
    print(f"  π/432 ≈ {pi/432:.10f}")


def investigate_4d_geometry_to_3d_projection():
    """
    Investiga come la geometria 4D si proietta in 3D.
    """
    print("\n" + "=" * 70)
    print("PROIEZIONE 4D → 3D")
    print("=" * 70)
    
    cell_24 = cell_24_properties()
    
    print(f"\nProprietà 24-cell:")
    print(f"  Vertici: {cell_24['n_vertices']}")
    print(f"  Volume 4D: {cell_24['volume']}")
    print(f"  Raggio circoscritto: {cell_24['circumradius']}")
    print(f"  Raggio inscritto: {cell_24['inradius']}")
    
    # Quando proiettiamo da 4D a 3D, i volumi cambiscono
    # Il rapporto tra volumi 4D e 3D è legato a fattori geometrici
    
    print("\nProiezione:")
    print("  La proiezione 4D → 3D riduce le dimensioni")
    print("  Il rapporto tra volumi dipende dalla direzione di proiezione")
    
    # Possibile connessione: il volume di una sfera 4D vs 3D
    sphere_4d = pi**2 / 2
    sphere_3d = 4 * pi / 3
    
    print(f"\nVolume sfera 4D unitaria: {sphere_4d:.10f}")
    print(f"Volume sfera 3D unitaria: {sphere_3d:.10f}")
    print(f"Rapporto 4D/3D: {sphere_4d/sphere_3d:.10f}")
    print(f"Questo è ≈ 3π/8 = {3*pi/8:.10f}")


def main():
    print("=" * 70)
    print("INVESTIGAZIONE: SCALA ELETTROMAGNETICA E GEOMETRIA HURWITZ")
    print("=" * 70)
    
    # Calcola scala elettromagnetica
    em_scale = electromagnetic_scale()
    print(f"\nScala elettromagnetica e²/(ε₀c): {em_scale:.6e} J·s")
    
    # Testa scale geometriche
    geometric_scales = test_geometric_scales(em_scale)
    
    # Investiga ponte dimensionale
    investigate_dimensional_bridge()
    
    # Investiga scale Hurwitz
    investigate_hurwitz_scales()
    
    # Investiga proiezione 4D→3D
    investigate_4d_geometry_to_3d_projection()
    
    # Conclusione
    print("\n\n" + "=" * 70)
    print("CONCLUSIONI")
    print("=" * 70)
    
    print("\nRisultati:")
    print("  1. Le scale geometriche del reticolo di Hurwitz sono adimensionali")
    print("  2. La scala elettromagnetica ha dimensioni [J·s]")
    print("  3. Per collegarle, serve un 'fattore di scala' con dimensioni")
    print("\nOsservazioni:")
    print("  - e²/(ε₀c) = 2α × h (identità esatta)")
    print("  - 2α è il fattore di conversione tra h e la scala elettromagnetica")
    print("  - Il problema si riduce a: perché α = 1/137.036?")
    print("\nProssima pista:")
    print("  Investigare se α (o 1/α) può emergere dalla geometria")
    print("  del reticolo di Hurwitz o del 24-cell, come numero puro.")
    print("  Questo eviterebbe il problema dimensionale.")


if __name__ == "__main__":
    main()
