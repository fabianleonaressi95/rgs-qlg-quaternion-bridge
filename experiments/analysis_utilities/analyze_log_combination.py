#!/usr/bin/env python3
"""
Analisi approfondita della scoperta: log(5) + log(π) ≈ 11/4
"""

from math import log, pi, sqrt
from itertools import combinations


def analyze_log_combination():
    """Analizza la combinazione log(5) + log(π)."""
    print("=" * 70)
    print("ANALISI: log(5) + log(π) ≈ 11/4")
    print("=" * 70)
    
    # Valore calcolato
    log_5 = log(5)
    log_pi = log(pi)
    sum_log = log_5 + log_pi
    
    target = 11/4
    
    print(f"\nlog(5) = {log_5:.10f}")
    print(f"log(π) = {log_pi:.10f}")
    print(f"Somma = {sum_log:.10f}")
    print(f"Target 11/4 = {target:.10f}")
    print(f"Differenza = {sum_log - target:.10f}")
    print(f"Errore relativo = {abs(sum_log - target) / target * 100:.6f}%")
    
    # Interpretazione
    print(f"\nInterpretazione:")
    print(f"  log(5π) = log(5) + log(π) = {log(5*pi):.10f}")
    print(f"  5π = {5*pi:.10f}")
    print(f"  (5π)^(1/4) = {(5*pi)**(1/4):.10f}")
    
    # Questo produce 5^(1/4) × π^(1/4)
    # Ma Wyler richiede 5^(1/4) × π^(11/4)
    # Quindi manca π^(10/4) = π^(5/2)
    
    print(f"\nConfronto con Wyler:")
    print(f"  Trovato: (5π)^(1/4)")
    print(f"  Wyler: 5^(1/4) × π^(11/4)")
    print(f"  Gap: π^(10/4) = π^(5/2)")
    
    # π^(5/2) = π^2 × sqrt(π)
    print(f"\n  π^(5/2) = π^2 × √π = {pi**2 * sqrt(pi):.10f}")
    
    # Cerca se π^(5/2) emerge da qualche parte
    print(f"\nRicerca π^(5/2):")
    print(f"  Volume S³ = 2π²")
    print(f"  (Volume S³)^(1/2) = √(2π²) = {sqrt(2*pi**2):.10f}")
    print(f"  Non è π^(5/2)")
    
    # Altre combinazioni per 19/4
    print(f"\n{'='*70}")
    print("RICERCA COMBINAZIONI PER 19/4")
    print(f"{'='*70}")
    
    target_19_4 = 19/4
    
    # Valori disponibili
    vals = {
        'log(2)': log(2),
        'log(3)': log(3),
        'log(5)': log(5),
        'log(π)': log(pi),
        'log(5!)': log(120),
        'log(24)': log(24),
        'log(48)': log(48),
        'log(64)': log(64),
        'log(216)': log(216),
        'log(Vol S³)': log(2*pi**2),
    }
    
    # Cerca combinazioni di 3 termini
    print(f"\nCombinazioni di 3 termini:")
    from itertools import combinations
    
    for combo in combinations(vals.keys(), 3):
        sum_val = sum(vals[name] for name in combo)
        if abs(sum_val - target_19_4) < 0.1:
            print(f"  {' + '.join(combo)} = {sum_val:.4f} (target: {target_19_4:.4f})")
    
    # Cerca combinazioni con coefficienti interi
    print(f"\nCombinazioni lineari con coefficienti interi:")
    
    # 19/4 = 4.75
    # Possibili: 2×log(2) + log(π) = 2×0.693 + 1.145 = 2.531
    # 3×log(2) + log(π) = 3×0.693 + 1.145 = 3.224
    # log(24) = 3.178
    # log(24) + log(π) = 3.178 + 1.145 = 4.323
    # log(48) = 3.871
    # log(48) + log(π) = 3.871 + 1.145 = 5.016
    
    combos_int = [
        ('log(24) + log(π)', log(24) + log(pi)),
        ('log(48) + log(π)', log(48) + log(pi)),
        ('log(5!) - log(2)', log(120) - log(2)),
        ('log(5!) - log(π)', log(120) - log(pi)),
        ('log(5!) - log(3)', log(120) - log(3)),
    ]
    
    for desc, val in combos_int:
        diff = abs(val - target_19_4)
        if diff < 1.0:
            print(f"  {desc} = {val:.4f} (diff: {diff:.4f})")
    
    # Analisi specifica per -7/4
    print(f"\n{'='*70}")
    print("RICERCA COMBINAZIONI PER -7/4")
    print(f"{'='*70}")
    
    target_minus_7_4 = -7/4
    
    # Possibili combinazioni negative
    combos_neg = [
        ('-log(2)', -log(2)),
        ('-log(3)', -log(3)),
        ('-log(5)', -log(5)),
        ('-log(π)', -log(pi)),
        ('log(1/2)', log(0.5)),
        ('log(1/3)', log(1/3)),
        ('log(1/5)', log(0.2)),
    ]
    
    for desc, val in combos_neg:
        diff = abs(val - target_minus_7_4)
        if diff < 1.0:
            print(f"  {desc} = {val:.4f} (target: {target_minus_7_4:.4f}, diff: {diff:.4f})")
    
    # log(1/3) = -1.099 ≈ -1.75 = -7/4
    # Gap: 0.65
    
    print(f"\n{'='*70}")
    print("CONCLUSIONE")
    print(f"{'='*70}")
    
    print(f"\n✓ Trovato: log(5) + log(π) ≈ 11/4")
    print(f"✗ Nessuna combinazione semplice per 19/4")
    print(f"✗ Nessuna combinazione semplice per -7/4")
    
    print(f"\nL'unica corrispondenza parziale è:")
    print(f"  11/4 ↔ log(5π)")
    print(f"  Ma Wyler richiede π^(11/4), non (5π)^(1/4)")


if __name__ == "__main__":
    analyze_log_combination()
