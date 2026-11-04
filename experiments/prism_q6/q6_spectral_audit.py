#!/usr/bin/env python3
"""
Audit reale dello spettro del Laplaciano combinatorio di Q6 (ipercubo 6D).
Costruisce la matrice, calcola autovalori e confronta con la formula esatta
λ_k = 2k con molteplicità C(6,k).
"""
import numpy as np
from itertools import product
from math import comb

print("=" * 60)
print("Q6 SPECTRAL AUDIT – calcolo reale")
print("=" * 60)

# Vertici di Q6: {0,1}^6
verts = list(product([0, 1], repeat=6))
n = len(verts)  # 64
assert n == 64

# Mappa vertice -> indice
idx = {v: i for i, v in enumerate(verts)}

# Laplaciano combinatorio: L = D - A
# Grado costante = 6
L = np.zeros((n, n), dtype=float)
for v in verts:
    i = idx[v]
    L[i, i] = 6.0
    for d in range(6):
        w = list(v)
        w[d] = 1 - w[d]
        j = idx[tuple(w)]
        L[i, j] = -1.0

# Autovalori numerici
evals = np.sort(np.linalg.eigvalsh(L))
print(f"Dimensione matrice : {n}×{n}")
print(f"Autovalori numerici (primi 10): {evals[:10].round(6)}")

# Formula esatta: λ = 2k, molteplicità C(6,k), k=0..6
exact = []
for k in range(7):
    exact.extend([2.0 * k] * comb(6, k))
exact = np.sort(exact)

residue = np.max(np.abs(evals - exact))
print(f"\nResiduo max |λ_num - λ_exact| = {residue:.3e}")

# Conteggio molteplicità
from collections import Counter
c_num = Counter(np.round(evals, 8))
c_ex  = Counter(exact)
print("\nMolteplicità:")
print(f"  {'k':>3}  {'λ=2k':>6}  {'exact':>6}  {'numerico':>8}")
for k in range(7):
    print(f"  {k:>3}  {2*k:>6}  {comb(6,k):>6}  {c_num[2.0*k]:>8}")

if residue < 1e-10:
    print("\nSPECTRUM_CERTIFIED (residuo numerico < 1e-10)")
else:
    print("\nSPECTRUM_FAILED")
print("=" * 60)
