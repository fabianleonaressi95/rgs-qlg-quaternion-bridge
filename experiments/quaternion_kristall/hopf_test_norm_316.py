#!/usr/bin/env python3
"""
Test Hopf corretto usando i 64 quaternioni autentici della norma 316.
Implementa h(q) = q*i*conj(q) / |q|²
"""

import numpy as np
from math import sqrt, acos, pi
import csv


def load_quaternions(csv_file):
    """Carica quaternioni dal file CSV."""
    quaternions = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row['t'])
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            group = row['group']
            quaternions.append((t, x, y, z, group))
    return quaternions


def quaternion_conjugate(q):
    """Coniugato di un quaternione."""
    t, x, y, z = q
    return (t, -x, -y, -z)


def quaternion_multiply(q1, q2):
    """Moltiplicazione quaternioni."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2
    )


def quaternion_norm(q):
    """Norma di un quaternione."""
    t, x, y, z = q
    return sqrt(t*t + x*x + y*y + z*z)


def normalize_quaternion(q):
    """Normalizza un quaternione a norma 1."""
    norm = quaternion_norm(q)
    t, x, y, z = q
    return (t/norm, x/norm, y/norm, z/norm)


def hopf_map(q):
    """
    Mappa di Hopf: h(q) = q*i*conj(q) / |q|²
    Proietta S³ → S²
    """
    # Quaternione i = (0, 1, 0, 0)
    i_quat = (0, 1, 0, 0)
    
    # Normalizza q
    q_norm = normalize_quaternion(q)
    
    # Calcola q * i * conj(q)
    conj_q = quaternion_conjugate(q_norm)
    temp = quaternion_multiply(q_norm, i_quat)
    result = quaternion_multiply(temp, conj_q)
    
    # Il risultato è un quaternione puro (parte reale = 0)
    # La parte immaginaria è il punto su S²
    return result[1:]  # (x, y, z)


def angular_distance(p1, p2):
    """Distanza angolare tra due punti su S²."""
    # Prodotto scalare
    dot = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
    # Clamp a [-1, 1] per evitare errori numerici
    dot = max(-1, min(1, dot))
    return acos(dot)


def build_graph(points, threshold=0.1):
    """
    Costruisce il grafo collegando punti con distanza angolare < threshold.
    """
    n = len(points)
    edges = []
    
    for i in range(n):
        for j in range(i+1, n):
            dist = angular_distance(points[i], points[j])
            if dist < threshold and dist > 1e-10:  # Evita auto-loop
                edges.append((i, j, dist))
    
    return edges


def compute_laplacian(n, edges):
    """Calcola il Laplaciano del grafo."""
    L = np.zeros((n, n))
    
    for i, j, _ in edges:
        L[i, j] = -1
        L[j, i] = -1
    
    # Diagonale = grado
    for i in range(n):
        L[i, i] = -sum(L[i, :])
    
    return L


def main():
    print("=" * 70)
    print("TEST HOPF CORRETTO: 64 QUATERNIONI NORMA 316")
    print("=" * 70)
    
    # Carica quaternioni
    quaternions = load_quaternions('norm_316_quaternions.csv')
    print(f"\nQuaternioni caricati: {len(quaternions)}")
    
    # Verifica distintezza
    unique_quats = set((round(q[0], 6), round(q[1], 6), round(q[2], 6), round(q[3], 6)) for q in quaternions)
    print(f"Quaternioni distinti: {len(unique_quats)}")
    
    # Applica mappa di Hopf
    print(f"\nApplicazione mappa di Hopf h(q) = q*i*conj(q)/|q|²...")
    hopf_images = []
    
    for i, q in enumerate(quaternions):
        image = hopf_map(q[:4])
        norm = sqrt(image[0]**2 + image[1]**2 + image[2]**2)
        hopf_images.append(image)
        
        if i < 5:
            print(f"  q_{i}: {q[:4]} → h(q) = {image} (norma = {norm:.10f})")
    
    # Verifica che le immagini siano su S²
    print(f"\nVerifica che le immagini siano su S²:")
    norms = [sqrt(p[0]**2 + p[1]**2 + p[2]**2) for p in hopf_images]
    print(f"  Norma minima: {min(norms):.10f}")
    print(f"  Norma massima: {max(norms):.10f}")
    
    if all(abs(n - 1) < 1e-6 for n in norms):
        print("  ✓ Tutte le immagini hanno norma 1 (sono su S²)")
    else:
        print("  ✗ Alcune immagini non hanno norma 1")
    
    # Verifica distintezza delle immagini
    unique_images = set((round(p[0], 6), round(p[1], 6), round(p[2], 6)) for p in hopf_images)
    print(f"\nImagini Hopf distinte: {len(unique_images)}")
    
    if len(unique_images) == len(hopf_images):
        print("  ✓ Tutte le immagini sono distinte")
    else:
        print(f"  ✗ Collasso: {len(hopf_images)} → {len(unique_images)}")
        print(f"  Rapporto di collasso: {len(unique_images)/len(hopf_images):.2f}")
    
    # Costruisci grafo
    print(f"\nCostruzione grafo delle immagini...")
    threshold = 0.5  # Soglia angolare in radianti (~28.6 gradi)
    edges = build_graph(hopf_images, threshold=threshold)
    
    print(f"  Soglia angolare: {threshold:.4f} rad ({threshold*180/pi:.2f}°)")
    print(f"  Numero archi: {len(edges)}")
    
    # Calcola Laplaciano
    n = len(hopf_images)
    L = compute_laplacian(n, edges)
    
    # Calcola autovalori
    eigenvalues = np.linalg.eigvalsh(L)
    eigenvalues = np.round(eigenvalues, 6)
    
    print(f"\nAutovalori del Laplaciano:")
    unique_eigenvals, counts = np.unique(eigenvalues, return_counts=True)
    for val, count in zip(unique_eigenvals, counts):
        print(f"  λ = {val:.6f} (molteplicità {count})")
    
    # Analisi risultati
    print(f"\n{'='*70}")
    print("ANALISI RISULTATI")
    print(f"{'='*70}")
    
    print(f"\nCollasso Hopf:")
    print(f"  64 quaternioni → {len(unique_images)} immagini distinte")
    
    print(f"\nGrafo:")
    print(f"  {n} vertici, {len(edges)} archi")
    
    print(f"\nAutovalori significativi:")
    nonzero_eigenvals = [val for val in eigenvalues if abs(val) > 1e-6]
    print(f"  Autovalori non nulli: {len(nonzero_eigenvals)}")
    
    if nonzero_eigenvals:
        print(f"  Autovalori unici: {np.unique(nonzero_eigenvals)}")
    
    # Cerca sequenza armonica 1, 3, 5, 7, ...
    print(f"\nRicerca sequenza armonica l(l+1):")
    expected = [2, 6, 12, 20, 30, 42, 56]  # l(l+1) per l=1,2,3,4,5,6,7
    found = []
    for val in expected:
        if any(abs(ev - val) < 0.1 for ev in nonzero_eigenvals):
            found.append(val)
    
    if found:
        print(f"  ✓ Trovati: {found}")
    else:
        print(f"  ✗ Nessun autovalore corrisponde a l(l+1)")
    
    # Analisi per gruppo
    print(f"\n{'='*70}")
    print("ANALISI PER GRUPPO (216 vs 300)")
    print(f"{'='*70}")
    
    group_216 = [i for i, q in enumerate(quaternions) if q[4] == '216']
    group_300 = [i for i, q in enumerate(quaternions) if q[4] == '300']
    
    print(f"\nGruppo 216: {len(group_216)} quaternioni")
    print(f"Gruppo 300: {len(group_300)} quaternioni")
    
    # Calcola immagini per gruppo
    images_216 = [hopf_images[i] for i in group_216]
    images_300 = [hopf_images[i] for i in group_300]
    
    unique_216 = set((round(p[0], 6), round(p[1], 6), round(p[2], 6)) for p in images_216)
    unique_300 = set((round(p[0], 6), round(p[1], 6), round(p[2], 6)) for p in images_300)
    
    print(f"\nImagini distinte:")
    print(f"  Gruppo 216: {len(unique_216)}")
    print(f"  Gruppo 300: {len(unique_300)}")
    
    print(f"\nConclusione:")
    if len(unique_images) == 64:
        print("  ✓ La mappa di Hopf è iniettiva su questo campione")
        print("  ✓ Possiamo testare la struttura armonica")
    else:
        print(f"  ✗ Collasso a {len(unique_images)} immagini")
        print("  La mappa non è iniettiva - la fibra S¹ identifica punti")


if __name__ == "__main__":
    main()
