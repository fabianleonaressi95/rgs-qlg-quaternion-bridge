#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
========================================================================
FIBONACCI SPECTRAL PREFactor AUDIT
========================================================================

Obiettivo:
    spiegare la discrepanza tra

        C_F,theory = 23.298975496882

e il limite numerico osservato

        phi^(2n) * lambda_2 ~ 8.9

Il programma verifica:

1. convenzione Fibonacci / lunghezza N_n
2. costruzione esatta della parola Fibonacci
3. costruzione del Laplaciano pesato
4. autovalori esatti
5. distinzione tra:
       - lambda_2
       - lambda_1 non nullo
       - normalized Laplacian
       - generalized eigenvalue
6. confronto N_n^2 lambda_2
7. confronto phi^(2n) lambda_2
8. formula teorica D_eff
9. formula teorica C_F
10. possibili normalizzazioni alternative
11. fit asintotico numerico
12. ricerca automatica del fattore di discrepanza
13. test del caso uniforme c_a = c_b
========================================================================
"""

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

# ============================================================================
# PARAMETRI
# ============================================================================

PHI = (1.0 + np.sqrt(5.0)) / 2.0

# Valore teorico riportato nel tuo audit
CF_REPORTED = 23.298975496882

# Se conosci già c_a e c_b, inseriscili qui.
# Altrimenti lo script può tentare una diagnosi dalla struttura numerica.
C_A = 1.0
C_B = 1.0

N_MIN = 2
N_MAX = 18


# ============================================================================
# 1. FIBONACCI SUBSTITUTION
# ============================================================================

def fibonacci_word(n):
    """
    Convenzione del documento:

        w_0 = b
        w_1 = a
        w_{n+1} = w_n w_{n-1}

    quindi:

        N_n = F_{n+1}
    """

    w0 = "b"
    w1 = "a"

    if n == 0:
        return w0
    if n == 1:
        return w1

    wm1 = w0
    w = w1

    for _ in range(1, n):
        wp1 = w + wm1
        wm1, w = w, wp1

    return w


# ============================================================================
# 2. LUNGHEZZE
# ============================================================================

def fibonacci_lengths(nmax):
    return [len(fibonacci_word(n)) for n in range(nmax + 1)]


# ============================================================================
# 3. LAPLACIANO PESATO
# ============================================================================

def weighted_path_laplacian(word, c_a, c_b):
    """
    Path con N archi e N+1 vertici.

    arco i:
        c_a se lettera = a
        c_b se lettera = b
    """

    N = len(word)
    L = np.zeros((N + 1, N + 1), dtype=float)

    for i, s in enumerate(word):
        c = c_a if s == "a" else c_b

        # vertici i e i+1
        L[i, i] += c
        L[i + 1, i + 1] += c

        L[i, i + 1] -= c
        L[i + 1, i] -= c

    return L


# ============================================================================
# 4. SPETTRO ESATTO
# ============================================================================

def spectrum(L):
    """
    Laplaciano simmetrico -> eigh.
    """

    vals = eigh(L, eigvals_only=True)

    # pulizia numerica del modo zero
    vals[np.abs(vals) < 1e-12] = 0.0

    return vals


# ============================================================================
# 5. EFFETTIVE DIFFUSIVITY
# ============================================================================

def D_eff(c_a, c_b):
    """
    Formula del documento:

    D_eff =
        [1/(phi*c_a) + 1/(phi^2*c_b)]^{-1}
    """

    return 1.0 / (
        1.0 / (PHI * c_a)
        +
        1.0 / (PHI**2 * c_b)
    )


def C_F_theory(c_a, c_b):
    """
    Formula dichiarata:

        C_F = 5*pi^2*c_a*c_b / (c_a + phi*c_b)
    """

    return (
        5.0
        * np.pi**2
        * c_a
        * c_b
        / (c_a + PHI * c_b)
    )


# ============================================================================
# 6. ALTERNATIVE NORMALIZZAZIONI
# ============================================================================

def C_F_from_D_eff(c_a, c_b):

    D = D_eff(c_a, c_b)

    return (
        5.0
        * np.pi**2
        / PHI**2
        * D
    )


def C_F_uniform(c):

    return 5.0 * np.pi**2 * c / PHI**2


# ============================================================================
# 7. ANALISI PRINCIPALE
# ============================================================================

rows = []

print("\n")
print("=" * 100)
print("FIBONACCI SPECTRAL PREFactor AUDIT")
print("=" * 100)

print(f"\nphi = {PHI:.15f}")
print(f"c_a = {C_A}")
print(f"c_b = {C_B}")

D = D_eff(C_A, C_B)
CF = C_F_theory(C_A, C_B)

print(f"\nD_eff               = {D:.15f}")
print(f"C_F theory formula = {CF:.15f}")
print(f"C_F reported       = {CF_REPORTED:.15f}")

print("\n")


for n in range(N_MIN, N_MAX + 1):

    word = fibonacci_word(n)
    N = len(word)

    L = weighted_path_laplacian(word, C_A, C_B)
    ev = spectrum(L)

    # ev[0] = zero mode
    lambda2 = ev[1]

    scaled_N = N**2 * lambda2
    scaled_phi = PHI**(2*n) * lambda2

    # teorico al livello N
    lambda_asym_N = np.pi**2 * D / N**2

    # rapporto rispetto all'asintotico N^-2
    ratio_N = scaled_N / (np.pi**2 * D)

    # rapporto rispetto a CF
    ratio_CF = scaled_phi / CF

    rows.append(
        (
            n,
            N,
            lambda2,
            scaled_N,
            scaled_phi,
            ratio_N,
            ratio_CF
        )
    )


# ============================================================================
# 8. STAMPA
# ============================================================================

header = (
    f"{'n':>3} "
    f"{'N':>7} "
    f"{'lambda2':>16} "
    f"{'N^2 lambda2':>16} "
    f"{'phi^(2n) lambda2':>20} "
    f"{'ratio N-asym':>15} "
    f"{'ratio CF':>12}"
)

print(header)
print("-" * len(header))

for r in rows:

    print(
        f"{r[0]:3d} "
        f"{r[1]:7d} "
        f"{r[2]:16.10e} "
        f"{r[3]:16.10f} "
        f"{r[4]:20.10f} "
        f"{r[5]:15.8f} "
        f"{r[6]:12.8f}"
    )


# ============================================================================
# 9. LIMITE NUMERICO
# ============================================================================

tail = np.array([r[4] for r in rows[-6:]])

CF_NUM = np.mean(tail)

print("\n")
print("=" * 100)
print("LIMIT ESTIMATE")
print("=" * 100)

print(f"\nCF numerical tail mean = {CF_NUM:.12f}")

print(f"CF theory             = {CF:.12f}")

print(f"CF reported           = {CF_REPORTED:.12f}")

print(
    f"\nNumerical / reported = "
    f"{CF_NUM / CF_REPORTED:.12f}"
)

print(
    f"Reported / numerical = "
    f"{CF_REPORTED / CF_NUM:.12f}"
)


# ============================================================================
# 10. D_eff RICAVATO DIRETTAMENTE DAI DATI
# ============================================================================

# Dal teorema:
#
# N^2 lambda2 -> pi^2 D_eff

scaled_N_values = np.array([r[3] for r in rows])

D_num = np.mean(
    scaled_N_values[-6:] / np.pi**2
)

print("\n")
print("=" * 100)
print("EFFECTIVE DIFFUSIVITY AUDIT")
print("=" * 100)

print(f"\nD_eff theoretical = {D:.12f}")
print(f"D_eff numerical   = {D_num:.12f}")

print(
    f"ratio numerical/theory = "
    f"{D_num / D:.12f}"
)


# ============================================================================
# 11. CONTROLLO DEL FATTORE 5/phi^2
# ============================================================================

conversion_factor = 5.0 / PHI**2

CF_from_numeric_D = conversion_factor * D_num

print("\n")
print("=" * 100)
print("PHI NORMALIZATION AUDIT")
print("=" * 100)

print(f"\n5/phi^2 = {conversion_factor:.15f}")

print(
    f"\n(5/phi^2) * D_eff numerical = "
    f"{CF_from_numeric_D:.12f}"
)

print(
    f"pi^2 * D_eff numerical = "
    f"{np.pi**2 * D_num:.12f}"
)

print(
    f"\nObserved CF numerical = "
    f"{CF_NUM:.12f}"
)


# ============================================================================
# 12. CHECK DIRETTO DEL RAPPORTO phi^(2n) / N^2
# ============================================================================

print("\n")
print("=" * 100)
print("BINET NORMALIZATION CHECK")
print("=" * 100)

for n in range(N_MIN, N_MAX + 1):

    N = len(fibonacci_word(n))

    empirical = PHI**(2*n) / N**2
    theoretical = 5.0 / PHI**2

    print(
        f"n={n:2d} "
        f"N={N:5d} "
        f"phi^(2n)/N^2={empirical:.12f} "
        f"target={theoretical:.12f} "
        f"ratio={empirical/theoretical:.12f}"
    )


# ============================================================================
# 13. CASO UNIFORME
# ============================================================================

print("\n")
print("=" * 100)
print("UNIFORM CASE TEST")
print("=" * 100)

c = 1.0

CF_uniform_theory = C_F_uniform(c)

print(f"\nc_a = c_b = {c}")

print(
    f"C_F uniform theory = "
    f"{CF_uniform_theory:.12f}"
)

uniform_values = []

for n in range(4, 18):

    word = fibonacci_word(n)
    N = len(word)

    L = weighted_path_laplacian(word, c, c)

    ev = spectrum(L)

    lam = ev[1]

    val = PHI**(2*n) * lam

    uniform_values.append(val)

    print(
        f"n={n:2d} "
        f"N={N:5d} "
        f"phi^(2n)*lambda2={val:.12f}"
    )

uniform_limit = np.mean(uniform_values[-5:])

print(
    f"\nUniform numerical limit = "
    f"{uniform_limit:.12f}"
)

print(
    f"Uniform theoretical limit = "
    f"{CF_uniform_theory:.12f}"
)


# ============================================================================
# 14. CONTROL: PATH UNIFORME ESATTO
# ============================================================================

print("\n")
print("=" * 100)
print("EXACT UNIFORM PATH FORMULA CHECK")
print("=" * 100)

for N in [3, 5, 8, 13, 21, 34]:

    # N archi -> N+1 vertici
    L = weighted_path_laplacian("a" * N, 1.0, 1.0)

    ev = spectrum(L)

    numerical = ev[1]

    exact = 2.0 * (
        1.0
        -
        np.cos(np.pi / (N + 1))
    )

    print(
        f"N={N:4d} "
        f"numeric={numerical:.14e} "
        f"exact={exact:.14e} "
        f"error={numerical-exact:.3e}"
    )


# ============================================================================
# 15. IDENTIFICAZIONE AUTOMATICA DEL FATTORE MANCANTE
# ============================================================================

print("\n")
print("=" * 100)
print("DISCREPANCY FACTOR DECOMPOSITION")
print("=" * 100)

ratio = CF_REPORTED / CF_NUM

print(f"\nCF_reported / CF_numeric = {ratio:.12f}")

candidates = {
    "phi": PHI,
    "phi^2": PHI**2,
    "sqrt(5)": np.sqrt(5),
    "5/phi": 5.0 / PHI,
    "phi^2/sqrt(5)": PHI**2 / np.sqrt(5),
    "pi": np.pi,
    "pi/2": np.pi / 2,
    "2": 2.0,
    "N+1 vs N asymptotic": 1.0,
}

print("\nPossible simple factors:")

for name, value in candidates.items():

    print(
        f"{name:25s} "
        f"{value:15.10f} "
        f"ratio/factor={ratio/value:15.10f}"
    )


# ============================================================================
# 16. STIMA ASINTOTICA
# ============================================================================

print("\n")
print("=" * 100)
print("ASYMPTOTIC FIT")
print("=" * 100)

#
# Modello:
#
# phi^(2n) lambda2
# =
# CF + A*(-phi^-2)^n + B*phi^(-2n)
#
# oppure semplicemente
#
# CF + A*r^n
#

ns = np.array([r[0] for r in rows], dtype=float)
ys = np.array([r[4] for r in rows], dtype=float)


def model(p):

    CF0, A, r = p

    return CF0 + A * r**ns


def residual(p):

    return model(p) - ys


fit = least_squares(
    residual,
    x0=[
        ys[-1],
        ys[-1] - ys[-2],
        -PHI**-2
    ],
    bounds=(
        [-np.inf, -np.inf, -0.999],
        [np.inf, np.inf, 0.999]
    )
)

CF_fit, A_fit, r_fit = fit.x

print(f"\nCF fit = {CF_fit:.12f}")
print(f"A      = {A_fit:.12f}")
print(f"r      = {r_fit:.12f}")

print(
    f"\nExpected subleading scale "
    f"-phi^-2 = {-PHI**-2:.12f}"
)


# ============================================================================
# 17. TEST NORMALIZED LAPLACIAN
# ============================================================================

print("\n")
print("=" * 100)
print("NORMALIZED LAPLACIAN CONTROL")
print("=" * 100)


def normalized_laplacian(L):

    d = np.diag(L)

    Dm12 = np.zeros_like(L)

    mask = d > 1e-14

    Dm12[mask, mask] = 1.0 / np.sqrt(d[mask])

    return (
        np.eye(len(L))
        -
        Dm12 @ L @ Dm12
    )


for n in range(5, 11):

    word = fibonacci_word(n)

    L = weighted_path_laplacian(word, C_A, C_B)

    Ln = normalized_laplacian(L)

    ev = spectrum(Ln)

    print(
        f"n={n:2d} "
        f"lambda2(normalized)={ev[1]:.12e}"
    )


# ============================================================================
# 18. CHECK COUNTEGGIO ARCHI vs VERTICI
# ============================================================================

print("\n")
print("=" * 100)
print("N vs N+1 NORMALIZATION CHECK")
print("=" * 100)

for n in range(5, 13):

    word = fibonacci_word(n)

    N_edges = len(word)
    N_vertices = N_edges + 1

    L = weighted_path_laplacian(word, C_A, C_B)

    ev = spectrum(L)

    lam = ev[1]

    c_edges = PHI**(2*n) * lam

    # normalizzazione artificiale usando N+1
    c_vertices = (
        PHI**(2*n)
        * lam
        * (N_vertices / N_edges)**2
    )

    print(
        f"n={n:2d} "
        f"N_edges={N_edges:5d} "
        f"N_vertices={N_vertices:5d} "
        f"CF_edges={c_edges:12.8f} "
        f"CF_vertices={c_vertices:12.8f}"
    )


# ============================================================================
# 19. REPORT FINALE AUTOMATICO
# ============================================================================

print("\n")
print("=" * 100)
print("FINAL DIAGNOSTIC")
print("=" * 100)

print(f"""
1. Numerical CF:
       {CF_NUM:.12f}

2. Reported CF:
       {CF_REPORTED:.12f}

3. Theory from current c_a,c_b:
       {CF:.12f}

4. Numerical D_eff:
       {D_num:.12f}

5. Theoretical D_eff:
       {D:.12f}

6. Discrepancy:
       reported / numerical =
       {ratio:.12f}

7. Binet factor:
       5/phi^2 =
       {5/PHI**2:.12f}

8. Fitted CF:
       {CF_fit:.12f}

9. Fitted subleading eigenvalue:
       {r_fit:.12f}

10. Expected renormalization eigenvalue:
       -phi^-2 =
       {-PHI**-2:.12f}
""")


# ============================================================================
# 20. CONCLUSIONE
# ============================================================================

print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

if abs(CF_NUM - CF) / abs(CF) < 0.05:

    print("""
[PASS]
Il limite numerico è compatibile con la formula teorica
per i valori correnti di c_a e c_b.
""")

else:

    print("""
[FAIL]
Il limite numerico NON è compatibile con la formula teorica
per i valori correnti di c_a e c_b.

Occorre controllare:

  A) c_a, c_b
  B) definizione di lambda_2
  C) costruzione del Laplaciano
  D) convenzione N_n
  E) normalizzazione spettro
  F) fattore di scala phi^(2n)
  G) formula D_eff
  H) formula del prefattore C_F
""")

print("=" * 100)
print("AUDIT COMPLETED")
print("=" * 100)
