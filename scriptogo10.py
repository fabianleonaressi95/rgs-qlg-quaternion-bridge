#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=======================================================================
DEFINITIVE FIBONACCI / GOLDEN SPECTRAL AUDIT
=======================================================================

Obiettivo
---------
Determinare senza ambiguità l'origine di:

    C_F = 23.298975496882

e confrontarlo con:

    C_F = 18.849267128503...

Test eseguiti
-------------

A) c_a = c_b = 1
B) c_a = 1, c_b = 2

C) Fibonacci indexing:
       N = F_(n+1)
       N = F_(n+2)

D) spectral scalings:
       N^2 lambda_2
       phi^(2n) lambda_2
       phi^(2(n+1)) lambda_2

E) verifica analitica:

       C_F(1,2)
       = 5 pi^2 * 2/(1 + 2 phi)
       = 10 pi^2 (sqrt(5)-2)

F) confronto con il dataset originale:

       n = 2,...,12
       N = 3,5,8,13,...,377
       lambda_2 forniti dall'audit precedente

=======================================================================
"""

import numpy as np
from scipy.linalg import eigh


# =====================================================================
# COSTANTI
# =====================================================================

phi = (1.0 + np.sqrt(5.0)) / 2.0

sqrt5 = np.sqrt(5.0)

CF_TARGET = 23.298975496882

CF_UNIFORM = 5.0 * np.pi**2 / phi**2


# =====================================================================
# FIBONACCI
# =====================================================================

def fibonacci(k):
    """
    F_0 = 0
    F_1 = 1
    """

    if k == 0:
        return 0

    if k == 1:
        return 1

    a, b = 0, 1

    for _ in range(2, k + 1):
        a, b = b, a + b

    return b


# =====================================================================
# FIBONACCI WORD
# =====================================================================

def fibonacci_word(n):
    """
    Convenzione:

        w_0 = b
        w_1 = a
        w_(n+1) = w_n + w_(n-1)
    """

    if n == 0:
        return "b"

    if n == 1:
        return "a"

    wm1 = "b"
    w = "a"

    for _ in range(1, n):
        wp1 = w + wm1
        wm1 = w
        w = wp1

    return w


# =====================================================================
# LAPLACIANO PESATO
# =====================================================================

def weighted_laplacian(word, c_a, c_b):
    """
    Path Laplacian.

    Ogni lettera rappresenta un arco:

        a -> c_a
        b -> c_b

    N archi -> N+1 vertici.
    """

    N_edges = len(word)

    L = np.zeros(
        (N_edges + 1, N_edges + 1),
        dtype=float
    )

    for i, letter in enumerate(word):

        if letter == "a":
            c = c_a
        elif letter == "b":
            c = c_b
        else:
            raise ValueError("Simbolo Fibonacci non valido.")

        L[i, i] += c
        L[i + 1, i + 1] += c

        L[i, i + 1] -= c
        L[i + 1, i] -= c

    return L


# =====================================================================
# LAMBDA_2
# =====================================================================

def lambda2(word, c_a, c_b):

    L = weighted_laplacian(
        word,
        c_a,
        c_b
    )

    ev = eigh(
        L,
        eigvals_only=True
    )

    ev[np.abs(ev) < 1e-13] = 0.0

    return ev[1]


# =====================================================================
# C_F TEORICO
# =====================================================================

def CF_theory(c_a, c_b):

    return (
        5.0
        * np.pi**2
        * c_a
        * c_b
        /
        (c_a + phi*c_b)
    )


# =====================================================================
# D_eff
# =====================================================================

def D_eff(c_a, c_b):

    return 1.0 / (
        1.0/(phi*c_a)
        +
        1.0/(phi**2*c_b)
    )


# =====================================================================
# HEADER
# =====================================================================

print()
print("=" * 100)
print("DEFINITIVE FIBONACCI / GOLDEN SPECTRAL AUDIT")
print("=" * 100)

print()
print(f"phi   = {phi:.15f}")
print(f"sqrt5 = {sqrt5:.15f}")

print()
print(f"CF target   = {CF_TARGET:.15f}")
print(f"CF uniform  = {CF_UNIFORM:.15f}")


# =====================================================================
# 1. IDENTITÀ ANALITICA DI 23.298975...
# =====================================================================

print()
print("=" * 100)
print("1. EXACT ANALYTIC RECONSTRUCTION")
print("=" * 100)

CF_12 = CF_theory(1.0, 2.0)

CF_exact = (
    10.0
    * np.pi**2
    * (sqrt5 - 2.0)
)

CF_alt = (
    5.0
    * np.pi**2
    * 2.0
    /
    (1.0 + 2.0*phi)
)

print()
print("c_a = 1")
print("c_b = 2")

print()
print(f"C_F formula             = {CF_12:.15f}")
print(f"C_F = 10*pi^2*(sqrt5-2) = {CF_exact:.15f}")
print(f"C_F alternative         = {CF_alt:.15f}")

print()
print("absolute differences:")

print(
    "|formula - exact| =",
    abs(CF_12 - CF_exact)
)

print(
    "|formula - target| =",
    abs(CF_12 - CF_TARGET)
)


# =====================================================================
# 2. FATTORE AUREO ESATTO
# =====================================================================

print()
print("=" * 100)
print("2. GOLDEN FACTOR")
print("=" * 100)

CF_RATIO = CF_TARGET / CF_UNIFORM

print()
print(
    "CF_target / CF_uniform =",
    f"{CF_RATIO:.15f}"
)

print(
    "sqrt(5)-1             =",
    f"{sqrt5 - 1.0:.15f}"
)

print(
    "2/phi                 =",
    f"{2.0/phi:.15f}"
)

print()
print(
    "errors:"
)

print(
    "|ratio - (sqrt5-1)| =",
    abs(CF_RATIO - (sqrt5 - 1.0))
)

print(
    "|ratio - 2/phi|     =",
    abs(CF_RATIO - 2.0/phi)
)


# =====================================================================
# 3. D_eff
# =====================================================================

print()
print("=" * 100)
print("3. EFFECTIVE DIFFUSIVITY")
print("=" * 100)

for ca, cb in [
    (1.0, 1.0),
    (1.0, 2.0),
]:

    D = D_eff(ca, cb)
    CF = CF_theory(ca, cb)

    print()
    print(
        f"c_a={ca:.6f}, c_b={cb:.6f}"
    )

    print(
        f"D_eff = {D:.15f}"
    )

    print(
        f"C_F   = {CF:.15f}"
    )


# =====================================================================
# 4. GENERAZIONE SPETTRO PER c_a=c_b=1
# =====================================================================

print()
print("=" * 100)
print("4. SPECTRAL TEST: c_a = c_b = 1")
print("=" * 100)

print()
print(
    f"{'n':>3} "
    f"{'N_edges':>8} "
    f"{'lambda2':>18} "
    f"{'N² lambda2':>18} "
    f"{'phi²n lambda2':>20}"
)

print("-" * 75)

for n in range(2, 19):

    word = fibonacci_word(n)

    N = len(word)

    lam = lambda2(
        word,
        1.0,
        1.0
    )

    a = N**2 * lam

    b = phi**(2*n) * lam

    print(
        f"{n:3d} "
        f"{N:8d} "
        f"{lam:18.12e} "
        f"{a:18.12f} "
        f"{b:20.12f}"
    )


# =====================================================================
# 5. GENERAZIONE SPETTRO PER c_a=1,c_b=2
# =====================================================================

print()
print("=" * 100)
print("5. SPECTRAL TEST: c_a = 1, c_b = 2")
print("=" * 100)

print()
print(
    f"{'n':>3} "
    f"{'N_edges':>8} "
    f"{'lambda2':>18} "
    f"{'N² lambda2':>18} "
    f"{'phi²n lambda2':>20} "
    f"{'phi²(n+1) lambda2':>24}"
)

print("-" * 100)

values = []

for n in range(2, 19):

    word = fibonacci_word(n)

    N = len(word)

    lam = lambda2(
        word,
        1.0,
        2.0
    )

    a = N**2 * lam

    b = phi**(2*n) * lam

    c = phi**(2*(n+1)) * lam

    values.append(
        (
            n,
            N,
            lam,
            a,
            b,
            c
        )
    )

    print(
        f"{n:3d} "
        f"{N:8d} "
        f"{lam:18.12e} "
        f"{a:18.12f} "
        f"{b:20.12f} "
        f"{c:24.12f}"
    )


# =====================================================================
# 6. LIMIT ESTIMATE
# =====================================================================

print()
print("=" * 100)
print("6. ASYMPTOTIC LIMIT: c_a=1,c_b=2")
print("=" * 100)

tail = np.array(
    [x[4] for x in values[-6:]]
)

tail_shifted = np.array(
    [x[5] for x in values[-6:]]
)

CF_num = np.mean(tail)

CF_num_shifted = np.mean(tail_shifted)

print()
print("phi^(2n) lambda2 tail:")

for x in tail:
    print(f"{x:.15f}")

print()
print(
    f"mean = {CF_num:.15f}"
)

print()
print("phi^(2(n+1)) lambda2 tail:")

for x in tail_shifted:
    print(f"{x:.15f}")

print()
print(
    f"mean = {CF_num_shifted:.15f}"
)

print()
print(
    f"theory = {CF_12:.15f}"
)

print()
print(
    "error unshifted =",
    abs(CF_num - CF_12)
)

print(
    "error shifted   =",
    abs(CF_num_shifted - CF_12)
)


# =====================================================================
# 7. N = F_(n+1) vs F_(n+2)
# =====================================================================

print()
print("=" * 100)
print("7. FIBONACCI INDEXING AUDIT")
print("=" * 100)

print()
print(
    f"{'n':>3} "
    f"{'F(n+1)':>10} "
    f"{'F(n+2)':>10} "
    f"{'word length':>15}"
)

print("-" * 50)

for n in range(2, 13):

    f1 = fibonacci(n + 1)
    f2 = fibonacci(n + 2)

    Nword = len(fibonacci_word(n))

    print(
        f"{n:3d} "
        f"{f1:10d} "
        f"{f2:10d} "
        f"{Nword:15d}"
    )


# =====================================================================
# 8. RIPRODUZIONE ESATTA DEL VECCHIO DATASET
# =====================================================================

print()
print("=" * 100)
print("8. ORIGINAL DATASET REPRODUCTION")
print("=" * 100)

n_old = np.arange(2, 13)

N_old = np.array(
    [
        3,
        5,
        8,
        13,
        21,
        34,
        55,
        89,
        144,
        233,
        377
    ],
    dtype=float
)

lambda_old = np.array(
    [
        0.763932,
        0.318669,
        0.150018,
        0.061806,
        0.025138,
        0.009954,
        0.003888,
        0.001506,
        0.000580,
        0.000223,
        0.000085
    ]
)

A = phi**(2*n_old) * lambda_old

B = phi**(2*(n_old+1)) * lambda_old

E = N_old**2 * lambda_old

print()
print(
    f"{'n':>3} "
    f"{'N':>6} "
    f"{'lambda2':>14} "
    f"{'phi²n λ2':>16} "
    f"{'phi²(n+1) λ2':>20} "
    f"{'N² λ2':>14}"
)

print("-" * 90)

for i in range(len(n_old)):

    print(
        f"{n_old[i]:3d} "
        f"{int(N_old[i]):6d} "
        f"{lambda_old[i]:14.8f} "
        f"{A[i]:16.10f} "
        f"{B[i]:20.10f} "
        f"{E[i]:14.8f}"
    )


# =====================================================================
# 9. CONFRONTO DEL VECCHIO DATASET CON C_F
# =====================================================================

print()
print("=" * 100)
print("9. OLD DATASET vs TARGET")
print("=" * 100)

for name, arr in [
    ("phi^(2n) lambda2", A),
    ("phi^(2(n+1)) lambda2", B),
]:

    tail = arr[-5:]

    mean = np.mean(tail)

    print()
    print(name)

    print(
        f"tail mean = {mean:.15f}"
    )

    print(
        f"target    = {CF_TARGET:.15f}"
    )

    print(
        f"ratio     = {mean/CF_TARGET:.15f}"
    )

    print(
        f"error     = {mean-CF_TARGET:.15f}"
    )


# =====================================================================
# 10. TEST N vs N+1 SUL VECCHIO DATASET
# =====================================================================

print()
print("=" * 100)
print("10. OLD DATASET: N vs N+1")
print("=" * 100)

old_N_factor_plus = (
    A *
    ((N_old + 1.0) / N_old)**2
)

old_N_factor_minus = (
    A *
    (N_old / (N_old + 1.0))**2
)

for name, arr in [
    (
        "A",
        A
    ),
    (
        "A*(N+1)^2/N^2",
        old_N_factor_plus
    ),
    (
        "A*N^2/(N+1)^2",
        old_N_factor_minus
    ),
]:

    mean = np.mean(arr[-5:])

    print()
    print(
        f"{name:30s} = {mean:.15f}"
    )

    print(
        f"error to target             = "
        f"{mean-CF_TARGET:.15f}"
    )


# =====================================================================
# 11. RICERCA DELLO SHIFT DI INDICE
# =====================================================================

print()
print("=" * 100)
print("11. INDEX SHIFT SEARCH")
print("=" * 100)

for shift in range(-3, 4):

    scaled = (
        phi**(2*(n_old + shift))
        *
        lambda_old
    )

    mean = np.mean(
        scaled[-5:]
    )

    error = abs(
        mean - CF_TARGET
    )

    print(
        f"shift={shift:+d} "
        f"mean={mean:18.12f} "
        f"error={error:18.12f}"
    )


# =====================================================================
# 12. RICERCA DELLA NORMALIZZAZIONE phi^(2(n+s))
# =====================================================================

print()
print("=" * 100)
print("12. CONTINUOUS INDEX SHIFT")
print("=" * 100)

shifts = np.linspace(
    -1.0,
    2.0,
    3001
)

errors = []

for s in shifts:

    value = (
        phi**(2*(n_old+s))
        *
        lambda_old
    )

    mean = np.mean(
        value[-5:]
    )

    errors.append(
        abs(mean - CF_TARGET)
    )

errors = np.array(errors)

best = np.argmin(errors)

best_shift = shifts[best]

best_value = np.mean(
    phi**(2*(n_old+best_shift))
    *
    lambda_old
)[-5:] if False else np.mean(
    (
        phi**(2*(n_old+best_shift))
        *
        lambda_old
    )[-5:]
)

print()
print(
    f"best continuous shift = {best_shift:.12f}"
)

print(
    f"resulting mean         = {best_value:.12f}"
)

print(
    f"target                 = {CF_TARGET:.12f}"
)

print(
    f"absolute error          = "
    f"{abs(best_value-CF_TARGET):.12f}"
)


# =====================================================================
# 13. CONTROL UNIFORM PATH
# =====================================================================

print()
print("=" * 100)
print("13. EXACT UNIFORM PATH CONTROL")
print("=" * 100)

print()
print(
    f"{'N':>6} "
    f"{'numeric':>20} "
    f"{'exact':>20} "
    f"{'error':>15}"
)

print("-" * 70)

for N in [
    2,
    3,
    5,
    8,
    13,
    21,
    34
]:

    word = "a" * N

    lam = lambda2(
        word,
        1.0,
        1.0
    )

    exact = (
        2.0
        *
        (
            1.0
            -
            np.cos(
                np.pi/(N+1)
            )
        )
    )

    print(
        f"{N:6d} "
        f"{lam:20.15e} "
        f"{exact:20.15e} "
        f"{lam-exact:15.3e}"
    )


# =====================================================================
# 14. CONCLUSIONE AUTOMATICA
# =====================================================================

print()
print("=" * 100)
print("14. FINAL DIAGNOSTIC")
print("=" * 100)

print()

if abs(CF_12 - CF_TARGET) < 1e-10:

    print(
        "[PASS] 23.298975496882 è esattamente "
        "la formula C_F(c_a=1,c_b=2)."
    )

else:

    print(
        "[FAIL] La formula non ricostruisce il target."
    )


if abs(
    CF_RATIO - (sqrt5 - 1.0)
) < 1e-10:

    print(
        "[PASS] Il rapporto con il caso uniforme è "
        "sqrt(5)-1 = 2/phi."
    )

else:

    print(
        "[FAIL] Il rapporto non è il fattore aureo atteso."
    )


print()
print("RISULTATO ANALITICO:")
print()
print(
    "C_F(1,2)"
)
print(
    " = 5*pi^2 * 2/(1+2*phi)"
)
print(
    " = 10*pi^2*(sqrt(5)-2)"
)
print(
    f" = {CF_12:.15f}"
)

print()
print(
    "Quindi:"
)

print(
    f"23.298975496882 = {CF_12:.15f}"
)

print()
print("=" * 100)
print("AUDIT COMPLETED")
print("=" * 100)
