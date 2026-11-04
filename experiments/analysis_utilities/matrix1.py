#!/usr/bin/env python3

# ============================================================
# QUASICRYSTAL / QNN GROUP-GEOMETRY AUDIT
#
# 64 STATES = F_2^6
#
# Tests:
#   1. XOR equivariance
#   2. S6 invariance
#   3. Hamming-shell structure
#   4. D(G,G') = f(d_H(G,G'))
#   5. null random-64
#
# IMPORTANT:
# The synthetic models are controls.
# They do NOT reconstruct an experimental QNN unless
# the actual P(a|G) distributions are supplied.
# ============================================================

import os
import itertools
import numpy as np
import pandas as pd

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr


# ============================================================
# CONFIGURATION
# ============================================================

SEED = 12345

N_STATES = 64
N_BITS = 6

OUTDIR = "qnn_hamming_audit"

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# 1. F_2^6
# ============================================================

states = np.arange(N_STATES, dtype=np.uint8)

bits = np.array(
    [
        [(g >> i) & 1 for i in range(N_BITS)]
        for g in states
    ],
    dtype=np.uint8
)


def hamming_state(a, b):
    return int(np.count_nonzero(bits[a] != bits[b]))


# ============================================================
# 2. HAMMING MATRIX
# ============================================================

D_H = np.zeros(
    (N_STATES, N_STATES),
    dtype=np.int8
)

for i in states:
    for j in states:
        D_H[i, j] = hamming_state(i, j)


pd.DataFrame(D_H).to_csv(
    os.path.join(OUTDIR, "D_Hamming_64x64.csv"),
    index=False
)


# ============================================================
# 3. S6
# ============================================================

S6 = list(itertools.permutations(range(N_BITS)))

print()
print("=" * 78)
print("S6 GROUP")
print("=" * 78)
print("Number of permutations:", len(S6))


def permute_state(g, permutation):
    """
    Natural S6 action on the six bits.
    """

    b = bits[g]

    bp = b[list(permutation)]

    out = 0

    for i, value in enumerate(bp):
        out |= int(value) << i

    return out


# Verify that S6 preserves Hamming weight
s6_weight_errors = 0

for g in states:

    w = int(bits[g].sum())

    for p in S6:

        gp = permute_state(g, p)

        if int(bits[gp].sum()) != w:
            s6_weight_errors += 1

print("S6 weight-preservation errors:", s6_weight_errors)


# ============================================================
# 4. SOFTMAX
# ============================================================

def softmax(x):

    x = np.asarray(x, dtype=float)

    x = x - np.max(x)

    e = np.exp(x)

    return e / np.sum(e)


# ============================================================
# 5. EXACT D2
#
# P(a|G) = K(a XOR G)
#
# K(x) = K(|x|)
#
# Therefore:
#
# P(a|G XOR v) = P(a XOR v | G)
#
# and
#
# K(delta) = K(pi delta)
#
# for every pi in S6.
#
# Therefore:
#
# D(G,G') = f(d_H(G,G'))
#
# exactly, up to floating point precision.
# ============================================================

def make_D2_exact(
    seed=SEED,
    beta=1.0,
    shell_scale=0.25
):

    rng = np.random.default_rng(seed)

    # ONLY seven independent parameters:
    #
    # q[0], ..., q[6]
    #
    # one for each Hamming shell
    q_shell = rng.normal(
        0.0,
        shell_scale,
        7
    )

    P = np.zeros(
        (N_STATES, N_STATES),
        dtype=float
    )

    for G in states:

        logits = np.zeros(N_STATES)

        for a in states:

            d = hamming_state(a, G)

            logits[a] = beta * q_shell[d]

        P[G] = softmax(logits)

    return P


# ============================================================
# 6. PERTURBED D2
#
# Start from the exact D2 kernel and add a controlled
# symmetry-breaking perturbation.
#
# epsilon = 0
#       exact D2
#
# epsilon > 0
#       increasingly breaks S6.
# ============================================================

def make_D2_perturbed(
    seed=SEED,
    beta=1.0,
    shell_scale=0.25,
    epsilon=1e-2
):

    rng = np.random.default_rng(seed)

    q_shell = rng.normal(
        0.0,
        shell_scale,
        7
    )

    # 64 arbitrary perturbations
    eta = rng.normal(
        0.0,
        1.0,
        N_STATES
    )

    P = np.zeros(
        (N_STATES, N_STATES),
        dtype=float
    )

    for G in states:

        logits = np.zeros(N_STATES)

        for a in states:

            delta = int(a ^ G)

            d = int(bits[delta].sum())

            logits[a] = beta * (
                q_shell[d]
                + epsilon * eta[delta]
            )

        P[G] = softmax(logits)

    return P


# ============================================================
# 7. NO-XOR CONTROL
#
# P(a|G) does not transform covariantly under XOR.
# ============================================================

def make_noXOR(
    seed=SEED,
    beta=1.0,
    scale=0.15
):

    rng = np.random.default_rng(seed)

    q = rng.normal(
        0.0,
        scale,
        N_STATES
    )

    P = np.zeros(
        (N_STATES, N_STATES)
    )

    for G in states:

        # G-dependent random perturbation
        local = rng.normal(
            0.0,
            scale,
            N_STATES
        )

        logits = beta * (q + local)

        P[G] = softmax(logits)

    return P


# ============================================================
# 8. NO-S6 CONTROL
#
# XOR covariance is retained:
#
# P(a|G)=K(a XOR G)
#
# but K(delta) is arbitrary over the 64 deltas.
#
# Thus XOR equivariance remains exact,
# while S6 invariance is broken.
# ============================================================

def make_noS6(
    seed=SEED,
    beta=1.0,
    scale=0.15
):

    rng = np.random.default_rng(seed)

    # arbitrary kernel on F2^6
    q = rng.normal(
        0.0,
        scale,
        N_STATES
    )

    P = np.zeros(
        (N_STATES, N_STATES)
    )

    for G in states:

        logits = np.zeros(N_STATES)

        for a in states:

            delta = int(a ^ G)

            logits[a] = beta * q[delta]

        P[G] = softmax(logits)

    return P


# ============================================================
# 9. NO-KERNEL
#
# Completely independent state-dependent distributions.
# ============================================================

def make_noKernel(
    seed=SEED,
    beta=1.0,
    scale=0.20
):

    rng = np.random.default_rng(seed)

    P = np.zeros(
        (N_STATES, N_STATES)
    )

    for G in states:

        logits = rng.normal(
            0.0,
            scale,
            N_STATES
        )

        P[G] = softmax(
            beta * logits
        )

    return P


# ============================================================
# 10. RANDOM-64 NULL
#
# Destroys the natural F2^6 / S6 correspondence.
# ============================================================

def make_random64(
    seed=SEED,
    beta=1.0,
    scale=0.20
):

    rng = np.random.default_rng(seed)

    P = np.zeros(
        (N_STATES, N_STATES)
    )

    base = rng.normal(
        0.0,
        scale,
        N_STATES
    )

    for G in states:

        permutation = rng.permutation(
            N_STATES
        )

        logits = base[permutation]

        P[G] = softmax(
            beta * logits
        )

    return P


# ============================================================
# 11. OPTIONAL: LOAD REAL QNN DATA
#
# If you have a CSV containing the real 64 probability
# distributions, it should contain:
#
# row = G
# columns = output probabilities
#
# exactly 64 rows x 64 columns.
#
# Set REAL_QNN_FILE to its filename.
# ============================================================

REAL_QNN_FILE = None

# Example:
#
# REAL_QNN_FILE = "P_QNN_real.csv"


def load_real_qnn(filename):

    df = pd.read_csv(filename)

    P = df.values.astype(float)

    if P.shape != (64, 64):

        raise ValueError(
            f"Expected 64x64 matrix, got {P.shape}"
        )

    # normalize defensively
    row_sums = P.sum(axis=1)

    if np.any(row_sums <= 0):

        raise ValueError(
            "At least one probability row has non-positive sum."
        )

    P = P / row_sums[:, None]

    return P


# ============================================================
# 12. JS DISTANCE MATRIX
# ============================================================

def js_distance_matrix(P):

    P = np.asarray(P, dtype=float)

    D = np.zeros(
        (N_STATES, N_STATES),
        dtype=float
    )

    for i in states:

        for j in range(
            int(i) + 1,
            N_STATES
        ):

            d = jensenshannon(
                P[i],
                P[j],
                base=2
            )

            D[i, j] = d
            D[j, i] = d

    return D


# ============================================================
# 13. CHECK PROBABILITY MATRIX
# ============================================================

def probability_audit(P):

    row_sums = P.sum(axis=1)

    return {
        "min_probability": float(P.min()),
        "max_probability": float(P.max()),
        "max_normalization_error":
            float(np.max(np.abs(row_sums - 1.0))),
        "mean_entropy_bits":
            float(
                -np.mean(
                    np.sum(
                        np.where(
                            P > 0,
                            P * np.log2(P),
                            0
                        ),
                        axis=1
                    )
                )
            )
    }


# ============================================================
# 14. PAIR TABLE
#
# 64 choose 2 = 2016
# ============================================================

def make_pair_table(D):

    rows = []

    for G in range(N_STATES):

        for Gp in range(
            G + 1,
            N_STATES
        ):

            delta = G ^ Gp

            d = int(
                bits[delta].sum()
            )

            rows.append({

                "G": G,

                "G_prime": Gp,

                "delta": delta,

                "d_H": d,

                "D": D[G, Gp]
            })

    df = pd.DataFrame(rows)

    assert len(df) == 2016

    return df


# ============================================================
# 15. SHELL AUDIT
# ============================================================

def shell_audit(pair_df):

    rows = []

    for d in range(1, 7):

        x = pair_df.loc[
            pair_df["d_H"] == d,
            "D"
        ].values

        shell_mean = np.mean(x)

        residuals = x - shell_mean

        rows.append({

            "d_H": d,

            "n_pairs": len(x),

            "mean_D": shell_mean,

            "std_D": np.std(x),

            "rmse_from_shell_mean":
                np.sqrt(
                    np.mean(
                        residuals ** 2
                    )
                ),

            "mae_from_shell_mean":
                np.mean(
                    np.abs(residuals)
                ),

            "min_D": np.min(x),

            "max_D": np.max(x),

            "range_D":
                np.max(x) - np.min(x)
        })

    return pd.DataFrame(rows)


# ============================================================
# 16. XOR AUDIT
#
# Test:
#
# D(G XOR v, G' XOR v) = D(G,G')
# ============================================================

def xor_audit(D):

    errors = []

    for G in states:

        for Gp in states:

            reference = D[G, Gp]

            for v in states:

                transformed = D[
                    int(G ^ v),
                    int(Gp ^ v)
                ]

                errors.append(
                    abs(
                        transformed
                        - reference
                    )
                )

    errors = np.asarray(errors)

    return {

        "xor_max_error":
            float(np.max(errors)),

        "xor_mean_error":
            float(np.mean(errors)),

        "xor_rms_error":
            float(
                np.sqrt(
                    np.mean(
                        errors ** 2
                    )
                )
            )
    }


# ============================================================
# 17. S6 AUDIT
#
# First reduce D to:
#
# F(delta)=D(0,delta)
#
# Then test:
#
# F(delta)=F(pi delta)
# ============================================================

def s6_orbit_audit(D):

    F = np.array([
        D[0, delta]
        for delta in states
    ])

    errors = []

    for delta in states:

        for p in S6:

            dp = permute_state(
                delta,
                p
            )

            errors.append(
                abs(
                    F[delta]
                    - F[dp]
                )
            )

    errors = np.asarray(errors)

    return {

        "S6_max_error":
            float(np.max(errors)),

        "S6_mean_error":
            float(np.mean(errors)),

        "S6_rms_error":
            float(
                np.sqrt(
                    np.mean(
                        errors ** 2
                    )
                )
            )
    }


# ============================================================
# 18. DIRECT S6 ORBIT TABLE
#
# For every Hamming shell, calculate variation among
# all S6-equivalent delta states.
# ============================================================

def s6_shell_audit(D):

    F = np.array([
        D[0, delta]
        for delta in states
    ])

    rows = []

    for d in range(7):

        deltas = [
            delta
            for delta in states
            if int(bits[delta].sum()) == d
        ]

        values = F[deltas]

        rows.append({

            "d_H": d,

            "orbit_size": len(values),

            "F_mean": np.mean(values),

            "F_std": np.std(values),

            "F_min": np.min(values),

            "F_max": np.max(values),

            "F_range":
                np.max(values)
                - np.min(values)
        })

    return pd.DataFrame(rows)


# ============================================================
# 19. HAMMING LAW
#
# Fit the six shell means:
#
# D_pred(G,G') = f(d_H(G,G'))
#
# Then measure residual.
# ============================================================

def hamming_law_audit(pair_df):

    shell_means = (
        pair_df
        .groupby("d_H")["D"]
        .mean()
        .to_dict()
    )

    predicted = np.array([
        shell_means[int(d)]
        for d in pair_df["d_H"]
    ])

    actual = pair_df["D"].values

    residual = (
        actual
        - predicted
    )

    r, p = pearsonr(
        pair_df["d_H"],
        actual
    )

    return {

        "pearson_r":
            float(r),

        "pearson_p":
            float(p),

        "RMSE":
            float(
                np.sqrt(
                    np.mean(
                        residual ** 2
                    )
                )
            ),

        "MAE":
            float(
                np.mean(
                    np.abs(
                        residual
                    )
                )
            ),

        "max_abs_error":
            float(
                np.max(
                    np.abs(
                        residual
                    )
                )
            )
    }


# ============================================================
# 20. STRONGER HAMMING TEST
#
# Pearson is NOT the main criterion.
#
# The decisive quantity is the within-shell dispersion.
# ============================================================

def hamming_geometry_score(shell_df):

    rmse = shell_df[
        "rmse_from_shell_mean"
    ].values

    return {

        "mean_shell_RMSE":
            float(np.mean(rmse)),

        "max_shell_RMSE":
            float(np.max(rmse)),

        "mean_shell_STD":
            float(
                np.mean(
                    shell_df["std_D"]
                )
            ),

        "max_shell_STD":
            float(
                np.max(
                    shell_df["std_D"]
                )
            )
    }


# ============================================================
# 21. COMPLETE AUDIT
# ============================================================

def full_audit(
    model_name,
    P
):

    print()
    print("=" * 78)
    print(model_name)
    print("=" * 78)

    # -----------------------------------------
    # Probability checks
    # -----------------------------------------

    prob = probability_audit(P)

    print()
    print("Probability audit")
    print(prob)

    # -----------------------------------------
    # Distance matrix
    # -----------------------------------------

    D = js_distance_matrix(P)

    # -----------------------------------------
    # Pair table
    # -----------------------------------------

    pairs = make_pair_table(D)

    # -----------------------------------------
    # Shells
    # -----------------------------------------

    shells = shell_audit(pairs)

    # -----------------------------------------
    # XOR
    # -----------------------------------------

    xor = xor_audit(D)

    print()
    print("XOR")
    print(xor)

    # -----------------------------------------
    # S6
    # -----------------------------------------

    s6 = s6_orbit_audit(D)

    print()
    print("S6")
    print(s6)

    # -----------------------------------------
    # S6 shell/orbit audit
    # -----------------------------------------

    s6_shells = s6_shell_audit(D)

    # -----------------------------------------
    # Hamming law
    # -----------------------------------------

    hamm = hamming_law_audit(
        pairs
    )

    print()
    print("Hamming law")
    print(hamm)

    # -----------------------------------------
    # Geometry score
    # -----------------------------------------

    geom = hamming_geometry_score(
        shells
    )

    print()
    print("Hamming geometry score")
    print(geom)

    # -----------------------------------------
    # Shell table
    # -----------------------------------------

    print()
    print("Shells")
    print(
        shells.to_string(
            index=False,
            float_format=lambda x:
                f"{x:.10g}"
        )
    )

    # -----------------------------------------
    # S6 shell table
    # -----------------------------------------

    print()
    print("S6 orbit shells")
    print(
        s6_shells.to_string(
            index=False,
            float_format=lambda x:
                f"{x:.10g}"
        )
    )

    # -----------------------------------------
    # Save
    # -----------------------------------------

    safe_name = model_name.replace(
        " ",
        "_"
    )

    pd.DataFrame(P).to_csv(
        os.path.join(
            OUTDIR,
            f"P_{safe_name}_64x64.csv"
        ),
        index=False
    )

    pd.DataFrame(D).to_csv(
        os.path.join(
            OUTDIR,
            f"D_{safe_name}_64x64.csv"
        ),
        index=False
    )

    pairs.to_csv(
        os.path.join(
            OUTDIR,
            f"pairs_{safe_name}_2016.csv"
        ),
        index=False
    )

    shells.to_csv(
        os.path.join(
            OUTDIR,
            f"shells_{safe_name}.csv"
        ),
        index=False
    )

    s6_shells.to_csv(
        os.path.join(
            OUTDIR,
            f"S6_shells_{safe_name}.csv"
        ),
        index=False
    )

    return {

        "P": P,

        "D": D,

        "pairs": pairs,

        "shells": shells,

        "S6_shells": s6_shells,

        "probability": prob,

        "xor": xor,

        "s6": s6,

        "hamming": hamm,

        "geometry": geom
    }


# ============================================================
# 22. BUILD ALL CONTROLS
# ============================================================

models = {}

# -----------------------------------------
# A. Exact D2
# -----------------------------------------

models["D2_EXACT"] = make_D2_exact(
    seed=SEED
)

# -----------------------------------------
# B. D2 with controlled symmetry breaking
# -----------------------------------------

models["D2_EPS_1E-4"] = make_D2_perturbed(
    seed=SEED,
    epsilon=1e-4
)

models["D2_EPS_1E-3"] = make_D2_perturbed(
    seed=SEED,
    epsilon=1e-3
)

models["D2_EPS_1E-2"] = make_D2_perturbed(
    seed=SEED,
    epsilon=1e-2
)

models["D2_EPS_1E-1"] = make_D2_perturbed(
    seed=SEED,
    epsilon=1e-1
)

# -----------------------------------------
# C. Controls
# -----------------------------------------

models["NO_XOR"] = make_noXOR(
    seed=SEED
)

models["NO_S6"] = make_noS6(
    seed=SEED
)

models["NO_KERNEL"] = make_noKernel(
    seed=SEED
)

models["RANDOM_64"] = make_random64(
    seed=SEED
)


# ============================================================
# 23. REAL QNN, IF PROVIDED
# ============================================================

if REAL_QNN_FILE is not None:

    print()
    print(
        "Loading real QNN:",
        REAL_QNN_FILE
    )

    models["REAL_QNN"] = load_real_qnn(
        REAL_QNN_FILE
    )


# ============================================================
# 24. RUN EVERYTHING
# ============================================================

results = {}

for name, P in models.items():

    results[name] = full_audit(
        name,
        P
    )


# ============================================================
# 25. SUMMARY TABLE
# ============================================================

summary_rows = []

for name, R in results.items():

    summary_rows.append({

        "model": name,

        "XOR_RMS":
            R["xor"]["xor_rms_error"],

        "XOR_MAX":
            R["xor"]["xor_max_error"],

        "S6_RMS":
            R["s6"]["S6_rms_error"],

        "S6_MAX":
            R["s6"]["S6_max_error"],

        "Hamming_RMSE":
            R["hamming"]["RMSE"],

        "Hamming_MAE":
            R["hamming"]["MAE"],

        "Hamming_MAX":
            R["hamming"]["max_abs_error"],

        "Shell_RMSE_mean":
            R["geometry"]["mean_shell_RMSE"],

        "Shell_RMSE_max":
            R["geometry"]["max_shell_RMSE"],

        "Pearson_r":
            R["hamming"]["pearson_r"],

        "Pearson_p":
            R["hamming"]["pearson_p"],

        "Entropy_bits":
            R["probability"]["mean_entropy_bits"]
    })


summary = pd.DataFrame(
    summary_rows
)

print()
print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print(
    summary.to_string(
        index=False,
        float_format=lambda x:
            f"{x:.8e}"
    )
)

summary.to_csv(
    os.path.join(
        OUTDIR,
        "FINAL_SUMMARY.csv"
    ),
    index=False
)


# ============================================================
# 26. FINAL THEORETICAL CHECK
# ============================================================

print()
print("=" * 78)
print("THEORETICAL INTERPRETATION")
print("=" * 78)

print(
    """
Exact group structure requires:

    XOR:
        D(G XOR v, G' XOR v) = D(G,G')

    S6:
        F(delta) = F(pi delta)

Together:

        D(G,G') = f(|G XOR G'|)

therefore:

        D(G,G') = f(d_H(G,G'))

The decisive numerical diagnostic is NOT Pearson r alone.

The strongest diagnostics are:

    XOR_RMS
    S6_RMS
    max within-shell STD
    Hamming RMSE

For an exact orbit construction these should be
close to floating-point precision (~1e-16).

A real QNN should be compared against these controls
without imposing the final Hamming distance directly.
"""
)


# ============================================================
# 27. CHECK 2016 PAIRS
# ============================================================

print()
print("=" * 78)
print("PAIR COUNT CHECK")
print("=" * 78)

print(
    "Expected pairs:",
    64 * 63 // 2
)

for name, R in results.items():

    print(
        f"{name:20s}:",
        len(R["pairs"])
    )


# ============================================================
# 28. CHECK HAMMING SHELL CARDINALITIES
# ============================================================

expected_shells = {
    1: 64 * 6,
    2: 64 * 15,
    3: 64 * 20,
    4: 64 * 15,
    5: 64 * 6,
    6: 64 * 1
}

print()
print("=" * 78)
print("HAMMING SHELL CARDINALITIES")
print("=" * 78)

for d, expected in expected_shells.items():

    print(
        f"d={d}: expected {expected}"
    )


# ============================================================
# 29. DONE
# ============================================================

print()
print("=" * 78)
print("AUDIT COMPLETE")
print("=" * 78)

print(
    "Output directory:",
    os.path.abspath(OUTDIR)
)

print()
print("Generated files include:")
print("  D_Hamming_64x64.csv")
print("  P_<model>_64x64.csv")
print("  D_<model>_64x64.csv")
print("  pairs_<model>_2016.csv")
print("  shells_<model>.csv")
print("  S6_shells_<model>.csv")
print("  FINAL_SUMMARY.csv")
