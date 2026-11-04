"""
Dockerfile56.py
============================================================
Q-PROOF L4.6 — NESTED / COUPLED MONTE CARLO CRITICALITY
============================================================

L4.5 found approximately:

    fc ~ 15%

L4.6 improves statistical robustness by using COUPLED
NESTED FAILURE REALIZATIONS.

For every Monte Carlo seed:

    V_1, V_2, ..., V_N

is one random permutation of validators.

Fault levels are then nested:

    0%  -> {}
    1%  -> {V_1}
    2%  -> {V_1,V_2}
    3%  -> {V_1,V_2,V_3}
    ...

Therefore:

    F(f_i) subset F(f_j)       for f_i < f_j

This removes much of the artificial point-to-point noise
introduced by independently sampled fault sets.

Outputs
-------

    lambda2(f)
    P_partition(f)
    bootstrap 95% CI
    f_c at P_partition >= 0.50
    lambda2 warning threshold
    monotonicity diagnostics
    Pearson correlation
    Spearman correlation
    shard vulnerability
    CSV export
    PNG plots

============================================================
"""

import csv
import random
from dataclasses import dataclass

import numpy as np
import networkx as nx

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# ============================================================
# CONFIGURATION
# ============================================================

NUM_SHARDS = 4
VALIDATORS_PER_SHARD = 12

TOTAL_VALIDATORS = (
    NUM_SHARDS * VALIDATORS_PER_SHARD
)

MASTER_SEED = 56

# Monte Carlo
N_SEEDS = 1000

# Fault scan
FAULT_MIN = 0.00
FAULT_MAX = 0.30
FAULT_STEP = 0.01

FAULT_LEVELS = np.arange(
    FAULT_MIN,
    FAULT_MAX + 1e-9,
    FAULT_STEP
)

# Bootstrap
N_BOOTSTRAP = 5000

CONFIDENCE = 0.95

# Topological warning threshold
LAMBDA2_WARNING = 0.05

# Critical partition probability
CRITICAL_PROBABILITY = 0.50

# Output
CSV_FILE = (
    "qproof_l46_nested_results.csv"
)

PLOT_FILE = (
    "qproof_l46_criticality.png"
)


# ============================================================
# VALIDATOR
# ============================================================

@dataclass
class Validator:

    validator_id: int
    shard_id: int


# ============================================================
# NETWORK
# ============================================================

class QProofNetwork:

    def __init__(self):

        self.validators = []

        validator_id = 0

        for shard_id in range(
            NUM_SHARDS
        ):

            for _ in range(
                VALIDATORS_PER_SHARD
            ):

                self.validators.append(
                    Validator(
                        validator_id=validator_id,
                        shard_id=shard_id
                    )
                )

                validator_id += 1

        self.graph = nx.Graph()

        self.bridge_edges = []

        self.build_topology()

    # ========================================================
    # TOPOLOGY
    # ========================================================

    def build_topology(self):

        G = nx.Graph()

        # Nodes
        for v in self.validators:

            G.add_node(
                v.validator_id,
                shard=v.shard_id
            )

        # ----------------------------------------------------
        # Complete intra-shard connectivity
        # ----------------------------------------------------

        for shard_id in range(
            NUM_SHARDS
        ):

            nodes = [
                v.validator_id
                for v in self.validators
                if v.shard_id == shard_id
            ]

            for i in range(
                len(nodes)
            ):

                for j in range(
                    i + 1,
                    len(nodes)
                ):

                    G.add_edge(
                        nodes[i],
                        nodes[j],
                        weight=1.0,
                        edge_type="intra"
                    )

        # ----------------------------------------------------
        # Cross-shard bridges
        # ----------------------------------------------------

        for shard_a in range(
            NUM_SHARDS
        ):

            for shard_b in range(
                shard_a + 1,
                NUM_SHARDS
            ):

                A = [
                    v.validator_id
                    for v in self.validators
                    if v.shard_id == shard_a
                ]

                B = [
                    v.validator_id
                    for v in self.validators
                    if v.shard_id == shard_b
                ]

                edge = (
                    A[0],
                    B[0]
                )

                G.add_edge(
                    edge[0],
                    edge[1],
                    weight=0.5,
                    edge_type="bridge"
                )

                self.bridge_edges.append(
                    edge
                )

        self.graph = G

    # ========================================================
    # ACTIVE GRAPH
    # ========================================================

    def active_graph(
        self,
        failed_ids
    ):

        failed_ids = set(
            failed_ids
        )

        active = [
            v.validator_id
            for v in self.validators
            if v.validator_id
            not in failed_ids
        ]

        return self.graph.subgraph(
            active
        ).copy()

    # ========================================================
    # CONNECTIVITY
    # ========================================================

    def is_connected(
        self,
        failed_ids
    ):

        G = self.active_graph(
            failed_ids
        )

        if G.number_of_nodes() < 2:

            return False

        return nx.is_connected(G)

    # ========================================================
    # LAMBDA 2
    # ========================================================

    def lambda2(
        self,
        failed_ids
    ):

        G = self.active_graph(
            failed_ids
        )

        if G.number_of_nodes() < 2:

            return 0.0

        if not nx.is_connected(G):

            return 0.0

        L = nx.laplacian_matrix(
            G,
            weight="weight"
        ).toarray()

        eigenvalues = np.linalg.eigvalsh(
            L
        )

        eigenvalues.sort()

        return float(
            eigenvalues[1]
        )

    # ========================================================
    # SHARD FAILURE COUNTS
    # ========================================================

    def shard_failure_counts(
        self,
        failed_ids
    ):

        failed_ids = set(
            failed_ids
        )

        counts = [
            0
            for _ in range(
                NUM_SHARDS
            )
        ]

        for v in self.validators:

            if v.validator_id in failed_ids:

                counts[v.shard_id] += 1

        return counts


# ============================================================
# RESULT STRUCTURE
# ============================================================

@dataclass
class LevelResult:

    fault_fraction: float

    failed_count: int

    active_count: float

    lambda_mean: float

    lambda_ci_low: float

    lambda_ci_high: float

    lambda_std: float

    partition_probability: float

    partition_ci_low: float

    partition_ci_high: float

    warning_probability: float

    partition_count: int

    total_runs: int


# ============================================================
# BOOTSTRAP
# ============================================================

def bootstrap_ci(
    values,
    statistic=np.mean,
    n_bootstrap=N_BOOTSTRAP,
    confidence=CONFIDENCE,
    seed=56
):

    values = np.asarray(
        values,
        dtype=float
    )

    n = len(values)

    if n == 0:

        return (
            np.nan,
            np.nan,
            np.nan
        )

    rng = np.random.default_rng(
        seed
    )

    indices = rng.integers(
        0,
        n,
        size=(
            n_bootstrap,
            n
        )
    )

    samples = values[
        indices
    ]

    statistics = np.array([
        statistic(sample)
        for sample in samples
    ])

    alpha = (
        1.0
        - confidence
    )

    low = np.quantile(
        statistics,
        alpha / 2
    )

    high = np.quantile(
        statistics,
        1.0 - alpha / 2
    )

    return (
        float(statistic(values)),
        float(low),
        float(high)
    )


# ============================================================
# NESTED FAULT LEVELS
# ============================================================

def failed_count_for_fraction(
    fraction
):

    return int(
        round(
            fraction
            * TOTAL_VALIDATORS
        )
    )


# ============================================================
# COUPLED MONTE CARLO
# ============================================================

def run_nested_monte_carlo():

    network = QProofNetwork()

    # --------------------------------------------------------
    # Containers
    # --------------------------------------------------------

    lambda_samples = {
        float(f): []
        for f in FAULT_LEVELS
    }

    partition_samples = {
        float(f): []
        for f in FAULT_LEVELS
    }

    active_samples = {
        float(f): []
        for f in FAULT_LEVELS
    }

    shard_failure_samples = {
        shard: []
        for shard in range(
            NUM_SHARDS
        )
    }

    # --------------------------------------------------------
    # Seed loop
    # --------------------------------------------------------

    for seed_index in range(
        N_SEEDS
    ):

        rng = random.Random(
            MASTER_SEED
            + seed_index
        )

        # ----------------------------------------------------
        # ONE permutation per seed
        # ----------------------------------------------------

        permutation = [
            v.validator_id
            for v in network.validators
        ]

        rng.shuffle(
            permutation
        )

        # ----------------------------------------------------
        # Nested fault states
        # ----------------------------------------------------

        for fault in FAULT_LEVELS:

            fault = float(
                fault
            )

            failed_count = (
                failed_count_for_fraction(
                    fault
                )
            )

            failed_ids = permutation[
                :failed_count
            ]

            lam = network.lambda2(
                failed_ids
            )

            connected = (
                network.is_connected(
                    failed_ids
                )
            )

            partition = (
                not connected
            )

            warning = (
                lam
                < LAMBDA2_WARNING
            )

            active = (
                TOTAL_VALIDATORS
                - failed_count
            )

            lambda_samples[
                fault
            ].append(
                lam
            )

            partition_samples[
                fault
            ].append(
                float(partition)
            )

            active_samples[
                fault
            ].append(
                active
            )

        # ----------------------------------------------------
        # Shard vulnerability
        # ----------------------------------------------------

        # Evaluate at 15%
        critical_count = (
            failed_count_for_fraction(
                0.15
            )
        )

        critical_failed = permutation[
            :critical_count
        ]

        shard_counts = (
            network.shard_failure_counts(
                critical_failed
            )
        )

        for shard_id in range(
            NUM_SHARDS
        ):

            shard_failure_samples[
                shard_id
            ].append(
                shard_counts[shard_id]
            )

    # ========================================================
    # AGGREGATION
    # ========================================================

    results = []

    for fault in FAULT_LEVELS:

        fault = float(
            fault
        )

        lambdas = np.asarray(
            lambda_samples[fault]
        )

        partitions = np.asarray(
            partition_samples[fault]
        )

        actives = np.asarray(
            active_samples[fault]
        )

        (
            lambda_mean,
            lambda_ci_low,
            lambda_ci_high
        ) = bootstrap_ci(
            lambdas,
            statistic=np.mean,
            seed=MASTER_SEED
            + int(fault * 1000)
        )

        (
            partition_probability,
            partition_ci_low,
            partition_ci_high
        ) = bootstrap_ci(
            partitions,
            statistic=np.mean,
            seed=MASTER_SEED
            + 10000
            + int(fault * 1000)
        )

        warning_probability = float(
            np.mean(
                lambdas
                < LAMBDA2_WARNING
            )
        )

        result = LevelResult(

            fault_fraction=fault,

            failed_count=
                failed_count_for_fraction(
                    fault
                ),

            active_count=float(
                np.mean(
                    actives
                )
            ),

            lambda_mean=lambda_mean,

            lambda_ci_low=lambda_ci_low,

            lambda_ci_high=lambda_ci_high,

            lambda_std=float(
                np.std(
                    lambdas,
                    ddof=1
                )
            ),

            partition_probability=
                partition_probability,

            partition_ci_low=
                partition_ci_low,

            partition_ci_high=
                partition_ci_high,

            warning_probability=
                warning_probability,

            partition_count=int(
                partitions.sum()
            ),

            total_runs=N_SEEDS
        )

        results.append(
            result
        )

    return (
        results,
        shard_failure_samples
    )


# ============================================================
# CRITICAL THRESHOLD
# ============================================================

def interpolate_threshold(
    results,
    target=0.50
):

    for i in range(
        len(results) - 1
    ):

        r1 = results[i]
        r2 = results[i + 1]

        p1 = r1.partition_probability
        p2 = r2.partition_probability

        if (
            p1 < target
            and p2 >= target
        ):

            f1 = r1.fault_fraction
            f2 = r2.fault_fraction

            if abs(
                p2 - p1
            ) < 1e-12:

                return f2

            fraction = (
                f1
                + (
                    target - p1
                )
                / (
                    p2 - p1
                )
                * (
                    f2 - f1
                )
            )

            return float(
                fraction
            )

    return None


# ============================================================
# WARNING THRESHOLD
# ============================================================

def lambda_warning_threshold(
    results,
    target=0.50
):

    for i in range(
        len(results) - 1
    ):

        r1 = results[i]
        r2 = results[i + 1]

        p1 = r1.warning_probability
        p2 = r2.warning_probability

        if (
            p1 < target
            and p2 >= target
        ):

            f1 = r1.fault_fraction
            f2 = r2.fault_fraction

            if abs(
                p2 - p1
            ) < 1e-12:

                return f2

            return float(
                f1
                + (
                    target - p1
                )
                / (
                    p2 - p1
                )
                * (
                    f2 - f1
                )
            )

    return None


# ============================================================
# MONOTONICITY TEST
# ============================================================

def monotonicity_diagnostics(
    results
):

    probabilities = np.array([
        r.partition_probability
        for r in results
    ])

    lambdas = np.array([
        r.lambda_mean
        for r in results
    ])

    # --------------------------------------------------------
    # Count decreases in partition probability
    # --------------------------------------------------------

    probability_deltas = np.diff(
        probabilities
    )

    probability_violations = int(
        np.sum(
            probability_deltas < -1e-12
        )
    )

    # --------------------------------------------------------
    # Count increases in lambda2
    # --------------------------------------------------------

    lambda_deltas = np.diff(
        lambdas
    )

    lambda_violations = int(
        np.sum(
            lambda_deltas > 1e-12
        )
    )

    # --------------------------------------------------------
    # Spearman
    # --------------------------------------------------------

    try:

        from scipy.stats import spearmanr

        spearman = float(
            spearmanr(
                lambdas,
                probabilities
            ).statistic
        )

    except Exception:

        # rank correlation fallback
        lambda_rank = (
            np.argsort(
                np.argsort(
                    lambdas
                )
            )
        )

        probability_rank = (
            np.argsort(
                np.argsort(
                    probabilities
                )
            )
        )

        spearman = float(
            np.corrcoef(
                lambda_rank,
                probability_rank
            )[0, 1]
        )

    pearson = float(
        np.corrcoef(
            lambdas,
            probabilities
        )[0, 1]
    )

    return {
        "partition_probability_violations":
            probability_violations,

        "lambda2_monotonicity_violations":
            lambda_violations,

        "pearson":
            pearson,

        "spearman":
            spearman
    }


# ============================================================
# PRINT TABLE
# ============================================================

def print_results(
    results
):

    print()
    print("=" * 120)
    print(
        "Q-PROOF L4.6 — NESTED MONTE CARLO"
    )
    print("=" * 120)

    print(
        "Fault | Active | Lambda2 | "
        "95% CI Lambda2 | "
        "P(partition) | 95% CI | "
        "P(lambda2<thr)"
    )

    print("-" * 120)

    for r in results:

        print(
            f"{r.fault_fraction*100:5.1f}% | "
            f"{r.active_count:6.1f} | "
            f"{r.lambda_mean:7.4f} | "
            f"[{r.lambda_ci_low:.4f},"
            f"{r.lambda_ci_high:.4f}] | "
            f"{r.partition_probability:12.4f} | "
            f"[{r.partition_ci_low:.3f},"
            f"{r.partition_ci_high:.3f}] | "
            f"{r.warning_probability:.4f}"
        )


# ============================================================
# SHARD REPORT
# ============================================================

def print_shard_report(
    shard_failure_samples
):

    print()
    print("=" * 80)
    print(
        "SHARD VULNERABILITY @ 15% FAULT"
    )
    print("=" * 80)

    for shard_id in range(
        NUM_SHARDS
    ):

        values = np.asarray(
            shard_failure_samples[
                shard_id
            ]
        )

        (
            mean,
            low,
            high
        ) = bootstrap_ci(
            values,
            statistic=np.mean,
            seed=1000 + shard_id
        )

        print(
            f"Shard {shard_id}: "
            f"mean failed = {mean:.3f}, "
            f"95% CI = "
            f"[{low:.3f}, {high:.3f}]"
        )


# ============================================================
# CSV EXPORT
# ============================================================

def export_csv(
    results,
    filename=CSV_FILE
):

    with open(
        filename,
        "w",
        newline=""
    ) as f:

        writer = csv.writer(
            f
        )

        writer.writerow([
            "fault_fraction",
            "failed_count",
            "active_count",
            "lambda2_mean",
            "lambda2_ci_low",
            "lambda2_ci_high",
            "lambda2_std",
            "partition_probability",
            "partition_ci_low",
            "partition_ci_high",
            "warning_probability",
            "partition_count",
            "total_runs"
        ])

        for r in results:

            writer.writerow([
                r.fault_fraction,
                r.failed_count,
                r.active_count,
                r.lambda_mean,
                r.lambda_ci_low,
                r.lambda_ci_high,
                r.lambda_std,
                r.partition_probability,
                r.partition_ci_low,
                r.partition_ci_high,
                r.warning_probability,
                r.partition_count,
                r.total_runs
            ])

    print()
    print(
        f"[CSV Export] {filename}"
    )


# ============================================================
# PLOT
# ============================================================

def create_plot(
    results,
    filename=PLOT_FILE
):

    if not MATPLOTLIB_AVAILABLE:

        print(
            "[Plot] matplotlib not installed."
        )

        return

    fault = np.array([
        r.fault_fraction * 100
        for r in results
    ])

    probability = np.array([
        r.partition_probability
        for r in results
    ])

    probability_low = np.array([
        r.partition_ci_low
        for r in results
    ])

    probability_high = np.array([
        r.partition_ci_high
        for r in results
    ])

    lambda2 = np.array([
        r.lambda_mean
        for r in results
    ])

    lambda_low = np.array([
        r.lambda_ci_low
        for r in results
    ])

    lambda_high = np.array([
        r.lambda_ci_high
        for r in results
    ])

    fig, ax1 = plt.subplots(
        figsize=(11, 7)
    )

    ax1.plot(
        fault,
        probability,
        marker="o",
        label="P(partition)"
    )

    ax1.fill_between(
        fault,
        probability_low,
        probability_high,
        alpha=0.20
    )

    ax1.axhline(
        0.50,
        linestyle="--",
        label="Critical probability = 0.50"
    )

    ax1.set_xlabel(
        "Failed validators (%)"
    )

    ax1.set_ylabel(
        "Partition probability"
    )

    ax1.set_ylim(
        -0.02,
        1.02
    )

    ax2 = ax1.twinx()

    ax2.plot(
        fault,
        lambda2,
        marker="s",
        label="Lambda_2"
    )

    ax2.fill_between(
        fault,
        lambda_low,
        lambda_high,
        alpha=0.15
    )

    ax2.axhline(
        LAMBDA2_WARNING,
        linestyle=":",
        label="Lambda_2 warning threshold"
    )

    ax2.set_ylabel(
        "Algebraic connectivity Lambda_2"
    )

    ax1.set_title(
        "Q-PROOF L4.6 — Nested Monte Carlo Criticality"
    )

    fig.tight_layout()

    fig.savefig(
        filename,
        dpi=200
    )

    plt.close(
        fig
    )

    print(
        f"[Plot Export] {filename}"
    )


# ============================================================
# FINAL REPORT
# ============================================================

def final_report(
    results
):

    diagnostics = (
        monotonicity_diagnostics(
            results
        )
    )

    fc = interpolate_threshold(
        results,
        CRITICAL_PROBABILITY
    )

    fw = lambda_warning_threshold(
        results,
        CRITICAL_PROBABILITY
    )

    print()
    print("=" * 80)
    print(
        "L4.6 CRITICALITY REPORT"
    )
    print("=" * 80)

    if fc is not None:

        print(
            f"Estimated f_c "
            f"(P_partition=0.50): "
            f"{fc*100:.3f}%"
        )

    else:

        print(
            "Estimated f_c: N/A"
        )

    if fw is not None:

        print(
            f"Lambda2 warning threshold "
            f"(P_warning=0.50): "
            f"{fw*100:.3f}%"
        )

    else:

        print(
            "Lambda2 warning threshold: N/A"
        )

    print()

    print(
        "Partition probability "
        "monotonicity violations: "
        f"{diagnostics['partition_probability_violations']}"
    )

    print(
        "Lambda2 monotonicity violations: "
        f"{diagnostics['lambda2_monotonicity_violations']}"
    )

    print()

    print(
        f"Pearson(lambda2, P_partition): "
        f"{diagnostics['pearson']:.8f}"
    )

    print(
        f"Spearman(lambda2, P_partition): "
        f"{diagnostics['spearman']:.8f}"
    )

    print()

    if (
        diagnostics[
            "partition_probability_violations"
        ] == 0
    ):

        print(
            "[Nested Monotonicity PASS]"
        )

    else:

        print(
            "[Nested Monotonicity INFO]"
        )

    if (
        diagnostics[
            "lambda2_monotonicity_violations"
        ] == 0
    ):

        print(
            "[Lambda2 Monotonicity PASS]"
        )

    else:

        print(
            "[Lambda2 Monotonicity INFO]"
        )

    if fc is not None:

        print(
            "[Criticality Estimation PASS]"
        )

    else:

        print(
            "[Criticality Estimation WARNING]"
        )

    print("=" * 80)

    return fc, fw, diagnostics


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 80)
    print(
        "Q-PROOF L4.6 — NESTED / COUPLED MONTE CARLO"
    )
    print("=" * 80)

    print(
        f"Validators: "
        f"{TOTAL_VALIDATORS}"
    )

    print(
        f"Shards: "
        f"{NUM_SHARDS}"
    )

    print(
        f"Validators/shard: "
        f"{VALIDATORS_PER_SHARD}"
    )

    print(
        f"Seeds: "
        f"{N_SEEDS}"
    )

    print(
        f"Bootstrap: "
        f"{N_BOOTSTRAP}"
    )

    print(
        f"Fault range: "
        f"{FAULT_MIN*100:.1f}%"
        f" -> "
        f"{FAULT_MAX*100:.1f}%"
    )

    print(
        f"Fault step: "
        f"{FAULT_STEP*100:.1f}%"
    )

    # --------------------------------------------------------
    # Baseline
    # --------------------------------------------------------

    baseline_network = (
        QProofNetwork()
    )

    baseline_lambda = (
        baseline_network.lambda2(
            []
        )
    )

    print()
    print(
        f"Baseline Lambda_2: "
        f"{baseline_lambda:.6f}"
    )

    # --------------------------------------------------------
    # Monte Carlo
    # --------------------------------------------------------

    results, shard_samples = (
        run_nested_monte_carlo()
    )

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    print_results(
        results
    )

    # --------------------------------------------------------
    # Criticality
    # --------------------------------------------------------

    fc, fw, diagnostics = (
        final_report(
            results
        )
    )

    # --------------------------------------------------------
    # Shards
    # --------------------------------------------------------

    print_shard_report(
        shard_samples
    )

    # --------------------------------------------------------
    # Export
    # --------------------------------------------------------

    export_csv(
        results
    )

    create_plot(
        results
    )

    # --------------------------------------------------------
    # Final
    # --------------------------------------------------------

    print()
    print("=" * 80)
    print(
        "Q-PROOF L4.6 FINAL RESULT"
    )
    print("=" * 80)

    print(
        f"Baseline Lambda_2: "
        f"{baseline_lambda:.6f}"
    )

    if fc is not None:

        print(
            f"Critical fault fraction: "
            f"{fc*100:.3f}%"
        )

    else:

        print(
            "Critical fault fraction: N/A"
        )

    if fw is not None:

        print(
            f"Lambda2 warning fraction: "
            f"{fw*100:.3f}%"
        )

    else:

        print(
            "Lambda2 warning fraction: N/A"
        )

    print(
        f"Pearson correlation: "
        f"{diagnostics['pearson']:.8f}"
    )

    print(
        f"Spearman correlation: "
        f"{diagnostics['spearman']:.8f}"
    )

    print()

    print(
        "[L4.6 Monte Carlo PASS]"
    )

    print(
        "Coupled nested fault realizations "
        "completed successfully."
    )

    print("=" * 80)


if __name__ == "__main__":

    main()
