# -*- coding: utf-8 -*-
"""
NARESSI–ZULIANI PHASE I — cleaned / corrected version

Obiettivo:
training -> pesi W -> spettro -> coarse-graining -> ricorrenza -> test Phi

Nota metodologica:
- Le quantità spettrali sono calcolate direttamente dai valori singolari.
- Il test S^2 = a S + b I è separato dal test dell'ipotesi Fibonacci:
  a = 1, b = 1 NON viene assunto dal fit.
- Il coarse-graining usa una SVD troncata:
      W_r = U_r^T W V_r
  che è dimensionalmente ben definita anche per matrici rettangolari.
- Per l'RG confrontiamo matrici normalizzate, così la semplice crescita/decrescita
  della norma non viene scambiata per struttura universale.
"""

import os
import copy
import random
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
from torchvision import datasets, transforms

from scipy.linalg import svd
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")

# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

SEED = 42
N_SEEDS = 3          # aumentare a 10/30 solo dopo il test
BATCH_SIZE = 256
EPOCHS = 10
LR = 1e-3

RG_STEPS = 3
RG_REDUCTION = 0.5
EPS = 1e-12

DATA_DIR = "/kaggle/working/data"
OUT_DIR = "/kaggle/working/phase1_results_fixed"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

PHI = (1 + np.sqrt(5.0)) / 2
PHI_INV = 1.0 / PHI

print("=" * 78)
print("NARESSI–ZULIANI PHASE I — FIXED")
print("=" * 78)
print("Device:", DEVICE)
print("Phi:", PHI)
print("Phi^-1:", PHI_INV)
print("Output:", OUT_DIR)


# ============================================================
# 2. RIPRODUCIBILITÀ
# ============================================================

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # Determinismo dove supportato.
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


set_seed(SEED)


# ============================================================
# 3. DATASET MNIST
# ============================================================

transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.1307,), (0.3081,))
])

train_full = datasets.MNIST(
    root=DATA_DIR,
    train=True,
    download=True,
    transform=transform
)

test_dataset = datasets.MNIST(
    root=DATA_DIR,
    train=False,
    download=True,
    transform=transform
)

train_size = 55_000
val_size = len(train_full) - train_size

train_dataset, val_dataset = random_split(
    train_full,
    [train_size, val_size],
    generator=torch.Generator().manual_seed(SEED)
)

# num_workers=0 evita differenze di comportamento fra ambienti Kaggle/locali.
train_loader = DataLoader(
    train_dataset,
    batch_size=BATCH_SIZE,
    shuffle=True,
    num_workers=0,
    pin_memory=torch.cuda.is_available()
)

val_loader = DataLoader(
    val_dataset,
    batch_size=BATCH_SIZE,
    shuffle=False,
    num_workers=0,
    pin_memory=torch.cuda.is_available()
)

test_loader = DataLoader(
    test_dataset,
    batch_size=BATCH_SIZE,
    shuffle=False,
    num_workers=0,
    pin_memory=torch.cuda.is_available()
)

print("Training:", len(train_dataset))
print("Validation:", len(val_dataset))
print("Test:", len(test_dataset))


# ============================================================
# 4. MODELLO
# ============================================================

class MLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(784, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 10)

    def forward(self, x):
        x = x.view(x.size(0), -1)
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return self.fc3(x)


criterion = nn.CrossEntropyLoss()


# ============================================================
# 5. EVALUATION / TRAINING
# ============================================================

def evaluate(model, loader):
    model.eval()
    total_loss = 0.0
    correct = 0
    total = 0

    with torch.no_grad():
        for x, y in loader:
            x = x.to(DEVICE, non_blocking=True)
            y = y.to(DEVICE, non_blocking=True)

            logits = model(x)
            loss = criterion(logits, y)

            total_loss += loss.item() * x.size(0)
            correct += (logits.argmax(dim=1) == y).sum().item()
            total += x.size(0)

    return total_loss / total, correct / total


def train_model(seed, epochs=EPOCHS):
    set_seed(seed)

    model = MLP().to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)

    history = []
    checkpoints = []

    for epoch in range(1, epochs + 1):
        model.train()
        running_loss = 0.0
        total = 0

        for x, y in train_loader:
            x = x.to(DEVICE, non_blocking=True)
            y = y.to(DEVICE, non_blocking=True)

            optimizer.zero_grad(set_to_none=True)

            logits = model(x)
            loss = criterion(logits, y)

            loss.backward()
            optimizer.step()

            running_loss += loss.item() * x.size(0)
            total += x.size(0)

        train_loss = running_loss / total
        val_loss, val_acc = evaluate(model, val_loader)

        checkpoints.append({
            "epoch": epoch,
            "state_dict": copy.deepcopy(model.state_dict())
        })

        history.append({
            "seed": seed,
            "epoch": epoch,
            "train_loss": train_loss,
            "val_loss": val_loss,
            "val_accuracy": val_acc
        })

        print(
            f"Seed {seed:03d} | Epoch {epoch:02d}/{epochs} | "
            f"loss={train_loss:.4f} | val={val_acc:.4f}"
        )

    test_loss, test_acc = evaluate(model, test_loader)
    print(f"Seed {seed:03d} | TEST accuracy={test_acc:.4f}")

    return model, pd.DataFrame(history), checkpoints


# ============================================================
# 6. REPLICHE
# ============================================================

all_runs = []

for seed in range(SEED, SEED + N_SEEDS):
    print("\n" + "=" * 78)
    print("TRAINING SEED:", seed)
    print("=" * 78)

    model, history, checkpoints = train_model(seed, EPOCHS)

    all_runs.append({
        "seed": seed,
        "model": model,
        "history": history,
        "checkpoints": checkpoints
    })


# ============================================================
# 7. ESTRAZIONE PESI
# ============================================================

def get_weight_matrices(state_dict):
    return {
        "fc1": state_dict["fc1.weight"].detach().cpu().numpy(),
        "fc2": state_dict["fc2.weight"].detach().cpu().numpy(),
        "fc3": state_dict["fc3.weight"].detach().cpu().numpy(),
    }


# ============================================================
# 8. METRICHE SPETTRALI
# ============================================================

def safe_normalize(values):
    values = np.asarray(values, dtype=np.float64)
    total = values.sum()

    if total < EPS:
        return np.ones_like(values) / len(values)

    return values / total


def spectral_metrics(W):
    """
    Per W rettangolare:
        sigma_i = singular values
        lambda_i = sigma_i^2

    rho_gram = lambda_max = sigma_max^2
    spectral_norm = sigma_max
    pr_sigma = (sum sigma^2)^2 / sum sigma^4

    effective_rank:
        exp(H(p)), p_i = sigma_i / sum_j sigma_j
    """

    W = np.asarray(W, dtype=np.float64)

    singular_values = svd(W, compute_uv=False)
    singular_values = np.maximum(singular_values, 0.0)

    sigma2 = singular_values ** 2
    sigma4 = singular_values ** 4

    if len(singular_values) == 0:
        raise ValueError("Empty matrix.")

    rho = float(sigma2[0])

    if len(sigma2) >= 2 and sigma2[0] > EPS:
        ratio = max(float(sigma2[1] / sigma2[0]), EPS)
        kappa = float(-np.log(ratio))
    else:
        kappa = np.nan

    if sigma4.sum() > EPS:
        pr_sigma = float((sigma2.sum() ** 2) / sigma4.sum())
    else:
        pr_sigma = np.nan

    p = safe_normalize(singular_values)
    p = p[p > 0]
    effective_rank = float(np.exp(-np.sum(p * np.log(p))))

    return {
        "rho_gram": rho,
        "kappa_gram": kappa,
        "pr_sigma": pr_sigma,
        "effective_rank": effective_rank,
        "spectral_norm": float(singular_values[0]),
        "frobenius_norm": float(np.linalg.norm(W, ord="fro")),
        "stable_rank": float(
            np.sum(sigma2) / max(sigma2[0], EPS)
        ),
        "rank": int(np.linalg.matrix_rank(W)),
        "n_rows": W.shape[0],
        "n_cols": W.shape[1],
    }


def analyze_checkpoints(checkpoints, seed, model_label):
    rows = []

    for checkpoint in checkpoints:
        epoch = checkpoint["epoch"]
        weights = get_weight_matrices(checkpoint["state_dict"])

        for layer_name, W in weights.items():
            metrics = spectral_metrics(W)
            metrics.update({
                "seed": seed,
                "model": model_label,
                "epoch": epoch,
                "layer": layer_name
            })
            rows.append(metrics)

    return pd.DataFrame(rows)


# ============================================================
# 9. TRAINED / UNTRAINED / SHUFFLED CONTROLS
# ============================================================

trained_tables = [
    analyze_checkpoints(
        run["checkpoints"], run["seed"], "trained"
    )
    for run in all_runs
]

trained_metrics = pd.concat(trained_tables, ignore_index=True)
trained_metrics.to_csv(
    os.path.join(OUT_DIR, "trained_spectral_metrics.csv"),
    index=False
)

untrained_tables = []

for run in all_runs:
    set_seed(run["seed"])
    model = MLP().to(DEVICE)

    checkpoints = [{
        "epoch": 0,
        "state_dict": copy.deepcopy(model.state_dict())
    }]

    untrained_tables.append(
        analyze_checkpoints(
            checkpoints, run["seed"], "untrained"
        )
    )

untrained_metrics = pd.concat(
    untrained_tables, ignore_index=True
)

untrained_metrics.to_csv(
    os.path.join(OUT_DIR, "untrained_spectral_metrics.csv"),
    index=False
)


def shuffled_state_dict(state_dict, seed=123):
    """
    Permutazione indipendente degli elementi di ciascun kernel.
    I bias non vengono modificati.
    """
    rng = np.random.default_rng(seed)
    new_state = copy.deepcopy(state_dict)

    for key, tensor in state_dict.items():
        if not key.endswith(".weight"):
            continue

        arr = tensor.detach().cpu().numpy().copy()
        flat = arr.reshape(-1)
        rng.shuffle(flat)

        new_state[key] = torch.as_tensor(
            flat.reshape(arr.shape),
            dtype=tensor.dtype
        )

    return new_state


shuffled_tables = []

for run in all_runs:
    state = run["checkpoints"][-1]["state_dict"]

    shuffled_state = shuffled_state_dict(
        state,
        seed=run["seed"] + 10_000
    )

    shuffled_tables.append(
        analyze_checkpoints(
            [{
                "epoch": EPOCHS,
                "state_dict": shuffled_state
            }],
            run["seed"],
            "shuffled"
        )
    )

shuffled_metrics = pd.concat(
    shuffled_tables, ignore_index=True
)

shuffled_metrics.to_csv(
    os.path.join(OUT_DIR, "shuffled_spectral_metrics.csv"),
    index=False
)


# ============================================================
# 10. CONFRONTO
# ============================================================

final_trained = trained_metrics[
    trained_metrics["epoch"] == EPOCHS
]

comparison = pd.concat(
    [untrained_metrics, shuffled_metrics, final_trained],
    ignore_index=True
)

comparison.to_csv(
    os.path.join(OUT_DIR, "three_controls_comparison.csv"),
    index=False
)

summary = comparison.groupby(
    ["model", "layer"]
)[[
    "spectral_norm",
    "effective_rank",
    "pr_sigma",
    "stable_rank",
    "kappa_gram"
]].agg(["mean", "std"])

print("\n=== THREE-CONTROL SUMMARY ===")
print(summary)


# ============================================================
# 11. TRAINING CURVES
# ============================================================

history_all = pd.concat(
    [run["history"] for run in all_runs],
    ignore_index=True
)

history_all.to_csv(
    os.path.join(OUT_DIR, "training_history.csv"),
    index=False
)

plt.figure(figsize=(8, 5))

for seed in history_all["seed"].unique():
    sub = history_all[history_all["seed"] == seed]
    plt.plot(
        sub["epoch"],
        sub["val_accuracy"],
        marker="o",
        label=f"seed {seed}"
    )

plt.xlabel("Epoch")
plt.ylabel("Validation accuracy")
plt.title("Training performance")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(
    os.path.join(OUT_DIR, "training_accuracy.png"),
    dpi=160
)
plt.show()


# ============================================================
# 12. EVOLUZIONE DELLO SPETTRO
# ============================================================

for metric, ylabel, filename in [
    ("spectral_norm", "Spectral norm", "spectral_evolution"),
    ("effective_rank", "Effective rank", "effective_rank"),
    ("pr_sigma", "Participation ratio", "participation_ratio"),
]:
    for layer in ["fc1", "fc2", "fc3"]:
        plt.figure(figsize=(8, 5))

        sub = trained_metrics[
            trained_metrics["layer"] == layer
        ]

        for seed in sub["seed"].unique():
            s = sub[sub["seed"] == seed]
            plt.plot(
                s["epoch"],
                s[metric],
                marker="o",
                label=f"seed {seed}"
            )

        plt.xlabel("Epoch")
        plt.ylabel(ylabel)
        plt.title(f"{ylabel} — {layer}")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                OUT_DIR,
                f"{layer}_{filename}.png"
            ),
            dpi=160
        )
        plt.show()


# ============================================================
# 13. OPERATORI FINALI
# ============================================================

final_operators = {}

for run in all_runs:
    seed = run["seed"]
    state = run["checkpoints"][-1]["state_dict"]
    final_operators[seed] = get_weight_matrices(state)


# ============================================================
# 14. PRINCIPAL MODES
# ============================================================

def principal_modes(W, k=10):
    W = np.asarray(W, dtype=np.float64)

    # Per evitare di costruire WW^T esplicitamente:
    # gli autovalori sono sigma^2 e i left singular vectors sono U.
    U, s, _ = svd(W, full_matrices=False)

    k = min(k, len(s))

    return s[:k] ** 2, U[:, :k]


principal_rows = []

for seed, layers in final_operators.items():
    for layer, W in layers.items():
        eigenvalues, _ = principal_modes(W, k=10)

        for i, value in enumerate(eigenvalues):
            principal_rows.append({
                "seed": seed,
                "layer": layer,
                "mode": i + 1,
                "eigenvalue_gram": float(value),
                "normalized": (
                    float(value / eigenvalues[0])
                    if eigenvalues[0] > 0 else np.nan
                )
            })

principal_df = pd.DataFrame(principal_rows)

principal_df.to_csv(
    os.path.join(OUT_DIR, "principal_modes.csv"),
    index=False
)


# ============================================================
# 15. TEST S^2 - a S - b I
# ============================================================

def polynomial_operator_test(W, normalize=True):
    """
    Fit:
        S^2 = a S + b I

    IMPORTANT:
    Il fit di a,b è descrittivo.
    La compatibilità con Fibonacci viene testata separatamente
    imponendo a=b=1.

    Normalizziamo S per rho(S) quando richiesto, perché altrimenti
    a,b dipendono fortemente dalla scala dei pesi.
    """

    W = np.asarray(W, dtype=np.float64)
    S = W @ W.T

    scale = np.linalg.norm(S, ord=2)

    if normalize and scale > EPS:
        S = S / scale

    n = S.shape[0]
    I = np.eye(n)
    S2 = S @ S

    def residual(theta):
        a, b = theta
        return (S2 - a * S - b * I).ravel()

    result = least_squares(
        residual,
        x0=np.array([1.0, 1.0]),
        method="trf"
    )

    a, b = result.x
    R_fit = S2 - a * S - b * I
    R_phi = S2 - S - I

    denom = np.linalg.norm(S2, ord="fro") + EPS

    fit_rel = np.linalg.norm(R_fit, ord="fro") / denom
    phi_rel = np.linalg.norm(R_phi, ord="fro") / denom

    roots = np.roots([1.0, -a, -b])

    return {
        "scale_normalized": bool(normalize),
        "a_fit": float(a),
        "b_fit": float(b),
        "fit_relative_residual": float(fit_rel),
        "phi_fixed_relative_residual": float(phi_rel),
        "delta_a_from_1": float(a - 1.0),
        "delta_b_from_1": float(b - 1.0),
        "root_plus": complex(roots[0]),
        "root_minus": complex(roots[1]),
        "dominant_root_abs": float(
            np.max(np.abs(roots))
        ),
    }


poly_rows = []

for seed, layers in final_operators.items():
    for layer, W in layers.items():
        result = polynomial_operator_test(W, normalize=True)
        poly_rows.append({
            "seed": seed,
            "layer": layer,
            **result
        })

polynomial_df = pd.DataFrame(poly_rows)

polynomial_df.to_csv(
    os.path.join(OUT_DIR, "polynomial_operator_test.csv"),
    index=False
)

print("\n=== POLYNOMIAL OPERATOR TEST ===")
print(polynomial_df)


# ============================================================
# 16. CORRECT COARSE-GRAINING VIA TRUNCATED SVD
# ============================================================

def coarse_grain_svd(W, rank):
    """
    W: m x n

    SVD:
        W = U Sigma V^T

    Coarse-grained operator:
        W_r = U_r^T W V_r = Sigma_r

    Restituisce un operatore quadrato r x r.
    """

    W = np.asarray(W, dtype=np.float64)
    U, s, Vt = svd(W, full_matrices=False)

    r = min(int(rank), len(s))

    if r < 1:
        raise ValueError("Coarse-graining rank must be >= 1.")

    W_r = U[:, :r].T @ W @ Vt[:r, :].T

    return W_r, s[:r]


def rg_spectrum(W, steps=RG_STEPS, reduction=RG_REDUCTION):
    """
    Costruisce una sequenza di operatori a risoluzione decrescente.

    La scala viene scelta geometricamente, ma mai sotto rank=1.
    """

    W = np.asarray(W, dtype=np.float64)

    max_rank = min(W.shape)
    rows = []

    current = W.copy()

    for step in range(steps + 1):
        metrics = spectral_metrics(current)
        metrics.update({
            "rg_step": step,
            "rank_requested": min(
                max_rank,
                max(1, int(round(max_rank * reduction ** step)))
            ),
            "rows": current.shape[0],
            "cols": current.shape[1],
        })

        rows.append(metrics)

        if step == steps or min(current.shape) <= 1:
            break

        next_rank = max(
            1,
            int(round(min(current.shape) * reduction))
        )

        current, _ = coarse_grain_svd(
            current,
            next_rank
        )

    return pd.DataFrame(rows)


rg_tables = []

for seed, layers in final_operators.items():
    for layer, W in layers.items():
        rg = rg_spectrum(W)

        rg["seed"] = seed
        rg["layer"] = layer

        rg_tables.append(rg)

rg_df = pd.concat(rg_tables, ignore_index=True)

rg_df.to_csv(
    os.path.join(OUT_DIR, "rg_spectral_flow.csv"),
    index=False
)

print("\n=== RG SPECTRAL FLOW ===")
print(
    rg_df[
        [
            "seed", "layer", "rg_step",
            "rank_requested", "spectral_norm",
            "effective_rank", "pr_sigma"
        ]
    ]
)


# ============================================================
# 17. SCALING / RECURRENCE TEST
# ============================================================

def recurrence_fit(values):
    """
    Testa:
        s_{n+1} = a s_n + b s_{n-1}

    Non si assume a=b=1.
    """

    y = np.asarray(values, dtype=np.float64)

    if len(y) < 3:
        return {
            "a_fit": np.nan,
            "b_fit": np.nan,
            "rss": np.nan,
            "phi_fixed_rss": np.nan
        }

    X = np.column_stack([y[1:-1], y[:-2]])
    target = y[2:]

    coef, *_ = np.linalg.lstsq(X, target, rcond=None)
    a, b = coef

    pred = X @ coef
    rss = np.sum((target - pred) ** 2)

    phi_pred = X @ np.array([1.0, 1.0])
    phi_rss = np.sum((target - phi_pred) ** 2)

    return {
        "a_fit": float(a),
        "b_fit": float(b),
        "rss": float(rss),
        "phi_fixed_rss": float(phi_rss)
    }


rec_rows = []

for (seed, layer), sub in rg_df.groupby(["seed", "layer"]):
    sub = sub.sort_values("rg_step")

    # Con pochi RG step il test è necessariamente esplorativo.
    result = recurrence_fit(
        sub["effective_rank"].to_numpy()
    )

    rec_rows.append({
        "seed": seed,
        "layer": layer,
        **result
    })

recurrence_df = pd.DataFrame(rec_rows)

recurrence_df.to_csv(
    os.path.join(OUT_DIR, "rg_recurrence_test.csv"),
    index=False
)

print("\n=== RG RECURRENCE TEST ===")
print(recurrence_df)


# ============================================================
# 18. GOLDEN-RATIO TEST ESPLICITO
# ============================================================

golden_rows = []

for (seed, layer), sub in rg_df.groupby(["seed", "layer"]):
    sub = sub.sort_values("rg_step")

    values = sub["effective_rank"].to_numpy()

    if len(values) >= 2:
        ratios = values[1:] / np.maximum(values[:-1], EPS)

        golden_rows.append({
            "seed": seed,
            "layer": layer,
            "mean_ratio": float(np.mean(ratios)),
            "std_ratio": float(np.std(ratios)),
            "phi_inverse": PHI_INV,
            "abs_error_to_phi_inverse": float(
                abs(np.mean(ratios) - PHI_INV)
            )
        })
    else:
        golden_rows.append({
            "seed": seed,
            "layer": layer,
            "mean_ratio": np.nan,
            "std_ratio": np.nan,
            "phi_inverse": PHI_INV,
            "abs_error_to_phi_inverse": np.nan
        })

golden_df = pd.DataFrame(golden_rows)

golden_df.to_csv(
    os.path.join(OUT_DIR, "golden_ratio_test.csv"),
    index=False
)

print("\n=== GOLDEN RATIO TEST ===")
print(golden_df)


# ============================================================
# 19. FINAL AUDIT
# ============================================================

audit = {
    "n_seeds": N_SEEDS,
    "epochs": EPOCHS,
    "phi": PHI,
    "phi_inverse": PHI_INV,
    "polynomial_test_is_fixed_phi": False,
    "coarse_graining": "truncated SVD",
    "recurrence_is_descriptive": True,
    "golden_ratio_is_not_assumed": True,
}

with open(
    os.path.join(OUT_DIR, "phase1_audit.txt"),
    "w",
    encoding="utf-8"
) as f:
    for k, v in audit.items():
        f.write(f"{k}: {v}\n")

print("\n" + "=" * 78)
print("PHASE I COMPLETED")
print("=" * 78)
print("Results:", OUT_DIR)
