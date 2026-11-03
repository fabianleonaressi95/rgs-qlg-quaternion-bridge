#!/usr/bin/env python3
"""
CUORE Pulse Shape Challenge — reproducible comparison
Models:
  1) Baseline: logistic regression on fixed pulse-shape summary features
  2) CNN1D: compact convolutional network (PyTorch)
  3) RGS-PRISM candidate: fixed multiscale temporal-graph summaries + logistic regression

Important:
- The RGS-PRISM model here is an explicit, testable candidate feature map, NOT a claim
  that the full theoretical RGS–PRISM operator has been derived or validated.
- Test labels are not used. The script creates test probabilities for class 1 (pile-up).
- All validation metrics are calculated only on a stratified split of the labelled training set.
"""

import argparse
import json
import random
import warnings
from pathlib import Path

import h5py
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
from sklearn.ensemble import HistGradientBoostingClassifier

SEED = 20261015


def set_seed(seed=SEED):
    random.seed(seed)
    np.random.seed(seed)


def find_key(h5, candidates):
    keys = list(h5.keys())
    # exact match first, then case-insensitive match
    for c in candidates:
        if c in h5:
            return c
    lower = {k.lower().replace("_", "").replace(" ", ""): k for k in keys}
    for c in candidates:
        norm = c.lower().replace("_", "").replace(" ", "")
        if norm in lower:
            return lower[norm]
    raise KeyError(f"Could not find any of {candidates}. Available HDF5 keys: {keys}")


def load_h5(path, need_labels):
    with h5py.File(path, "r") as f:
        waveform_key = find_key(f, ["Waveform", "Waveforms", "pulse", "pulses"])
        x = np.asarray(f[waveform_key], dtype=np.float32)
        if need_labels:
            label_key = find_key(f, ["Labels", "Label", "labels", "label"])
            y = np.asarray(f[label_key]).reshape(-1).astype(np.int64)
        else:
            y = None

        def optional(names):
            try:
                return np.asarray(f[find_key(f, names)], dtype=np.float32).reshape(-1)
            except KeyError:
                return None

        offset = optional(["normalizationOffset", "normalization_offset"])
        scale = optional(["normalizationScale", "normalization_scale"])
        maximum = optional(["normalizationMaximum", "normalization_maximum"])
        ids = optional(["id", "IDs", "event_id"])
    if x.ndim != 2:
        raise ValueError(f"Expected waveform matrix (N,T), got {x.shape}")
    if need_labels and len(y) != len(x):
        raise ValueError(f"Labels length {len(y)} != number of waveforms {len(x)}")
    return x, y, offset, scale, maximum, ids


def normalize_waveforms(x, offset, scale):
    # Follow the data-release definition where per-event parameters are available.
    x = np.asarray(x, dtype=np.float32)
    if offset is not None and scale is not None and len(offset) == len(x) and len(scale) == len(x):
        z = (x - offset[:, None]) / (scale[:, None] + 1e-8)
    else:
        # Fallback equivalent for this release: baseline mean over first 2500 samples,
        # scale = waveform maximum - baseline. Logged so this is auditable.
        warnings.warn("Normalization fields absent: calculating offset/scale from waveforms.")
        mu = np.mean(x[:, :2500], axis=1)
        s = np.max(x, axis=1) - mu
        z = (x - mu[:, None]) / (s[:, None] + 1e-8)
    return np.clip(z, 0.0, 1.0).astype(np.float32)


def downsample_mean(z, n_bins):
    n, t = z.shape
    usable = (t // n_bins) * n_bins
    return z[:, :usable].reshape(n, n_bins, usable // n_bins).mean(axis=2)


def physical_features(z):
    """Fixed, interpretable pulse-shape summaries; no learned preprocessing."""
    n, t = z.shape
    idx = np.arange(t, dtype=np.float32)
    feats = []
    # Global waveform shape and coarse time occupancy.
    feats.append(downsample_mean(z, 64))
    feats.append(downsample_mean(z, 256))
    # Derivative summaries.
    d1 = np.diff(z, axis=1)
    d2 = np.diff(z, n=2, axis=1)
    for d in (d1, d2):
        feats.append(np.stack([
            np.mean(d, axis=1), np.std(d, axis=1),
            np.mean(np.abs(d), axis=1), np.max(d, axis=1),
            np.min(d, axis=1), np.quantile(np.abs(d), .95, axis=1)
        ], axis=1))
    # Robust amplitude/shape statistics.
    feats.append(np.stack([
        np.max(z, axis=1), np.mean(z, axis=1), np.std(z, axis=1),
        np.quantile(z, .90, axis=1), np.quantile(z, .99, axis=1),
        np.sum(z > .2, axis=1) / t, np.sum(z > .5, axis=1) / t,
        np.sum(z > .8, axis=1) / t,
        np.mean(z[:, :2500], axis=1), np.std(z[:, :2500], axis=1),
        np.mean(z[:, 2500:5000], axis=1), np.mean(z[:, 5000:7500], axis=1),
        np.mean(z[:, 7500:], axis=1)
    ], axis=1))
    return np.concatenate(feats, axis=1).astype(np.float32)


def rgs_prism_candidate_features(z):
    """
    Transparent candidate geometric map:
    fixed temporal scales (16, 32, 64, 128, 256 bins), first/second differences,
    and adjacent-bin contrasts. This operationalizes multiscale temporal geometry
    without claiming to implement the entire RGS–PRISM theory.
    """
    blocks = []
    for bins in (16, 32, 64, 128, 256):
        pooled = downsample_mean(z, bins)
        blocks.append(pooled)
        # Local edge/transition structure on the same fixed grid.
        blocks.append(np.diff(pooled, axis=1))
        if pooled.shape[1] >= 3:
            blocks.append(np.diff(pooled, n=2, axis=1))
        # Window-to-window contrast (a simple discrete graph-gradient summary).
        blocks.append(np.abs(np.diff(pooled, axis=1)))
    # Add full-resolution coarse physical summaries to preserve interpretability.
    blocks.append(physical_features(z))
    return np.concatenate(blocks, axis=1).astype(np.float32)


def fpr_at_tpr(y_true, scores, target_tpr=0.90):
    fpr, tpr, _ = roc_curve(y_true, scores)
    eligible = np.flatnonzero(tpr >= target_tpr)
    return float(fpr[eligible[0]]) if len(eligible) else float("nan")


def report_metrics(name, y, p, threshold=0.5):
    out = {
        "model": name,
        "auc_roc": float(roc_auc_score(y, p)),
        "accuracy_at_0.5": float(accuracy_score(y, p >= threshold)),
        "fpr_at_tpr_0.90": fpr_at_tpr(y, p, 0.90),
    }
    print(json.dumps(out, indent=2))
    return out


def fit_cnn(x_train, y_train, x_val, y_val, x_test, epochs=12, batch_size=64):
    try:
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as e:
        raise RuntimeError("PyTorch is not installed. Install it to run CNN1D, or use --skip-cnn.") from e

    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(SEED)
    # Downsample deterministically to 2500 points to reduce compute while retaining
    # the 10 s event window and millisecond-scale structure at 4 ms per point.
    def prep(a):
        n, t = a.shape
        usable = (t // 2500) * 2500
        return a[:, :usable].reshape(n, 2500, usable // 2500).mean(axis=2).astype(np.float32)

    xt, xv, xe = prep(x_train), prep(x_val), prep(x_test)
    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.features = nn.Sequential(
                nn.Conv1d(1, 16, 9, padding=4), nn.BatchNorm1d(16), nn.ReLU(),
                nn.MaxPool1d(4),
                nn.Conv1d(16, 32, 9, padding=4), nn.BatchNorm1d(32), nn.ReLU(),
                nn.MaxPool1d(4),
                nn.Conv1d(32, 64, 7, padding=3), nn.BatchNorm1d(64), nn.ReLU(),
                nn.AdaptiveAvgPool1d(32),
            )
            self.head = nn.Sequential(nn.Flatten(), nn.Linear(64*32, 64), nn.ReLU(),
                                      nn.Dropout(0.25), nn.Linear(64, 1))
        def forward(self, x):
            return self.head(self.features(x))

    def tensor(a):
        return torch.from_numpy(a[:, None, :]).float()
    train_ds = TensorDataset(tensor(xt), torch.from_numpy(y_train.astype(np.float32)))
    loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    model = Net().to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()
    xv_t = tensor(xv).to(device)
    best_auc, best_state, patience = -1.0, None, 3
    for epoch in range(epochs):
        model.train()
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            opt.zero_grad()
            loss = loss_fn(model(xb).squeeze(1), yb)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            pv = torch.sigmoid(model(xv_t).squeeze(1)).cpu().numpy()
        auc = roc_auc_score(y_val, pv)
        print(f"CNN1D epoch {epoch+1}/{epochs}: validation AUC={auc:.6f}")
        if auc > best_auc:
            best_auc, best_state, patience = auc, {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}, 3
        else:
            patience -= 1
            if patience == 0:
                break
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        pv = torch.sigmoid(model(xv_t).squeeze(1)).cpu().numpy()
        pt = torch.sigmoid(model(tensor(xe).to(device)).squeeze(1)).cpu().numpy()
    return pv, pt, model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, help="Path to cuoreTraining.h5")
    parser.add_argument("--test", required=True, help="Path to cuoreTest.h5")
    parser.add_argument("--outdir", default="cuore_challenge_output")
    parser.add_argument("--epochs", type=int, default=12)
    parser.add_argument("--skip-cnn", action="store_true")
    args = parser.parse_args()
    set_seed()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    x_raw, y, off, scale, maximum, train_ids = load_h5(args.train, True)
    xt_raw, _, toff, tscale, tmax, test_ids = load_h5(args.test, False)
    print(f"Train raw shape={x_raw.shape}; test raw shape={xt_raw.shape}")
    print(f"Train labels/counts={np.unique(y, return_counts=True)}")
    if x_raw.shape[1] != 10000 or xt_raw.shape[1] != 10000:
        warnings.warn("Expected 10000 samples per waveform according to the release paper.")
    if not set(np.unique(y)).issubset({0, 1}):
        raise ValueError(f"Unexpected labels: {np.unique(y)}")
    if len(xt_raw) != 1000:
        warnings.warn(f"Expected 1000 test events from the paper, found {len(xt_raw)}.")

    x = normalize_waveforms(x_raw, off, scale)
    x_test = normalize_waveforms(xt_raw, toff, tscale)
    idx_tr, idx_val = train_test_split(
        np.arange(len(y)), test_size=0.20, random_state=SEED, stratify=y
    )
    y_tr, y_val = y[idx_tr], y[idx_val]
    results = []

    # Baseline: physically interpretable fixed summaries.
    f = physical_features(x)
    baseline = make_pipeline(
        SimpleImputer(strategy="median"), StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
    )
    baseline.fit(f[idx_tr], y_tr)
    p_val = baseline.predict_proba(f[idx_val])[:, 1]
    results.append(report_metrics("baseline_physical_logreg", y_val, p_val))
    p_test_base = baseline.predict_proba(physical_features(x_test))[:, 1]

    # RGS–PRISM candidate: fixed multiscale temporal geometry, no hand-tuned coefficients.
    fg = rgs_prism_candidate_features(x)
    rgs = make_pipeline(
        SimpleImputer(strategy="median"), StandardScaler(),
        LogisticRegression(max_iter=2500, class_weight="balanced", C=0.1, random_state=SEED)
    )
    rgs.fit(fg[idx_tr], y_tr)
    pg_val = rgs.predict_proba(fg[idx_val])[:, 1]
    results.append(report_metrics("rgs_prism_multiscale_candidate", y_val, pg_val))
    pg_test = rgs.predict_proba(rgs_prism_candidate_features(x_test))[:, 1]

    # Optional non-linear tabular reference model on the geometric features.
    hgb = make_pipeline(
        SimpleImputer(strategy="median"),
        HistGradientBoostingClassifier(max_iter=150, learning_rate=0.06,
                                       l2_regularization=1.0, random_state=SEED)
    )
    hgb.fit(fg[idx_tr], y_tr)
    ph_val = hgb.predict_proba(fg[idx_val])[:, 1]
    results.append(report_metrics("multiscale_hist_gradient_boosting", y_val, ph_val))
    ph_test = hgb.predict_proba(rgs_prism_candidate_features(x_test))[:, 1]

    predictions = {
        "baseline_physical_logreg": p_test_base,
        "rgs_prism_multiscale_candidate": pg_test,
        "multiscale_hist_gradient_boosting": ph_test,
    }

    if not args.skip_cnn:
        try:
            pc_val, pc_test, cnn = fit_cnn(
                x[idx_tr], y_tr, x[idx_val], y_val, x_test, epochs=args.epochs
            )
            results.append(report_metrics("cnn1d", y_val, pc_val))
            predictions["cnn1d"] = pc_test
            import torch
            torch.save(cnn.state_dict(), outdir / "cnn1d_state_dict.pt")
        except RuntimeError as e:
            warnings.warn(str(e))
            print("CNN1D skipped; other models still completed.")

    # Save validation metrics and all test predictions, one probability per event.
    with open(outdir / "validation_metrics.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    for name, prob in predictions.items():
        if len(prob) != len(xt_raw) or not np.all(np.isfinite(prob)):
            raise ValueError(f"Invalid predictions for {name}")
        arr = np.column_stack([test_ids, prob]) if test_ids is not None else prob.reshape(-1, 1)
        np.savetxt(outdir / f"predictions_{name}.csv", arr, delimiter=",",
                   header="id,probability_class1" if test_ids is not None else "probability_class1",
                   comments="", fmt=["%d", "%.9g"] if test_ids is not None else "%.9g")

    # Primary submission candidate is selected only after validation AUC comparison.
    best = max(results, key=lambda r: r["auc_roc"])
    print("\nValidation-only best model:", best["model"], "AUC=", best["auc_roc"])
    print("Predictions saved. IMPORTANT: inspect the official submission requirements before emailing.")
    print("Output directory:", outdir.resolve())


if __name__ == "__main__":
    main()
