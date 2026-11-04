#!/usr/bin/env python3
"""
osm_audit.py -- descrittori locali dalla densita' rho(r) su griglia periodica ortorombica.

Quantita' (unita' atomiche di default, hbar = m = 1):
    p_osm = +(hbar/2) grad ln rho          (momento osmotico, convenzione Nelson: u = +hbar/(2m) grad ln rho)
    Q     = -(hbar^2/2m) lap(sqrt rho)/sqrt rho   (potenziale "quantistico" dalla densita')
    t_W   = (hbar^2/8m) |grad rho|^2 / rho        (densita' di energia cinetica di von Weizsacker)

ATTENZIONE (interpretazione):
  * Per molti elettroni sqrt(rho) NON e' la funzione d'onda: Q e p_osm sono descrittori di
    disomogeneita' della densita', non il potenziale quantistico fisico del sistema.
  * p_phase = grad S NON si ricava da rho. Non e' calcolato qui.
  * Vicino ai nuclei grad ln rho -> -2Z (a.u.) (cuspide): usare maschere sui nuclei.
  * Dove rho -> 0 ln rho diverge: usare soglia rho_floor.

Uso:
    python3 osm_audit.py --selftest
    python3 osm_audit.py --cube density.cube --core-radius 1.0 --outdir out
    python3 osm_audit.py --cube A.cube --compare B.cube        # confronto KS tra due strutture
"""
import argparse
import json
import os
import sys

import numpy as np

HBAR = 1.0
MASS = 1.0


# ----------------------------------------------------------------------------- derivate
def _d1(f, axis, h):
    return (np.roll(f, -1, axis) - np.roll(f, 1, axis)) / (2.0 * h)


def _d2(f, axis, h):
    return (np.roll(f, -1, axis) - 2.0 * f + np.roll(f, 1, axis)) / (h * h)


def gradient(f, h):
    return np.stack([_d1(f, i, h[i]) for i in range(3)], axis=-1)


def laplacian(f, h):
    return sum(_d2(f, i, h[i]) for i in range(3))


# ----------------------------------------------------------------------------- descrittori
def descriptors(rho, h, rho_floor=1e-8, hbar=HBAR, m=MASS):
    """rho: array 3D (periodico). h: passi (hx,hy,hz). Restituisce dict di array; NaN dove rho<rho_floor."""
    rho = np.asarray(rho, dtype=float)
    valid = rho > rho_floor
    r = np.where(valid, rho, rho_floor)
    L = np.log(r)
    gL = gradient(L, h)
    sq = np.sqrt(r)
    p_osm = 0.5 * hbar * gL
    Q = -(hbar**2 / (2 * m)) * laplacian(sq, h) / sq
    gr = gradient(r, h)
    t_w = (hbar**2 / (8 * m)) * np.sum(gr**2, axis=-1) / r
    out = {
        "p_osm": p_osm,
        "p_osm_norm": np.linalg.norm(p_osm, axis=-1),
        "Q": Q,
        "t_W": t_w,
        "lap_lnrho": laplacian(L, h),
    }
    for k, v in out.items():
        if v.ndim == 3:
            out[k] = np.where(valid, v, np.nan)
        else:
            out[k] = np.where(valid[..., None], v, np.nan)
    out["valid"] = valid
    return out


def identity_residual(d, h, hbar=HBAR, m=MASS):
    """Verifica Q = -(hbar^2/4m) lap(ln rho) - (hbar^2/8m)|grad ln rho|^2 (identita' esatta)."""
    rhs = -(hbar**2 / (4 * m)) * d["lap_lnrho"] - (hbar**2 / (8 * m)) * (
        (2 * d["p_osm"] / hbar) ** 2
    ).sum(-1)
    diff = np.abs(d["Q"] - rhs)
    return np.nanmax(diff), np.nanmedian(diff)


# ----------------------------------------------------------------------------- maschere e statistiche
def grid_coords(shape, h):
    return [np.arange(n) * h[i] for i, n in enumerate(shape)]


def core_mask(shape, h, cell, atoms, radius):
    """True fuori da sfere di raggio `radius` attorno agli atomi (minimum image, cella ortorombica)."""
    xs, ys, zs = grid_coords(shape, h)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    keep = np.ones(shape, bool)
    for pos in atoms:
        d2 = 0
        for C, p, L in zip((X, Y, Z), pos, cell):
            dd = np.abs(C - p)
            dd = np.minimum(dd, L - dd)
            d2 = d2 + dd**2
        keep &= d2 > radius**2
    return keep


def summarize(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0}
    q = np.quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95])
    return {"n": int(x.size), "mean": float(x.mean()), "std": float(x.std()),
            "q05": float(q[0]), "q25": float(q[1]), "median": float(q[2]),
            "q75": float(q[3]), "q95": float(q[4])}


def ks_2samp(a, b, max_n=200000, seed=0):
    """Statistica KS a due campioni + p-value asintotico. Sottocampiona: i voxel sono correlati,
    quindi il p-value e' solo indicativo; confrontare sempre con un controllo nullo."""
    rng = np.random.default_rng(seed)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size > max_n:
        a = rng.choice(a, max_n, replace=False)
    if b.size > max_n:
        b = rng.choice(b, max_n, replace=False)
    a = np.sort(a)
    b = np.sort(b)
    allv = np.concatenate([a, b])
    cdfa = np.searchsorted(a, allv, side="right") / a.size
    cdfb = np.searchsorted(b, allv, side="right") / b.size
    D = float(np.max(np.abs(cdfa - cdfb)))
    n = a.size * b.size / (a.size + b.size)
    lam = (np.sqrt(n) + 0.12 + 0.11 / np.sqrt(n)) * D
    p = 2 * sum((-1) ** (k - 1) * np.exp(-2 * (k * lam) ** 2) for k in range(1, 101))
    return D, float(min(max(p, 0.0), 1.0))


# ----------------------------------------------------------------------------- I/O cube
def read_cube(path):
    """Gaussian cube (VASP/QE/pp.x). Assume vettori voxel ortogonali. Unita' Bohr. Ritorna rho, h, cell, atoms."""
    with open(path) as f:
        lines = f.readlines()
    nat = int(lines[2].split()[0])
    vox = []
    n = []
    for i in range(3):
        t = lines[3 + i].split()
        n.append(int(t[0]))
        vox.append(np.array(t[1:4], float))
    vox = np.array(vox)
    off = vox - np.diag(np.diag(vox))
    if np.abs(off).max() > 1e-8:
        raise ValueError("Cella non ortorombica: non supportata da questo modulo.")
    if any(x < 0 for x in n):
        raise ValueError("Cube in Angstrom (n<0) non gestito: convertire in Bohr.")
    h = np.diag(vox)
    atoms = []
    for i in range(nat):
        t = lines[6 + i].split()
        atoms.append(np.array(t[2:5], float))
    data = np.array(" ".join(lines[6 + nat:]).split(), float)
    rho = data.reshape(n)  # ordine cube: x esterno, z interno
    cell = h * np.array(n)
    return rho, h, cell, atoms


# ----------------------------------------------------------------------------- self-test
def selftest():
    ok = True
    # 1) H 1s: rho = exp(-2r)/pi ; u = +grad ln rho /2 -> p_osm = -r_hat ; Q = -1/2 + 1/r
    L, N = 24.0, 240
    h = np.array([L / N] * 3)
    x = np.arange(N) * h[0] - L / 2
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.exp(-2 * R) / np.pi
    d = descriptors(rho, h, rho_floor=1e-10)
    region = (R > 1.0) & (R < 5.0)
    err_p = np.nanmedian(np.abs(d["p_osm_norm"][region] - 1.0))
    err_Q = np.nanmedian(np.abs(d["Q"][region] - (-0.5 + 1.0 / R[region])))
    print(f"[H 1s] |p_osm|=1 atteso: errore mediano {err_p:.2e}")
    print(f"[H 1s] Q=-1/2+1/r atteso: errore mediano {err_Q:.2e}")
    ok &= err_p < 5e-3 and err_Q < 5e-3
    # von Weizsacker per H 1s: t_W = rho/2 (|grad ln rho|^2/8 = 4/8 = 1/2)
    err_t = np.nanmedian(np.abs(d["t_W"][region] - 0.5 * rho[region]) / rho[region])
    print(f"[H 1s] t_W = rho/2: errore relativo mediano {err_t:.2e}")
    ok &= err_t < 5e-3
    # 2) identita' esatta
    mx, md = identity_residual(d, h)
    print(f"[identita' Q] max {mx:.2e}  mediana {md:.2e}  (discretizzazione, non deve essere ~1)")
    ok &= md < 1e-3
    # 3) Gaussiana 3D: rho ~ exp(-a r^2): grad ln rho = -2 a r ; Q = -(3a) + a^2 r^2 * ...  (verifica via identita' analitica)
    a = 0.5
    rho2 = np.exp(-a * R**2)
    d2 = descriptors(rho2, h, rho_floor=1e-12)
    reg2 = R < 4.0
    Qan = -0.5 * (a**2 * R**2 - 3 * a)
    err_g = np.nanmedian(np.abs(d2["Q"][reg2] - Qan[reg2]))
    print(f"[Gaussiana] Q analitico: errore mediano {err_g:.2e}")
    ok &= err_g < 5e-3
    print("SELFTEST:", "OK" if ok else "FALLITO")
    return 0 if ok else 1


# ----------------------------------------------------------------------------- CLI
def analyse_cube(path, core_radius, rho_floor, tag):
    rho, h, cell, atoms = read_cube(path)
    d = descriptors(rho, h, rho_floor=rho_floor)
    mask = d["valid"] & core_mask(rho.shape, h, cell, atoms, core_radius)
    res = {"file": path, "grid": list(rho.shape), "cell_bohr": cell.tolist(),
           "core_radius_bohr": core_radius, "frac_voxels_used": float(mask.mean()),
           "p_osm_norm": summarize(d["p_osm_norm"][mask]),
           "Q": summarize(d["Q"][mask]), "t_W": summarize(d["t_W"][mask])}
    # anisotropia: rapporto autovalori del tensore <p_i p_j>
    P = d["p_osm"][mask]
    T = (P[:, :, None] * P[:, None, :]).mean(0)
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    res["p_osm_tensor_eigs"] = ev.tolist()
    res["p_osm_anisotropy_ev1_over_ev3"] = float(ev[0] / ev[2]) if ev[2] > 0 else None
    return d, mask, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--cube")
    ap.add_argument("--compare")
    ap.add_argument("--core-radius", type=float, default=1.0, help="Bohr, esclude sfere attorno agli atomi")
    ap.add_argument("--rho-floor", type=float, default=1e-6)
    ap.add_argument("--outdir", default="osm_out")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.cube:
        ap.error("serve --cube oppure --selftest")
    os.makedirs(a.outdir, exist_ok=True)
    dA, mA, rA = analyse_cube(a.cube, a.core_radius, a.rho_floor, "A")
    out = {"A": rA}
    if a.compare:
        dB, mB, rB = analyse_cube(a.compare, a.core_radius, a.rho_floor, "B")
        out["B"] = rB
        D, p = ks_2samp(dA["p_osm_norm"][mA], dB["p_osm_norm"][mB])
        out["KS_p_osm_norm"] = {"D": D, "p_indicativo": p,
                                "nota": "voxel correlati: p-value ottimistico; confrontare con controllo nullo"}
    with open(os.path.join(a.outdir, "osm_report.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
