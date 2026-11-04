#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dns_original_extracted.py

DNS 3D Navier-Stokes JAX implementation extracted from Untitled2.ipynb
Source: /home/fabianleo/Scaricati/Untitled2.ipynb
Status: ORIGINAL KAGGLE CODE
"""

import time
import numpy as np
import jax
import jax.numpy as jnp


def run_dns_rk2_filtered(
    N=128,
    steps=25000,
    dt=1e-4,
    nu=2e-3,       # Viscosità ridotta per estendere l'intervallo inerziale
    F=0.05,
    sample_every=10,
    transient=15000,
    seed=42
):
    """
    Original DNS implementation from Kaggle notebook.
    3D periodic box, pseudo-spectral Fourier, N=128.
    """
    t0 = time.time()

    # 1. Griglia e operatori base
    kvals = jnp.fft.fftfreq(N) * N
    kx, ky, kz = jnp.meshgrid(kvals, kvals, kvals, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    k2_safe = jnp.where(k2 == 0, 1.0, k2)
    kmag = jnp.sqrt(k2)
    kcut = N / 3.0

    k2_np = np.asarray(k2)
    kmag_np = np.asarray(kmag)

    # 2. Dealiasing 2/3
    dealias = (
        (jnp.abs(kx) <= kcut) &
        (jnp.abs(ky) <= kcut) &
        (jnp.abs(kz) <= kcut)
    ).astype(jnp.float32)

    # 3. Proiettore di Leray
    def project_hat(hx, hy, hz):
        kdotu = kx * hx + ky * hy + kz * hz
        return (
            jnp.where(k2 == 0, 0.0, hx - kdotu * kx / k2_safe),
            jnp.where(k2 == 0, 0.0, hy - kdotu * ky / k2_safe),
            jnp.where(k2 == 0, 0.0, hz - kdotu * kz / k2_safe)
        )

    # 4. Forcing
    x = jnp.arange(N) * (2.0 * jnp.pi / N)
    y = jnp.arange(N) * (2.0 * jnp.pi / N)
    fx = jnp.broadcast_to(F * jnp.sin(y[None, :, None]), (N, N, N))
    fy = jnp.broadcast_to(F * jnp.sin(x[:, None, None]), (N, N, N))
    fz = jnp.zeros((N, N, N))

    hfx, hfy, hfz = project_hat(jnp.fft.fftn(fx), jnp.fft.fftn(fy), jnp.fft.fftn(fz))
    fx_proj = jnp.real(jnp.fft.ifftn(hfx * dealias))
    fy_proj = jnp.real(jnp.fft.ifftn(hfy * dealias))
    fz_proj = jnp.real(jnp.fft.ifftn(hfz * dealias))

    # 5. Condizione iniziale
    key = jax.random.PRNGKey(seed)
    k1, k2key, k3 = jax.random.split(key, 3)
    ux = jnp.real(jnp.fft.ifftn(jax.random.normal(k1, (N, N, N)) * dealias))
    uy = jnp.real(jnp.fft.ifftn(jax.random.normal(k2key, (N, N, N)) * dealias))
    uz = jnp.real(jnp.fft.ifftn(jax.random.normal(k3, (N, N, N)) * dealias))

    # 6. Shells spettrali
    shell_edges = np.arange(1, int(np.floor(kcut)) + 1)
    shell_masks = [
        (kmag_np >= shell_edges[n]) & (kmag_np < shell_edges[n + 1])
        for n in range(len(shell_edges) - 1)
    ]
    num_shells = len(shell_masks)

    # Funzione RHS per RK2
    def get_rhs(ux, uy, uz):
        hux, huy, huz = jnp.fft.fftn(ux), jnp.fft.fftn(uy), jnp.fft.fftn(uz)

        dux = [jnp.real(jnp.fft.ifftn(1j * kx * hux)), jnp.real(jnp.fft.ifftn(1j * ky * hux)), jnp.real(jnp.fft.ifftn(1j * kz * hux))]
        duy = [jnp.real(jnp.fft.ifftn(1j * kx * huy)), jnp.real(jnp.fft.ifftn(1j * ky * huy)), jnp.real(jnp.fft.ifftn(1j * kz * huy))]
        duz = [jnp.real(jnp.fft.ifftn(1j * kx * huz)), jnp.real(jnp.fft.ifftn(1j * ky * huz)), jnp.real(jnp.fft.ifftn(1j * kz * huz))]

        cx = ux * dux[0] + uy * dux[1] + uz * dux[2]
        cy = ux * duy[0] + uy * duy[1] + uz * duy[2]
        cz = ux * duz[0] + uy * duz[1] + uz * duz[2]

        hcx, hcy, hcz = project_hat(jnp.fft.fftn(cx), jnp.fft.fftn(cy), jnp.fft.fftn(cz))
        hcx, hcy, hcz = hcx * dealias, hcy * dealias, hcz * dealias

        lapx = jnp.real(jnp.fft.ifftn(-k2 * hux))
        lapy = jnp.real(jnp.fft.ifftn(-k2 * huy))
        lapz = jnp.real(jnp.fft.ifftn(-k2 * huz))

        rhs_x = -jnp.real(jnp.fft.ifftn(hcx)) + nu * lapx + fx_proj
        rhs_y = -jnp.real(jnp.fft.ifftn(hcy)) + nu * lapy + fy_proj
        rhs_z = -jnp.real(jnp.fft.ifftn(hcz)) + nu * lapz + fz_proj

        return rhs_x, rhs_y, rhs_z, hux, huy, huz, hcx, hcy, hcz

    # 7. Schema RK2 JIT
    @jax.jit
    def step_rk2(ux, uy, uz):
        r1_x, r1_y, r1_z, hux, huy, huz, hcx, hcy, hcz = get_rhs(ux, uy, uz)

        # Predictor
        u_tilde_x = ux + dt * r1_x
        u_tilde_y = uy + dt * r1_y
        u_tilde_z = uz + dt * r1_z

        r2_x, r2_y, r2_z = get_rhs(u_tilde_x, u_tilde_y, u_tilde_z)[:3]

        # Corrector
        ux_new = ux + 0.5 * dt * (r1_x + r2_x)
        uy_new = uy + 0.5 * dt * (r1_y + r2_y)
        uz_new = uz + 0.5 * dt * (r1_z + r2_z)

        # Proiezione finale
        hnux, hnuy, hnuz = project_hat(jnp.fft.fftn(ux_new), jnp.fft.fftn(uy_new), jnp.fft.fftn(uz_new))
        ux_out = jnp.real(jnp.fft.ifftn(hnux * dealias))
        uy_out = jnp.real(jnp.fft.ifftn(hnuy * dealias))
        uz_out = jnp.real(jnp.fft.ifftn(hnuz * dealias))

        return ux_out, uy_out, uz_out, hux, huy, huz, hcx, hcy, hcz

    print("Compilazione JIT RK2...")
    ux, uy, uz = step_rk2(ux, uy, uz)[:3]
    ux.block_until_ready()
    print("Compilazione completata.")

    T_history = []

    for s in range(steps):
        res = step_rk2(ux, uy, uz)
        ux, uy, uz, hux, huy, huz, hcx, hcy, hcz = res

        if s >= transient and s % sample_every == 0:
            ux_h, uy_h, uz_h = np.asarray(hux), np.asarray(huy), np.asarray(huz)
            cx_h, cy_h, cz_h = np.asarray(hcx), np.asarray(hcy), np.asarray(hcz)

            T = np.zeros(num_shells)
            norm = float(N**6)

            for i, mask in enumerate(shell_masks):
                if np.any(mask):
                    transfer = (
                        np.real(np.conjugate(ux_h[mask]) * cx_h[mask]) +
                        np.real(np.conjugate(uy_h[mask]) * cy_h[mask]) +
                        np.real(np.conjugate(uz_h[mask]) * cz_h[mask])
                    )
                    T[i] = -np.sum(transfer) / norm

            T_history.append(T)

        if s % 2500 == 0:
            print(f"Step {s:6d} / {steps}")

    T_history = np.asarray(T_history)

    print(f"\nTempo di esecuzione: {time.time() - t0:.2f} s")

    return {
        "T_history": T_history,
        "shell_edges": shell_edges,
        "shell_masks": shell_masks,
        "k2_np": k2_np,
        "kmag_np": kmag_np,
        "dealias": dealias,
        "final_state": (ux, uy, uz)
    }


if __name__ == "__main__":
    result = run_dns_rk2_filtered()
    print(f"T_history shape: {result['T_history'].shape}")
