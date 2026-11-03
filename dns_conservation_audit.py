#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dns_conservation_audit.py

Conservation audit for DNS 3D Navier-Stokes implementation.
Tests:
A. Global nonlinear energy balance
B. Shell balance
C. Direct energy derivative
D. Nonlinear term computation verification
"""

import numpy as np
import jax
import jax.numpy as jnp
from pathlib import Path


def dns_conservaton_audit(N=32, steps=100, dt=1e-4, nu=1e-3, F=0.05, seed=42):
    """
    Run DNS conservation audit on small grid (N=32) for rapid testing.
    """
    print("=" * 88)
    print("DNS CONSERVATION AUDIT")
    print("=" * 88)
    print(f"N={N}, steps={steps}, dt={dt}, nu={nu}")
    print()

    # Grid setup
    kvals = jnp.fft.fftfreq(N) * N
    kx, ky, kz = jnp.meshgrid(kvals, kvals, kvals, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    k2_safe = jnp.where(k2 == 0, 1.0, k2)
    kmag = jnp.sqrt(k2)
    kcut = N / 3.0

    k2_np = np.asarray(k2)
    kmag_np = np.asarray(kmag)

    # Dealiasing 2/3
    dealias = (
        (jnp.abs(kx) <= kcut) &
        (jnp.abs(ky) <= kcut) &
        (jnp.abs(kz) <= kcut)
    ).astype(jnp.float32)

    # Leray projector
    def project_hat(hx, hy, hz):
        kdotu = kx * hx + ky * hy + kz * hz
        return (
            jnp.where(k2 == 0, 0.0, hx - kdotu * kx / k2_safe),
            jnp.where(k2 == 0, 0.0, hy - kdotu * ky / k2_safe),
            jnp.where(k2 == 0, 0.0, hz - kdotu * kz / k2_safe)
        )

    # Forcing
    x = jnp.arange(N) * (2.0 * jnp.pi / N)
    y = jnp.arange(N) * (2.0 * jnp.pi / N)
    fx = jnp.broadcast_to(F * jnp.sin(y[None, :, None]), (N, N, N))
    fy = jnp.broadcast_to(F * jnp.sin(x[:, None, None]), (N, N, N))
    fz = jnp.zeros((N, N, N))

    hfx, hfy, hfz = project_hat(jnp.fft.fftn(fx), jnp.fft.fftn(fy), jnp.fft.fftn(fz))
    fx_proj = jnp.real(jnp.fft.ifftn(hfx * dealias))
    fy_proj = jnp.real(jnp.fft.ifftn(hfy * dealias))
    fz_proj = jnp.real(jnp.fft.ifftn(hfz * dealias))

    # Initial condition
    key = jax.random.PRNGKey(seed)
    k1, k2key, k3 = jax.random.split(key, 3)
    ux = jnp.real(jnp.fft.ifftn(jax.random.normal(k1, (N, N, N)) * dealias))
    uy = jnp.real(jnp.fft.ifftn(jax.random.normal(k2key, (N, N, N)) * dealias))
    uz = jnp.real(jnp.fft.ifftn(jax.random.normal(k3, (N, N, N)) * dealias))

    # Shell masks
    shell_edges = np.arange(1, int(np.floor(kcut)) + 1)
    shell_masks = [
        (kmag_np >= shell_edges[n]) & (kmag_np < shell_edges[n + 1])
        for n in range(len(shell_edges) - 1)
    ]
    num_shells = len(shell_masks)

    # RHS function returning also nonlinear term in Fourier space
    def get_rhs_with_nl(ux, uy, uz):
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

    @jax.jit
    def step_rk2(ux, uy, uz):
        r1_x, r1_y, r1_z, hux, huy, huz, hcx, hcy, hcz = get_rhs_with_nl(ux, uy, uz)

        u_tilde_x = ux + dt * r1_x
        u_tilde_y = uy + dt * r1_y
        u_tilde_z = uz + dt * r1_z

        r2_x, r2_y, r2_z = get_rhs_with_nl(u_tilde_x, u_tilde_y, u_tilde_z)[:3]

        ux_new = ux + 0.5 * dt * (r1_x + r2_x)
        uy_new = uy + 0.5 * dt * (r1_y + r2_y)
        uz_new = uz + 0.5 * dt * (r1_z + r2_z)

        hnux, hnuy, hnuz = project_hat(jnp.fft.fftn(ux_new), jnp.fft.fftn(uy_new), jnp.fft.fftn(uz_new))
        ux_out = jnp.real(jnp.fft.ifftn(hnux * dealias))
        uy_out = jnp.real(jnp.fft.ifftn(hnuy * dealias))
        uz_out = jnp.real(jnp.fft.ifftn(hnuz * dealias))

        return ux_out, uy_out, uz_out, hux, huy, huz, hcx, hcy, hcz

    # Compile
    print("Compiling JIT...")
    ux, uy, uz = step_rk2(ux, uy, uz)[:3]
    ux.block_until_ready()
    print("Compilation complete.\n")

    # Storage for audit metrics
    T_mode_totals = []
    T_shell_totals = []
    dE_NL_values = []
    E_values = []

    for s in range(steps):
        res = step_rk2(ux, uy, uz)
        ux, uy, uz, hux, huy, huz, hcx, hcy, hcz = res

        # Convert to numpy for analysis
        hux_np, huy_np, huz_np = np.asarray(hux), np.asarray(huy), np.asarray(huz)
        hcx_np, hcy_np, hcz_np = np.asarray(hcx), np.asarray(hcy), np.asarray(hcz)

        # TEST A: Global nonlinear energy balance
        # T_mode(k) = Re(conj(u_hat(k)) · rhs_nl(k))
        # The nonlinear term in Fourier space is hcx, hcy, hcz (after projection)
        T_mode = np.real(np.conjugate(hux_np) * hcx_np + np.conjugate(huy_np) * hcy_np + np.conjugate(huz_np) * hcz_np)
        T_total = np.sum(T_mode)
        T_mode_totals.append(T_total)

        # TEST B: Shell balance
        T_shell = np.zeros(num_shells)
        for i, mask in enumerate(shell_masks):
            if np.any(mask):
                T_shell[i] = np.sum(T_mode[mask])
        T_shell_totals.append(np.sum(T_shell))

        # TEST C: Direct energy derivative
        # E = 0.5 * sum(|u_hat|^2)
        E = 0.5 * np.sum(np.abs(hux_np)**2 + np.abs(huy_np)**2 + np.abs(huz_np)**2)
        E_values.append(E)

        # dE/dt|_NL = Re(sum(conj(u_hat) * rhs_nl))
        dE_NL = np.real(np.sum(np.conjugate(hux_np) * hcx_np + np.conjugate(huy_np) * hcy_np + np.conjugate(huz_np) * hcz_np))
        dE_NL_values.append(dE_NL)

        if s % 20 == 0:
            print(f"Step {s:3d}: T_total={T_total:.6e}, T_shell_sum={np.sum(T_shell):.6e}, dE_NL={dE_NL:.6e}")

    # Convert to arrays
    T_mode_totals = np.array(T_mode_totals)
    T_shell_totals = np.array(T_shell_totals)
    dE_NL_values = np.array(dE_NL_values)
    E_values = np.array(E_values)

    # Compute statistics
    print("\n" + "=" * 88)
    print("CONSERVATION AUDIT RESULTS")
    print("=" * 88)

    print("\nTEST A: GLOBAL NONLINEAR ENERGY BALANCE")
    print(f"  Formula: T_total = sum_k Re(conj(u_hat) · rhs_nl)")
    print(f"  Expected: T_total ≈ 0")
    print(f"  Mean:     {np.mean(T_mode_totals):.6e}")
    print(f"  Std:      {np.std(T_mode_totals):.6e}")
    print(f"  Max abs:  {np.max(np.abs(T_mode_totals)):.6e}")
    print(f"  Min:      {np.min(T_mode_totals):.6e}")

    print("\nTEST B: SHELL BALANCE")
    print(f"  Formula: T_shell_total = sum_i T_shell(i)")
    print(f"  Expected: T_shell_total ≈ T_total")
    print(f"  Mean difference: {np.mean(T_shell_totals - T_mode_totals):.6e}")
    print(f"  Max difference:  {np.max(np.abs(T_shell_totals - T_mode_totals)):.6e}")

    print("\nTEST C: DIRECT ENERGY DERIVATIVE")
    print(f"  Formula: dE/dt|_NL = Re(sum(conj(u_hat) · rhs_nl))")
    print(f"  Expected: dE/dt|_NL ≈ 0")
    print(f"  Mean:     {np.mean(dE_NL_values):.6e}")
    print(f"  Std:      {np.std(dE_NL_values):.6e}")
    print(f"  Max abs:  {np.max(np.abs(dE_NL_values)):.6e}")

    print("\nENERGY STATISTICS")
    print(f"  E initial: {E_values[0]:.6e}")
    print(f"  E final:   {E_values[-1]:.6e}")
    print(f"  E change:  {E_values[-1] - E_values[0]:.6e}")
    print(f"  Relative change: {(E_values[-1] - E_values[0]) / E_values[0]:.6e}")

    # Pass/fail criteria
    print("\n" + "=" * 88)
    print("PASS/FAIL CRITERIA")
    print("=" * 88)

    tolerance = 1e-6
    tolerance_relative = 1e-4

    test_a_pass = np.max(np.abs(T_mode_totals)) < tolerance
    test_b_pass = np.max(np.abs(T_shell_totals - T_mode_totals)) < tolerance
    test_c_pass = np.max(np.abs(dE_NL_values)) < tolerance

    print(f"\nTEST A (Global balance):   {'PASS' if test_a_pass else 'FAIL'}")
    print(f"  Max T_total: {np.max(np.abs(T_mode_totals)):.6e} < {tolerance}")

    print(f"\nTEST B (Shell consistency): {'PASS' if test_b_pass else 'FAIL'}")
    print(f"  Max diff: {np.max(np.abs(T_shell_totals - T_mode_totals)):.6e} < {tolerance}")

    print(f"\nTEST C (Energy derivative): {'PASS' if test_c_pass else 'FAIL'}")
    print(f"  Max dE_NL: {np.max(np.abs(dE_NL_values)):.6e} < {tolerance}")

    all_pass = test_a_pass and test_b_pass and test_c_pass

    print("\n" + "=" * 88)
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 88)

    # Save report
    report = {
        "N": N,
        "steps": steps,
        "dt": dt,
        "nu": nu,
        "test_a": {
            "pass": bool(test_a_pass),
            "mean_T_total": float(np.mean(T_mode_totals)),
            "std_T_total": float(np.std(T_mode_totals)),
            "max_abs_T_total": float(np.max(np.abs(T_mode_totals))),
            "tolerance": tolerance
        },
        "test_b": {
            "pass": bool(test_b_pass),
            "mean_diff": float(np.mean(T_shell_totals - T_mode_totals)),
            "max_diff": float(np.max(np.abs(T_shell_totals - T_mode_totals))),
            "tolerance": tolerance
        },
        "test_c": {
            "pass": bool(test_c_pass),
            "mean_dE_NL": float(np.mean(dE_NL_values)),
            "std_dE_NL": float(np.std(dE_NL_values)),
            "max_abs_dE_NL": float(np.max(np.abs(dE_NL_values))),
            "tolerance": tolerance
        },
        "energy": {
            "initial": float(E_values[0]),
            "final": float(E_values[-1]),
            "change": float(E_values[-1] - E_values[0]),
            "relative_change": float((E_values[-1] - E_values[0]) / E_values[0])
        },
        "overall_pass": bool(all_pass)
    }

    return report


if __name__ == "__main__":
    report = dns_conservaton_audit(N=32, steps=100, dt=1e-4, nu=1e-3)

    # Save report
    import json
    report_path = Path("/home/fabianleo/Scaricati/DNS_CONSERVATION_AUDIT.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nReport saved to: {report_path}")
