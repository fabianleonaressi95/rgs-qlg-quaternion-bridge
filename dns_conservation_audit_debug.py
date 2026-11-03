#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dns_conservation_audit_debug.py

Debug version to identify shell balance issue.
"""

import numpy as np


def dns_audit_debug(N=32):
    """Debug shell decomposition."""
    print("=" * 88)
    print("DNS SHELL DECOMPOSITION DEBUG")
    print("=" * 88)

    # Grid setup
    kvals = np.fft.fftfreq(N) * N
    kx, ky, kz = np.meshgrid(kvals, kvals, kvals, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    kmag = np.sqrt(k2)
    kcut = N / 3.0

    # Shell masks
    shell_edges = np.arange(1, int(np.floor(kcut)) + 1)
    print(f"Shell edges: {shell_edges}")
    print(f"Number of shells: {len(shell_edges) - 1}")

    shell_masks = []
    for n in range(len(shell_edges) - 1):
        mask = (kmag >= shell_edges[n]) & (kmag < shell_edges[n + 1])
        shell_masks.append(mask)
        print(f"Shell {n}: k in [{shell_edges[n]}, {shell_edges[n+1]}), count={np.sum(mask)}")

    # Check if shells are disjoint and cover the dealiased region
    all_shells = np.zeros_like(kmag, dtype=bool)
    for mask in shell_masks:
        all_shells = all_shells | mask

    dealias_mask = (np.abs(kx) <= kcut) & (np.abs(ky) <= kcut) & (np.abs(kz) <= kcut)

    print(f"\nDealiased modes count: {np.sum(dealias_mask)}")
    print(f"Total in shells: {np.sum(all_shells)}")
    print(f"Overlap between shells and dealias: {np.sum(all_shells & dealias_mask)}")
    print(f"Modes in dealias but NOT in shells: {np.sum(dealias_mask & ~all_shells)}")

    # Check the k=0 mode
    k0_mask = (k2 == 0)
    print(f"\nk=0 mode count: {np.sum(k0_mask)}")
    print(f"k=0 in dealias: {np.sum(k0_mask & dealias_mask)}")
    print(f"k=0 in shells: {np.sum(k0_mask & all_shells)}")

    # Create a test field and check decomposition
    test_field = np.random.randn(N, N, N) + 1j * np.random.randn(N, N, N)
    test_field_dealiased = test_field * dealias_mask

    total_sum = np.sum(test_field_dealiased)
    shell_sum = 0.0
    for mask in shell_masks:
        shell_sum += np.sum(test_field_dealiased[mask])

    print(f"\nTest field decomposition:")
    print(f"Total sum (dealiased): {total_sum}")
    print(f"Sum over shells: {shell_sum}")
    print(f"Difference: {abs(total_sum - shell_sum)}")

    # The issue: shells don't include k=0 mode, but dealias does
    # Also shells may not cover all dealiased modes perfectly


if __name__ == "__main__":
    dns_audit_debug(N=32)
