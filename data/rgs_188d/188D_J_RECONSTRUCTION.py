#!/usr/bin/env python3
"""188D J_A / J_B first-principles reconstruction gate.

Conclusion from verified local sources:
- J_A and J_B do not have a recoverable original mathematical definition.
- The only explicit formulas found are synthetic mock perturbations and random fallbacks.
- Therefore this script refuses to invent a formula.
"""

from pathlib import Path
import numpy as np

BASE = Path('/home/fabianleo/Scaricati')
LOCATIONS = [
    BASE / '188D_LOCAL_REPLAY_OUTPUT' / 'X_Oh_reconstructed.npy',
    BASE / 'X_Oh_reconstructed.npy',
]
X_Oh_path = next((p for p in LOCATIONS if p.exists()), None)
Q_path = next((p for p in [BASE / '188D_LOCAL_REPLAY_OUTPUT' / 'Q_rebuilt.npy', BASE / 'Q_rebuilt.npy'] if p.exists()), None)
K_path = next((p for p in [BASE / '188D_LOCAL_REPLAY_OUTPUT' / 'K_tangent_rebuilt.npy', BASE / 'K_tangent_rebuilt.npy'] if p.exists()), None)

if X_Oh_path is None or Q_path is None or K_path is None:
    print('Verified inputs not all present; verdict remains NOT_RECOVERABLE.')
    print('  X_Oh: ', X_Oh_path)
    print('  Q: ', Q_path)
    print('  K_tangent: ', K_path)
    raise SystemExit(0)

X_Oh = np.load(X_Oh_path)
Q = np.load(Q_path)
K_tangent = np.load(K_path)

print('Verified inputs:')
print('  X_Oh shape =', X_Oh.shape)
print('  Q shape =', Q.shape)
print('  K_tangent shape =', K_tangent.shape)

print("\nEvidence:")
print("  The recovered notebook contains synthetic mock formulas:")
print("    J_A = J + 0.01 * rng.normal(size=J.shape)")
print("    J_B = J + 0.01 * rng.normal(size=J.shape)")
print("  and a fallback random generation:")
print("    J_A = np.random.randn(N_COLL, N_DIM)")
print("    J_B = np.random.randn(N_COLL, N_DIM)")
print("\nVerdict:")
print("  J_A / J_B original mathematical definition = NOT_RECOVERABLE")
print("  No formula is invented here; all originals remain unverified.")

raise SystemExit(0)
