"""Mathematical constants for golden ratio geometry."""

import numpy as np

# Golden ratio
PHI = (1.0 + np.sqrt(5.0)) / 2.0

# Supergolden constant (solution to ψ³ = ψ + 1)
# ψ ≈ 1.324717957244746
PSI = ((9.0 + np.sqrt(69.0)) / 18.0) ** (1.0 / 3.0) + (
    (9.0 - np.sqrt(69.0)) / 18.0
) ** (1.0 / 3.0)

# Derived golden ratio constants
PHI_INVERSE = PHI - 1.0  # φ⁻¹ = φ - 1
PHI_SQUARED = PHI + 1.0  # φ² = φ + 1

# Verification constants
PHI_IDENTITY_ERROR = abs(PHI**2 - (PHI + 1.0))
PSI_IDENTITY_ERROR = abs(PSI**3 - (PSI + 1.0))


def verify_golden_identities(tol: float = 1e-10) -> dict:
    """Verify the defining identities of φ and ψ.

    Returns:
        Dictionary with verification results
    """
    return {
        "phi_squared_equals_phi_plus_one": PHI_IDENTITY_ERROR <= tol,
        "phi_squared_error": float(PHI_IDENTITY_ERROR),
        "psi_cubed_equals_psi_plus_one": PSI_IDENTITY_ERROR <= tol,
        "psi_cubed_error": float(PSI_IDENTITY_ERROR),
        "phi_value": float(PHI),
        "psi_value": float(PSI),
    }
