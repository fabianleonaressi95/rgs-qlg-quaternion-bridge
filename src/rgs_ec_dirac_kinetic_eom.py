"""Kinetic Euler operator from the repository's symmetrized EC Dirac form action.

The action convention is S_D=1/2 int [bar(psi) *gamma wedge D psi
+ D bar(psi) wedge *gamma psi], with eta=(-,+,+,+).  In an oriented
orthonormal frame, integration by parts and D gamma=0 give
delta S_D/delta bar(psi) = - gamma^a D_a psi (up to boundary terms).
No mass term is included because the repository's declared action above
does not specify one.  D_a psi is supplied by the caller, including the
chosen tetrad and spin connection.
"""
from __future__ import annotations

import numpy as np

from src.rgs_ec_algebraic_torsion import dirac_gamma_matrices_minus_plus
from src.rgs_ec_nonlinear_dirac_eom import (
    effective_dirac_contact_sources,
)


def dirac_kinetic_euler_operator(covariant_derivative_spinor):
    """Return K[psi]=-gamma^a D_a psi from the symmetrized form action.

    Input convention: derivative[:, a] is D_a psi in the local orthonormal
    frame. The caller is responsible for constructing the full covariant
    derivative (partial derivative plus spin connection).
    """
    derivative = np.asarray(covariant_derivative_spinor, dtype=complex)
    if derivative.shape != (4, 4) or not np.all(np.isfinite(derivative)):
        raise ValueError(
            "covariant_derivative_spinor must be a finite complex (4,4) array"
        )
    gamma = dirac_gamma_matrices_minus_plus()
    contracted = sum(gamma[a] @ derivative[:, a] for a in range(4))
    return {
        "kinetic_euler_operator": -contracted,
        "dirac_operator": contracted,
        "covariant_derivative_convention": "column a is D_a psi",
        "action_variation": "delta_S_D/delta_barpsi=-gamma^a D_a psi",
        "boundary_terms": "DISCARDED_WITH_COVARIANT_INTEGRATION_BY_PARTS",
        "mass_term_included": False,
        "status": "ACTION_DERIVED_MASSLESS_KINETIC_EULER_OPERATOR",
    }


def effective_ec_dirac_eom_residual(
    covariant_derivative_spinor, spinor, axial_current, kappa, *,
    bar_spinor=None,
):
    """Assemble the massless kinetic Euler term and EC axial contact source.

    The resulting residual is -gamma^a D_a psi + dL_4f/dbarpsi.
    It is a local equation residual, not a global solution or conservation
    audit. A mass term is intentionally not inferred from absent action data.
    """
    kinetic = dirac_kinetic_euler_operator(covariant_derivative_spinor)
    source = effective_dirac_contact_sources(
        spinor, axial_current, kappa, bar_spinor=bar_spinor
    )
    residual = kinetic["kinetic_euler_operator"] + source["dL_dbarpsi"]
    scale = max(
        float(np.linalg.norm(kinetic["kinetic_euler_operator"])),
        float(np.linalg.norm(source["dL_dbarpsi"])),
        1.0,
    )
    return {
        "residual": residual,
        "residual_norm": float(np.linalg.norm(residual)),
        "relative_residual": float(np.linalg.norm(residual) / scale),
        "kinetic_euler_operator": kinetic["kinetic_euler_operator"],
        "contact_euler_source": source["dL_dbarpsi"],
        "mass_term_included": False,
        "spin_connection_supplied_externally": True,
        "on_shell": bool(np.linalg.norm(residual) <= 1e-10 * scale),
        "conservation_status": "NOT_ESTABLISHED_BY_LOCAL_EOM_RESIDUAL",
    }
