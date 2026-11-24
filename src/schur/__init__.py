"""Schur complement and effective action.

Implements:
- Schur complement formula
- Effective Hamiltonian
- Memory kernel (placeholder)
"""

import numpy as np
from numpy.typing import NDArray


def schur_complement(
    K: NDArray,
    lambda_param: float,
    indices_o: list,
    indices_n: list,
) -> NDArray:
    """Compute Schur complement of matrix K.

    K_eff(λ) = K_oo - λI - K_on(K_nn - λI)^(-1)K_no

    Args:
        K: Full matrix
        lambda_param: Spectral parameter
        indices_o: Indices for observable block
        indices_n: Indices for noise block

    Returns:
        K_eff: Effective (N_o × N_o) matrix
    """
    K_oo = K[np.ix_(indices_o, indices_o)]
    K_on = K[np.ix_(indices_o, indices_n)]
    K_no = K[np.ix_(indices_n, indices_o)]
    K_nn = K[np.ix_(indices_n, indices_n)]

    N_o = len(indices_o)
    I = np.eye(N_o)

    K_nn_minus_lambda = K_nn - lambda_param * np.eye(len(indices_n))

    try:
        inv = np.linalg.inv(K_nn_minus_lambda)
    except np.linalg.LinAlgError:
        # Use pseudo-inverse if singular
        inv = np.linalg.pinv(K_nn_minus_lambda)

    K_eff = K_oo - lambda_param * I - K_on @ inv @ K_no

    return K_eff
