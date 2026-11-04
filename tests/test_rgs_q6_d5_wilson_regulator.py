import numpy as np
import pytest

from src.rgs_q6_d5_wilson_regulator import (
    q6_wilson_spectral_cutoff,
    audit_q6_wilson_d5_regulator_obstruction,
)
from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data
from src.rgs_q6_u1_character_electric_flux import build_q6_integer_cycle_basis


def test_exact_transfer_spectral_cutoff_equivalence():
    data = build_q6_euclidean_time_data(2)
    cycles = build_q6_integer_cycle_basis(data)
    states = np.vstack((np.zeros((1, 192), dtype=np.int64), cycles, -cycles))
    result = q6_wilson_spectral_cutoff(states, 80.0, 0.05, 0.25, data)
    assert result["selected_count"] >= 1
    assert result["selected_mask"][0]
    assert np.array_equal(
        result["selected_mask"],
        result["electric_energies"] <= result["equivalent_maximum_electric_energy"] + 1e-12,
    )
    assert np.all(result["transfer_weights"] > 0.0)
    assert np.all(result["transfer_weights"] <= 1.0 + 1e-12)


def test_q6_d5_projection_obstruction_audit():
    result = audit_q6_wilson_d5_regulator_obstruction()
    assert result["all_checks_passed"], result["checks"]
    assert result["five_mode_projection"]["rank"] == 5
    assert result["five_mode_projection"]["kernel_dimension"] == 124
    assert result["transfer"]["total_states"] == 259


@pytest.mark.parametrize(
    "beta,at,tau",
    [(0.0, 0.05, 0.2), (80.0, 0.0, 0.2), (80.0, 0.05, 0.0),
     (80.0, 0.05, 1.1), (float("nan"), 0.05, 0.2)],
)
def test_invalid_spectral_cutoff_parameters(beta, at, tau):
    with pytest.raises(ValueError):
        q6_wilson_spectral_cutoff(np.zeros((1, 192), dtype=np.int64), beta, at, tau)
