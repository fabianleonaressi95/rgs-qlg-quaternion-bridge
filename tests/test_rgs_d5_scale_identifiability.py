import numpy as np
import pytest

from src.rgs_d5_scale_identifiability import (
    audit_rgs_d5_scale_identifiability,
    local_covariance_coframe,
)
from src.rgs_ellis_refinement import load_rgs_seed_data


def test_fixed_kernel_length_breaks_naive_scale_covariance():
    x = load_rgs_seed_data().coordinates
    c, e = local_covariance_coframe(x)
    s = 1.7
    cs, es = local_covariance_coframe(s*x)
    # ell2 is a squared length. Holding it fixed while scaling X changes
    # exp(-|delta X|^2/ell2), so neither naive homogeneity law is valid.
    assert not np.allclose(cs, s*s*c, rtol=2e-12, atol=2e-12)
    assert not np.allclose(es, s*e, rtol=2e-12, atol=2e-12)


def test_coframe_is_invariant_when_all_length_scales_are_scaled_together():
    x = load_rgs_seed_data().coordinates
    c, e = local_covariance_coframe(x)
    s = 1.7
    cs, es = local_covariance_coframe(s*x, ell2=s*s*7.035514834154695,
                                     c0=s*s*43.96495130838031)
    assert np.allclose(cs, c, rtol=2e-12, atol=2e-12)
    assert np.allclose(es, e, rtol=2e-12, atol=2e-12)


def test_current_q6_embedding_cannot_induce_a_positive_definite_d5_metric():
    report = audit_rgs_d5_scale_identifiability()
    assert report["q6_to_r3_rank"]["rank_max"] <= 3
    assert report["d5_target"]["real_tangent_dimension"] == 10
    assert not report["d5_target"]["positive_definite_pullback_possible_from_current_q6_to_r3_map"]
    assert report["checks"] and all(report["checks"].values())


@pytest.mark.parametrize("bad", [0.0, -1.0, np.inf, np.nan])
def test_scale_probe_must_be_finite_positive(bad):
    with pytest.raises(ValueError):
        audit_rgs_d5_scale_identifiability(bad)
