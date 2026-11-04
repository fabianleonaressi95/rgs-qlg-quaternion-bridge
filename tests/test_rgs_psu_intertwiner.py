import numpy as np

from src.rgs_psu_intertwiner import audit_rgs_golden_intertwiner


def test_archived_rgs_and_golden_q6_are_not_affinely_equivalent():
    report = audit_rgs_golden_intertwiner()
    assert report["point_counts"] == {"rgs": 64, "golden_projected_q6": 64}
    assert report["population_covariance"]["rgs_isotropy_relative_residual"] < 1e-12
    assert report["population_covariance"]["q6_isotropy_relative_residual"] < 1e-12
    assert report["affine_no_go"]["covariance_argument_applies"]
    assert report["pair_distance_invariant"]["normalized_sorted_pair_distance_rms_residual"] > 1e-3
    assert report["affine_no_go"]["affine_equivalence_ruled_out"]
    assert report["epistemic_status"] == "NUMERICALLY_VERIFIED_NO_GO_FOR_AFFINE_IDENTIFICATION"


def test_rgs_golden_no_go_is_translation_and_rotation_invariant():
    from src.rgs_ellis_refinement import load_rgs_seed_data

    x = load_rgs_seed_data().coordinates
    q = np.linalg.qr(np.asarray([[1.0, 2.0, 3.0], [2.0, -1.0, 1.0], [1.0, 1.0, -1.0]]))[0]
    transformed = 2.7 * x @ q + np.asarray([4.0, -3.0, 2.0])
    report = audit_rgs_golden_intertwiner(transformed)
    assert report["affine_no_go"]["affine_equivalence_ruled_out"]
    assert np.isclose(
        report["pair_distance_invariant"]["normalized_sorted_pair_distance_rms_residual"],
        0.06920925965227274,
        atol=1e-10,
    )
