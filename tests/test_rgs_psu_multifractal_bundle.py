import numpy as np
from math import comb

from src.rgs_psu_multifractal_bundle import (
    audit_rgs_vector_bundle_connection,
    build_rgs_vector_bundle_connection,
    face_holonomies,
    gauge_transform_connection,
    multifractal_spectrum,
)


def test_seed_induced_rank_three_connection_is_metric_compatible_and_flat():
    report = audit_rgs_vector_bundle_connection()
    assert report["base_vertices"] == 8
    assert report["base_edges"] == 12
    assert report["base_square_faces"] == 6
    assert report["fibre_rank"] == 3
    assert report["metric_compatible"]
    assert report["links_are_invertible"]
    assert report["canonical_connection_is_flat"]
    assert report["max_face_holonomy_residual"] < 1e-12


def test_local_gl3_gauge_transform_preserves_flatness_and_metric_compatibility():
    connection = build_rgs_vector_bundle_connection()
    rng = np.random.default_rng(20261104)
    gauges = np.asarray([np.eye(3) + 0.08 * rng.normal(size=(3, 3))
                         for _ in range(8)])
    transformed = gauge_transform_connection(connection, gauges)
    report = audit_rgs_vector_bundle_connection(transformed)
    assert report["metric_compatible"]
    assert report["links_are_invertible"]
    assert report["canonical_connection_is_flat"]
    original = face_holonomies(connection)
    changed = face_holonomies(transformed)
    for k, (a, *_rest) in enumerate(connection.base_faces):
        expected = gauges[a] @ original[k] @ np.linalg.inv(gauges[a])
        assert np.allclose(changed[k], expected, atol=1e-11)


def test_uniform_nested_measure_has_one_dimensional_generalized_dimensions():
    scales = 0.5 ** np.arange(1, 8)
    measures = [np.full(2**n, 1.0 / 2**n) for n in range(1, 8)]
    report = multifractal_spectrum(scales, measures, q_values=(-1, 0, 0.5, 1, 2, 3))
    assert all(abs(value - 1.0) < 1e-12 for value in report["D_q"].values())


def test_binomial_cascade_has_nonconstant_generalized_dimensions():
    p = 0.73
    scales = 0.5 ** np.arange(1, 9)
    measures = []
    for n in range(1, 9):
        masses = np.asarray([
            p ** int(k.bit_count()) * (1.0 - p) ** (n - int(k.bit_count()))
            for k in range(2**n)
        ])
        measures.append(masses)
    report = multifractal_spectrum(scales, measures, q_values=(0, 0.5, 1, 2, 3))
    values = np.asarray(list(report["D_q"].values()))
    assert np.ptp(values) > 0.05
