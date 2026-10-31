import numpy as np

from src.binary_icosahedral import (
    IDENTITY, audit_binary_icosahedral, binary_icosahedral_elements,
    qconj, qmul, quaternion_to_so3,
)


def test_binary_icosahedral_exhaustive_audit():
    report = audit_binary_icosahedral()
    assert report["order_2I"] == 120
    assert report["closure_missing_products"] == 0
    assert report["inverse_missing"] == 0
    assert report["distinct_SO3_rotations"] == 60


def test_quaternion_inverse_and_identity():
    elements = binary_icosahedral_elements()
    for q in elements:
        assert np.allclose(qmul(q, qconj(q)), IDENTITY, atol=1e-12)
        assert np.allclose(qmul(IDENTITY, q), q, atol=1e-12)


def test_double_cover_q_and_minus_q_same_rotation():
    for q in binary_icosahedral_elements():
        assert np.allclose(quaternion_to_so3(q), quaternion_to_so3(-q), atol=1e-12)
