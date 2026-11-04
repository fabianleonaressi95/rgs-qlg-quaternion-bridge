import json
import numpy as np
import pytest

from scripts.rgs_psu_integrated_audit import _json_default


def test_json_default_preserves_complex_real_and_imaginary_parts():
    payload = {"scalar": 1 + 2j, "numpy_scalar": np.complex128(3 - 4j)}
    encoded = json.dumps(payload, default=_json_default)
    decoded = json.loads(encoded)
    assert decoded["scalar"] == {"real": 1.0, "imag": 2.0}
    assert decoded["numpy_scalar"] == {"real": 3.0, "imag": -4.0}


def test_json_default_serializes_numpy_arrays_and_scalars():
    payload = {"array": np.array([1.0, 2.0]), "integer": np.int64(7)}
    decoded = json.loads(json.dumps(payload, default=_json_default))
    assert decoded == {"array": [1.0, 2.0], "integer": 7}


def test_json_default_rejects_unknown_objects():
    with pytest.raises(TypeError):
        _json_default(object())
