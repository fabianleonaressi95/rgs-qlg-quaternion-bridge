from itertools import permutations

import pytest

from src.interconnection_networks_8x8 import (
    BENES_DIMS,
    OMEGA_DIMS,
    apply_settings,
    benes8,
    bit_reversal_permutation,
    omega8,
)


@pytest.fixture(scope="module")
def omega():
    return omega8()


@pytest.fixture(scope="module")
def benes():
    return benes8()


def test_exact_stage_architectures():
    assert OMEGA_DIMS == (0, 1, 2)
    assert BENES_DIMS == (0, 1, 2, 1, 0)
    assert len(OMEGA_DIMS) * 4 == 12
    assert len(BENES_DIMS) * 4 == 20


def test_omega_exhaustive_reachable_permutation_count(omega):
    assert omega.layer_counts == (1, 16, 256, 4096)
    assert omega.permutation_count == 4096
    assert omega.total_permutations == 40320


def test_benes_exhaustive_reachable_permutation_count(benes):
    assert benes.layer_counts == (1, 16, 256, 4096, 9216, 40320)
    assert benes.permutation_count == 40320
    assert benes.total_permutations == 40320


def test_benes_realizes_every_permutation(benes):
    for target in permutations(range(8)):
        assert benes.routes(target)


def test_bit_reversal_blocking_and_rearrangeability(omega, benes):
    reverse = bit_reversal_permutation()
    assert reverse == (0, 4, 2, 6, 1, 5, 3, 7)
    assert not omega.routes(reverse)
    assert benes.routes(reverse)


def test_returned_settings_reproduce_requested_permutation(benes):
    target = (7, 0, 5, 2, 1, 6, 3, 4)
    settings = benes.settings_for(target)
    assert apply_settings(range(8), BENES_DIMS, settings) == target


def test_omega_witnesses_reproduce_permutations(omega):
    for target, settings in omega.witnesses.items():
        assert apply_settings(range(8), OMEGA_DIMS, settings) == target


def test_invalid_permutation_rejected(omega):
    with pytest.raises(ValueError):
        omega.routes((0, 1, 2, 3, 4, 5, 6, 6))
