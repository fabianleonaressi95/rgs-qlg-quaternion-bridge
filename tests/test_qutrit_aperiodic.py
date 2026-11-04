import numpy as np
import pytest

from src.qutrit_aperiodic import (
    aperiodic_prefix_audit,
    audit_qutrit_aperiodic_model,
    build_fibonacci_qutrit_chain,
    fibonacci_word,
    ground_state_entanglement,
    qutrit_bell_state,
    schmidt_spectrum,
    von_neumann_entropy,
)


def test_fibonacci_substitution_prefix_and_finite_period_check():
    assert fibonacci_word(8) == "ABAABABA"
    report = aperiodic_prefix_audit(fibonacci_word(8))
    assert report["is_periodic_as_finite_word"] is False
    assert report["finite_word_periods"] == []


def test_qutrit_maximally_entangled_bell_state():
    probabilities = schmidt_spectrum(qutrit_bell_state(), (3, 3))
    assert np.allclose(probabilities, np.full(3, 1/3))
    assert von_neumann_entropy(probabilities) == pytest.approx(np.log2(3))
    assert np.linalg.norm(qutrit_bell_state()) == pytest.approx(1.0)


def test_fibonacci_qutrit_chain_is_hermitian_and_entangled():
    chain = build_fibonacci_qutrit_chain(4, j_a=1.0, j_b=0.37)
    assert chain.word == "ABA"
    assert chain.hamiltonian.shape == (81, 81)
    assert np.allclose(chain.hamiltonian, chain.hamiltonian.T)
    result = ground_state_entanglement(chain)
    assert result["gap"] >= -1e-10
    assert result["schmidt_rank"] >= 2
    assert result["entanglement_entropy_bits"] > 1e-8
    assert result["hamiltonian_hermiticity_residual"] < 1e-12


def test_combined_qutrit_aperiodic_audit():
    report = audit_qutrit_aperiodic_model(4)
    assert report["local_hilbert_dimension"] == 3
    assert report["bell_state_entropy_bits"] == pytest.approx(np.log2(3))
    assert report["ground_state"]["word"] == "ABA"


def test_invalid_chain_size_and_couplings_rejected():
    with pytest.raises(ValueError):
        build_fibonacci_qutrit_chain(1)
    with pytest.raises(ValueError):
        build_fibonacci_qutrit_chain(3, j_a=0.0)
