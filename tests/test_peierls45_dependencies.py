import numpy as np

from experiments.peierls import peierls45


def test_peierls45_resolves_helpers_from_peierls46():
    assert callable(peierls45.generate_amorphous_lattice)
    assert callable(peierls45.build_hamiltonian)
    assert callable(peierls45.find_local_optimal_gap)


def test_peierls45_small_lattice_hamiltonian_is_hermitian():
    coords, box_size = peierls45.generate_amorphous_lattice(
        48, density=1.0, steps=2, seed=123
    )
    assert coords.shape == (48, 2)
    assert np.isfinite(coords).all()
    assert np.all(coords >= 0.0)
    assert np.all(coords <= box_size)
    hamiltonian = peierls45.build_hamiltonian(
        coords, box_size, B_field=0.5, gamma=0.3
    )
    dense = hamiltonian.toarray()
    assert dense.shape == (48, 48)
    assert np.isfinite(dense).all()
    assert np.allclose(dense, dense.conj().T, atol=1e-12)
    assert np.count_nonzero(dense) > 0


def test_peierls45_gap_selector_returns_valid_occupied_count():
    eigenvalues = np.sort(
        np.r_[np.linspace(-2.0, -0.1, 39), 0.8, np.linspace(1.0, 2.0, 40)]
    )
    occupied = peierls45.find_local_optimal_gap(
        eigenvalues, target_fraction=0.5, window=0.2
    )
    assert isinstance(occupied, (int, np.integer))
    assert 1 <= occupied < len(eigenvalues)
