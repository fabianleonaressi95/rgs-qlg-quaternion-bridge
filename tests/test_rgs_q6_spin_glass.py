import numpy as np
import pytest

from src.rgs_q6_spin_glass import Q6SpinGlass, exact_ground_state, parallel_tempering


def test_q6_quenched_bonds_are_reproducible_and_have_expected_shape():
    a = Q6SpinGlass.random(seed=17, distribution="gaussian")
    b = Q6SpinGlass.random(seed=17, distribution="gaussian")
    assert a.edges.shape == (192, 2)
    assert a.faces.shape == (240, 4)
    assert a.couplings.shape == (192,)
    assert np.array_equal(a.couplings, b.couplings)
    assert not np.array_equal(a.couplings, Q6SpinGlass.random(seed=18).couplings)


@pytest.mark.parametrize("distribution", ["gaussian", "bimodal"])
def test_frustration_and_energy_are_z2_gauge_invariant(distribution):
    model = Q6SpinGlass.random(seed=31, distribution=distribution)
    rng = np.random.default_rng(32)
    spins = rng.choice([-1, 1], size=64).astype(np.int8)
    eps = rng.choice([-1, 1], size=64).astype(np.int8)
    transformed = model.gauge_transform(eps)
    assert model.energy(spins) == pytest.approx(transformed.energy(eps * spins), abs=1e-12)
    assert np.array_equal(model.frustration(), transformed.frustration())


def test_single_spin_flip_energy_increment_matches_hamiltonian():
    model = Q6SpinGlass.random(seed=41)
    rng = np.random.default_rng(42)
    spins = rng.choice([-1, 1], size=64).astype(np.int8)
    fields = model.local_fields(spins)
    for site in (0, 7, 31, 63):
        flipped = spins.copy()
        flipped[site] *= -1
        delta = model.energy(flipped) - model.energy(spins)
        assert delta == pytest.approx(2 * spins[site] * fields[site], abs=1e-12)


def test_exact_ground_state_works_on_small_induced_graph():
    model = Q6SpinGlass(
        edges=np.array([[0, 1], [1, 2], [0, 2]], dtype=np.int16),
        faces=np.empty((0, 4), dtype=np.int16),
        couplings=np.array([1.0, 1.0, -1.0]),
        n_vertices=3,
    )
    result = exact_ground_state(model)
    assert result["enumerated_states"] == 8
    assert result["energy"] == pytest.approx(-1.0)


def test_exact_ground_state_refuses_full_q6_enumeration():
    with pytest.raises(ValueError, match="exact enumeration refused"):
        exact_ground_state(Q6SpinGlass.random(seed=51))


def test_parallel_tempering_returns_finite_overlap_observables():
    model = Q6SpinGlass.random(seed=61, distribution="bimodal")
    report = parallel_tempering(model, temperatures=[0.8, 1.5, 2.5],
                                sweeps=24, burn_in=8, measure_every=4, seed=62)
    assert report["n_vertices"] == 64
    assert report["n_edges"] == 192
    assert report["n_plaquettes"] == 240
    assert report["measurements"] == 4
    for key in ("mean_energy_per_spin", "replica_overlap_second_moment",
                "spin_glass_susceptibility", "frustration_fraction"):
        assert np.all(np.isfinite(np.asarray(report[key], dtype=float)))
    assert np.all(np.asarray(report["spin_glass_susceptibility"]) >= 0)
