import numpy as np

from network_dynamics import (
    adjacency_matrix,
    build_networks,
    make_skill_profiles,
    neighbor_mixing_matrix,
    run_experiment,
    simulate_knowledge,
)


def test_networks_have_equal_node_and_edge_counts():
    profiles = make_skill_profiles(np.random.default_rng(7), 30, 8)
    networks = build_networks(profiles, seed=8, n_edges=45)
    assert set(networks) == {"hierarchical", "matrix", "practices_like"}
    assert all(len(edges) == 45 for edges in networks.values())
    assert all(i != j for edges in networks.values() for i, j in edges)


def test_mixing_matrix_is_row_stochastic():
    adjacency = adjacency_matrix(4, {(0, 1), (1, 2)})
    mixing = neighbor_mixing_matrix(adjacency)
    np.testing.assert_allclose(mixing.sum(axis=1), np.ones(4))
    assert np.all(mixing >= 0.0)


def test_knowledge_state_is_bounded_and_reproducible():
    profiles = make_skill_profiles(np.random.default_rng(11), 12, 4)
    initial = profiles * 0.25
    adjacency = adjacency_matrix(12, {(i, i + 1) for i in range(11)})
    mixing = neighbor_mixing_matrix(adjacency)
    first = simulate_knowledge(mixing, initial, profiles, steps=20)
    second = simulate_knowledge(mixing, initial, profiles, steps=20)
    np.testing.assert_allclose(first, second)
    assert np.all((first >= 0.0) & (first <= 1.0))


def test_experiment_is_reproducible_and_uses_same_budget():
    first = run_experiment(seed=123, steps=12)
    second = run_experiment(seed=123, steps=12)
    assert first == second
    assert {row["edges"] for row in first} == {45}
    assert {row["nodes"] for row in first} == {30}
