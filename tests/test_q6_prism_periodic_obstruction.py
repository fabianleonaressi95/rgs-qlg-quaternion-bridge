from itertools import product
from src.q6_prism_periodic_obstruction import (
    PRISM_CYCLE_LENGTHS, all_closed_walks, audit, winding_vector,
)

def test_all_closed_four_step_walks_have_zero_winding():
    closed = list(all_closed_walks(4))
    assert len(closed) == 270
    assert all(w == (0, 0, 0, 0, 0) for _, w in closed)

def test_no_four_step_walk_can_wind_any_factor():
    assert min(PRISM_CYCLE_LENGTHS) == 5
    # A winding number m_i != 0 requires at least n_i steps in factor i.
    for n in PRISM_CYCLE_LENGTHS:
        assert n > 4

def test_five_steps_can_wind_the_shortest_factor():
    walk = [(4, 1)] * 5
    assert winding_vector(walk) == (0, 0, 0, 0, 1)

def test_closed_walk_must_close_mod_each_factor():
    assert winding_vector([(4, 1)] * 4) is None
    assert winding_vector([(4, 1)] * 5) == (0, 0, 0, 0, 1)

def test_q6_graph_cycle_rank_and_audit():
    result = audit()
    assert result["q6_graph_cycle_rank"] == 129
    assert result["q6_square_plaquettes"] == 240
    assert result["four_step_closed_walks_with_nonzero_winding"] == 0
    assert result["minimum_steps_needed_for_nonzero_winding"] == 5
