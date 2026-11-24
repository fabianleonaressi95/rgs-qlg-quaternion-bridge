import pytest
from src.q3_q6_network_bridge import (
    audit, int_to_trigram, trigram_to_int, q6_demand_index,
    q6_index_to_demand, permutation_to_q6_demands,
    q6_demands_to_permutation, q3_coordinate_matching, q3_edges, route_trace,
)
from src.interconnection_networks_8x8 import omega8, benes8, bit_reversal_permutation

def test_all_trigrams_and_q6_demands_round_trip():
    for i in range(8):
        assert trigram_to_int(int_to_trigram(i)) == i
    for lower in range(8):
        for upper in range(8):
            index = q6_demand_index(int_to_trigram(lower), int_to_trigram(upper))
            assert q6_index_to_demand(index) == (int_to_trigram(lower), int_to_trigram(upper))
            assert index == lower + 8 * upper

def test_q3_edges_are_three_disjoint_coordinate_matchings():
    matchings = [q3_coordinate_matching(d) for d in range(3)]
    assert [len(m) for m in matchings] == [4, 4, 4]
    assert len(q3_edges()) == 12 and len(set(q3_edges())) == 12
    for d, matching in enumerate(matchings):
        assert all((a ^ b) == (1 << d) for a, b in matching)

def test_permutation_matrix_is_q6_perfect_matching():
    for p in (tuple(range(8)), bit_reversal_permutation(), (7, 0, 5, 2, 1, 6, 3, 4)):
        demands = permutation_to_q6_demands(p)
        assert len(demands) == 8
        assert q6_demands_to_permutation(demands) == p

def test_reject_nonmatching_q6_demands():
    identity = permutation_to_q6_demands(tuple(range(8)))
    with pytest.raises(ValueError): q6_demands_to_permutation(identity[:-1])
    with pytest.raises(ValueError): q6_demands_to_permutation(identity[:-1] + (identity[0],))
    with pytest.raises(ValueError): trigram_to_int("10x")

def test_all_network_permutations_have_distinct_q6_matching_encodings():
    for network in (omega8(), benes8()):
        encodings = {tuple(permutation_to_q6_demands(p)) for p in network.witnesses}
        assert len(encodings) == network.permutation_count

def test_route_trace_replays_witness_to_each_q6_demand():
    network = benes8()
    target = (7, 0, 5, 2, 1, 6, 3, 4)
    settings = network.settings_for(target)
    for source, destination in enumerate(target):
        trace = route_trace(source, network.dimensions, settings)
        assert trace[0] == source and trace[-1] == destination

def test_audit_records_prism_degree_obstruction():
    result = audit()
    assert result["q3"]["vertices"] == 8 and result["q3"]["edges"] == 12
    assert result["q6"]["vertices"] == 64
    assert result["omega"]["q6_demand_matchings"] == 4096
    assert result["benes"]["q6_demand_matchings"] == 40320
    assert result["prism_c64_factor"]["graph_isomorphism"] is False
