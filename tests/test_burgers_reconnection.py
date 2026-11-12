"""Tests for exact Burgers-boundary preservation under local reconnection."""
import pytest

from src.burgers_reconnection import (
    OrientedBurgersEdge,
    boundary_flux,
    junction_is_burgers_conservative,
    reconnect_equal_burgers,
)


def test_equal_burgers_head_exchange_preserves_vertex_boundary_exactly():
    first = OrientedBurgersEdge("a", "b", (1, 0))
    second = OrientedBurgersEdge("c", "d", (1, 0))
    after = reconnect_equal_burgers(first, second)
    assert after == (
        OrientedBurgersEdge("a", "d", (1, 0)),
        OrientedBurgersEdge("c", "b", (1, 0)),
    )
    assert boundary_flux((first, second)) == boundary_flux(after)


def test_unequal_burgers_head_exchange_is_rejected():
    first = OrientedBurgersEdge("a", "b", (1, 0))
    second = OrientedBurgersEdge("c", "d", (0, 1))
    with pytest.raises(ValueError, match="equal Burgers"):
        reconnect_equal_burgers(first, second)


def test_junction_conservation_uses_vector_sum_not_magnitudes():
    assert junction_is_burgers_conservative(
        incoming=((1, 0), (0, 1)),
        outgoing=((1, 1),),
    )
    assert not junction_is_burgers_conservative(
        incoming=((1, 0), (0, 1)),
        outgoing=((2, 0),),
    )


def test_junction_dimension_mismatch_is_rejected():
    with pytest.raises(ValueError, match="dimension"):
        junction_is_burgers_conservative(((1, 0),), ((1,),))


def test_boundary_flux_rejects_self_loops():
    with pytest.raises(ValueError, match="Self-loop"):
        boundary_flux((OrientedBurgersEdge("a", "a", (1, 0)),))
