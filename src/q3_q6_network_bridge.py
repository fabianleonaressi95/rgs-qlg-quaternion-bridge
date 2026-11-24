"""Exact finite bridge: I Ching Q3 trigrams, Q6 demands and 8x8 switch networks.

Bit convention: each trigram is a three-character string in bottom-to-top
line order. A Q6 word is lower_trigram || upper_trigram; its integer label is
lower + 8*upper, so the lower trigram occupies bits 0..2 and the upper bits 3..5.

A Q6 vertex is an input/output demand (source trigram, destination trigram).
A permutation is a set of eight demands forming a perfect matching of K_8,8.
This is a bijective encoding of permutation matrices, not an identification
of Q6 with the directed switching fabric.
"""
from __future__ import annotations
from typing import Iterable, Sequence, Tuple
from src.interconnection_networks_8x8 import (
    BITS, N, NetworkEnumeration, _stage_transform, omega8, benes8,
)

Trigram = str
Demand = Tuple[Trigram, Trigram]

def validate_trigram(value: str) -> Trigram:
    if not isinstance(value, str) or len(value) != 3 or set(value) - {"0", "1"}:
        raise ValueError("a trigram must be exactly three bits, in bottom-to-top order")
    return value

def trigram_to_int(value: str) -> int:
    return int(validate_trigram(value), 2)

def int_to_trigram(value: int) -> Trigram:
    if not isinstance(value, int) or not 0 <= value < 8:
        raise ValueError("trigram index must be an integer in 0..7")
    return f"{value:03b}"

def q6_demand_index(lower: str, upper: str) -> int:
    return trigram_to_int(lower) | (trigram_to_int(upper) << 3)

def q6_index_to_demand(index: int) -> Demand:
    if not isinstance(index, int) or not 0 <= index < 64:
        raise ValueError("Q6 index must be an integer in 0..63")
    return int_to_trigram(index & 7), int_to_trigram(index >> 3)

def permutation_to_q6_demands(permutation: Iterable[int]) -> Tuple[int, ...]:
    p = tuple(permutation)
    if len(p) != N or set(p) != set(range(N)):
        raise ValueError("expected a permutation of wires 0..7")
    return tuple(q6_demand_index(int_to_trigram(i), int_to_trigram(p[i])) for i in range(N))

def q6_demands_to_permutation(indices: Iterable[int]) -> Tuple[int, ...]:
    values = tuple(indices)
    if len(values) != N or len(set(values)) != N:
        raise ValueError("a permutation demand must contain eight distinct Q6 vertices")
    pairs = tuple(q6_index_to_demand(v) for v in values)
    inputs = tuple(trigram_to_int(a) for a, _ in pairs)
    outputs = tuple(trigram_to_int(b) for _, b in pairs)
    if set(inputs) != set(range(N)) or set(outputs) != set(range(N)):
        raise ValueError("Q6 demands must use every input and output trigram exactly once")
    result = [0] * N
    for a, b in pairs:
        result[trigram_to_int(a)] = trigram_to_int(b)
    return tuple(result)

def q3_coordinate_matching(dimension: int) -> Tuple[Tuple[int, int], ...]:
    if dimension not in range(BITS):
        raise ValueError("dimension must be 0, 1, or 2")
    return tuple((v, v ^ (1 << dimension)) for v in range(N)
                 if not (v >> dimension) & 1)

def q3_edges() -> Tuple[Tuple[int, int], ...]:
    return tuple(edge for d in range(BITS) for edge in q3_coordinate_matching(d))

def route_trace(input_wire: int, dimensions: Sequence[int],
                settings: Sequence[int]) -> Tuple[int, ...]:
    if not isinstance(input_wire, int) or input_wire not in range(N):
        raise ValueError("input wire must be in 0..7")
    dims, masks = tuple(dimensions), tuple(settings)
    if len(dims) != len(masks):
        raise ValueError("one settings mask is required per stage")
    position, trace = input_wire, [input_wire]
    for dimension, mask in zip(dims, masks):
        transform = _stage_transform(dimension, mask)
        position = transform[position]
        trace.append(position)
    return tuple(trace)

def network_q6_audit(network: NetworkEnumeration) -> dict:
    encoded = {tuple(permutation_to_q6_demands(p)) for p in network.witnesses}
    if len(encoded) != network.permutation_count:
        raise ArithmeticError("permutation-to-Q6 encoding is not injective")
    return {"network": network.name, "routable_permutations": network.permutation_count,
            "q6_demand_matchings": len(encoded),
            "bijection_verified": len(encoded) == network.permutation_count}

def audit() -> dict:
    omega, benes = omega8(), benes8()
    return {
        "bit_convention": "three-bit trigram, bottom-to-top; Q6=lower||upper",
        "q3": {"vertices": 8, "edges": len(q3_edges()), "degree": 3,
               "switch_pairs_per_coordinate": [
                   [list(e) for e in q3_coordinate_matching(d)] for d in range(3)]},
        "q6": {"vertices": 64, "meaning": "ordered input/output trigram demand",
               "permutation_matrix_ones": 8},
        "omega": network_q6_audit(omega),
        "benes": network_q6_audit(benes),
        "prism_c64_factor": {
            "label_bijection": "Q6 indices 0..63 may label C64 positions 0..63",
            "graph_isomorphism": False,
            "obstruction": "Q6 has degree 6; cycle C64 has degree 2",
            "status": "coordinate relabelling only; no edge-preserving Q6-to-C64 identification"},
        "scope": "finite combinatorial correspondence; no physical or dynamical equivalence",
    }
