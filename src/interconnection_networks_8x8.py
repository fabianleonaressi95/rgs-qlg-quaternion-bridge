"""Exhaustive finite models of 8x8 Omega and Beneš interconnection networks.

Wire labels are 3-bit integers. In shuffle-normalized coordinates, a stage
whose dimension is d contains four independent 2x2 switches. Switch IDs are
the two remaining wire bits, packed from low to high.

Omega-8 uses dimensions (0, 1, 2), equivalent to three adjacent-switch stages
with perfect-shuffle interstage wiring under stage-wise wire relabelling.
Beneš-8 uses (0, 1, 2, 1, 0).

This module enumerates distinct input-output permutations, not raw switch
configurations. Multiple configurations may realize the same permutation.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial
from typing import Dict, Iterable, Tuple

N = 8
BITS = 3
OMEGA_DIMS = (0, 1, 2)
BENES_DIMS = (0, 1, 2, 1, 0)
Permutation = Tuple[int, ...]
StageSettings = Tuple[int, ...]


def _switch_id(wire: int, dimension: int) -> int:
    packed = 0
    out_bit = 0
    for bit in range(BITS):
        if bit != dimension:
            packed |= ((wire >> bit) & 1) << out_bit
            out_bit += 1
    return packed


def _stage_transform(dimension: int, settings: int) -> Tuple[int, ...]:
    if dimension not in range(BITS):
        raise ValueError("dimension must be 0, 1, or 2")
    if not 0 <= settings < 16:
        raise ValueError("stage settings must be a 4-bit integer")
    return tuple(
        wire ^ (((settings >> _switch_id(wire, dimension)) & 1) << dimension)
        for wire in range(N)
    )


def validate_permutation(permutation: Iterable[int]) -> Permutation:
    result = tuple(permutation)
    if len(result) != N or set(result) != set(range(N)):
        raise ValueError("expected a permutation of wires 0..7")
    return result


@dataclass(frozen=True)
class NetworkEnumeration:
    name: str
    dimensions: Tuple[int, ...]
    witnesses: Dict[Permutation, StageSettings]
    layer_counts: Tuple[int, ...]

    @property
    def permutation_count(self) -> int:
        return len(self.witnesses)

    @property
    def total_permutations(self) -> int:
        return factorial(N)

    def routes(self, permutation: Iterable[int]) -> bool:
        return validate_permutation(permutation) in self.witnesses

    def settings_for(self, permutation: Iterable[int]) -> StageSettings:
        target = validate_permutation(permutation)
        try:
            return self.witnesses[target]
        except KeyError as exc:
            raise ValueError(f"{self.name} cannot realize permutation {target}") from exc


def enumerate_network(name: str, dimensions: Iterable[int]) -> NetworkEnumeration:
    dims = tuple(dimensions)
    if not dims or any(d not in range(BITS) for d in dims):
        raise ValueError("dimensions must be a nonempty sequence from {0,1,2}")

    identity = tuple(range(N))
    states: Dict[Permutation, StageSettings] = {identity: ()}
    counts = [1]

    for dimension in dims:
        next_states: Dict[Permutation, StageSettings] = {}
        transforms = tuple(_stage_transform(dimension, mask) for mask in range(16))
        for current, witness in states.items():
            for mask, transform in enumerate(transforms):
                output = tuple(transform[wire] for wire in current)
                if output not in next_states:
                    next_states[output] = witness + (mask,)
        states = next_states
        counts.append(len(states))

    return NetworkEnumeration(name, dims, states, tuple(counts))


def omega8() -> NetworkEnumeration:
    return enumerate_network("Omega-8", OMEGA_DIMS)


def benes8() -> NetworkEnumeration:
    return enumerate_network("Benes-8", BENES_DIMS)


def bit_reversal_permutation() -> Permutation:
    return tuple(int(f"{wire:03b}"[::-1], 2) for wire in range(N))


def apply_settings(
    permutation: Iterable[int], dimensions: Iterable[int], settings: Iterable[int]
) -> Permutation:
    current = validate_permutation(permutation)
    dims = tuple(dimensions)
    masks = tuple(settings)
    if len(dims) != len(masks):
        raise ValueError("one settings mask is required for each stage")
    for dimension, mask in zip(dims, masks):
        transform = _stage_transform(dimension, mask)
        current = tuple(transform[wire] for wire in current)
    return current
