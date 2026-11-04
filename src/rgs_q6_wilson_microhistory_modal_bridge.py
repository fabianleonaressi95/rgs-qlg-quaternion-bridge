"""Symmetric microhistory quotient and local modal semantics for Q6 Wilson fluxes.

For fixed signed edge occupancy n, the S_s-orbit of labelled event histories
has size s!/prod_e |n_e|!, and its stabilizer has size prod_e |n_e|!.
The Wilson Bessel coefficient is bounded by (beta/2)^s/prod_e |n_e|!.
This is a symmetrized counting identity, not a derivation of Gödel's
ontological argument or an exact quotient representation of the Wilson
measure.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from math import factorial
from typing import Callable, Sequence


def validate_signed_flux(flux: Sequence[int], edge_count: int = 192) -> tuple[int, ...]:
    values = tuple(flux)
    if len(values) != edge_count:
        raise ValueError(f"expected {edge_count} edge fluxes, got {len(values)}")
    if any(isinstance(x, bool) or not isinstance(x, int) for x in values):
        raise TypeError("edge fluxes must be integers")
    return values


def flux_l1_order(flux: Sequence[int], edge_count: int = 192) -> int:
    return sum(abs(x) for x in validate_signed_flux(flux, edge_count))


def labelled_microhistory(flux: Sequence[int], edge_count: int = 192) -> tuple[int, ...]:
    """Canonical event list; signed token +/- (edge_index+1) per unit flux."""
    n = validate_signed_flux(flux, edge_count)
    return tuple(token for edge, value in enumerate(n)
                 for token in ((edge + 1 if value > 0 else -(edge + 1),) * abs(value)))


def microhistory_orbit_stabilizer(flux: Sequence[int], edge_count: int = 192) -> dict:
    """Exact orbit-stabilizer data for relabelling the s event slots by S_s."""
    n = validate_signed_flux(flux, edge_count)
    multiplicities = tuple(abs(x) for x in n if x)
    s = sum(multiplicities)
    stabilizer = 1
    for multiplicity in multiplicities:
        stabilizer *= factorial(multiplicity)
    orbit = factorial(s) // stabilizer
    assert orbit * stabilizer == factorial(s)
    return {
        "order": s,
        "nonzero_multiplicities": multiplicities,
        "stabilizer_order": stabilizer,
        "orbit_order": orbit,
        "s_factorial": factorial(s),
        "orbit_over_s_factorial": Fraction(orbit, factorial(s)),
        "inverse_stabilizer": Fraction(1, stabilizer),
        "orbit_stabilizer_identity": True,
    }


def wilson_factorial_majorant(flux: Sequence[int], beta: Fraction | int | str,
                               edge_count: int = 192) -> Fraction:
    """Exact rational majorant (beta/2)^s / prod_e |n_e|!."""
    n = validate_signed_flux(flux, edge_count)
    b = Fraction(beta)
    if b < 0:
        raise ValueError("beta must be nonnegative")
    s = sum(abs(x) for x in n)
    denominator = 1
    for value in n:
        denominator *= factorial(abs(value))
    return (b / 2) ** s / denominator


def symmetrized_microhistory_weight(flux: Sequence[int], beta: Fraction | int | str,
                                    edge_count: int = 192) -> Fraction:
    """(beta/2)^s * |Orbit|/s!, exactly the factorial majorant."""
    n = validate_signed_flux(flux, edge_count)
    b = Fraction(beta)
    data = microhistory_orbit_stabilizer(n, edge_count)
    return (b / 2) ** data["order"] * data["orbit_over_s_factorial"]


def enumerate_orbit(history: Sequence[int]) -> frozenset[tuple[int, ...]]:
    """Finite exact orbit of a small labelled history under S_s."""
    return frozenset(permutations(tuple(history)))


def is_orbit_invariant(predicate: Callable[[tuple[int, ...]], bool],
                        history: Sequence[int]) -> bool:
    return len({bool(predicate(w)) for w in enumerate_orbit(history)}) <= 1


@dataclass(frozen=True)
class OrbitKripkeFrame:
    """Finite S_s-orbit frame; accessibility is permutation-orbit membership."""
    worlds: frozenset[tuple[int, ...]]

    @classmethod
    def from_history(cls, history: Sequence[int]) -> "OrbitKripkeFrame":
        return cls(enumerate_orbit(history))

    def accessible(self, world: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
        if world not in self.worlds:
            raise ValueError("world is not in this orbit frame")
        return self.worlds

    def satisfies_box(self, world: tuple[int, ...],
                      predicate: Callable[[tuple[int, ...]], bool]) -> bool:
        return all(predicate(v) for v in self.accessible(world))

    def satisfies_diamond(self, world: tuple[int, ...],
                          predicate: Callable[[tuple[int, ...]], bool]) -> bool:
        return any(predicate(v) for v in self.accessible(world))

    def is_equivalence_frame(self) -> bool:
        return bool(self.worlds) and all(self.accessible(w) == self.worlds
                                         for w in self.worlds)


def audit_q6_wilson_microhistory_modal_bridge() -> dict:
    """Small exact audit; no enumeration is attempted at large s."""
    flux = (2, -1, 1, 0)
    data = microhistory_orbit_stabilizer(flux, edge_count=4)
    beta = Fraction(3, 2)
    direct = wilson_factorial_majorant(flux, beta, edge_count=4)
    sym = symmetrized_microhistory_weight(flux, beta, edge_count=4)
    history = labelled_microhistory(flux, edge_count=4)
    frame = OrbitKripkeFrame.from_history(history)
    invariant = lambda w: Counter(w) == Counter(history)
    noninvariant = lambda w: w[0] == history[0]
    return {
        "flux": list(flux),
        "beta": str(beta),
        "orbit_stabilizer": {
            k: (str(v) if isinstance(v, Fraction) else v) for k, v in data.items()
        },
        "factorial_majorant": str(direct),
        "symmetrized_weight": str(sym),
        "symmetrized_identity_exact": direct == sym,
        "orbit_size_matches": len(frame.worlds) == data["orbit_order"],
        "frame_is_equivalence": frame.is_equivalence_frame(),
        "invariant_predicate_is_necessary": frame.satisfies_box(history, invariant),
        "noninvariant_predicate_is_necessary": frame.satisfies_box(history, noninvariant),
        "modal_scope": "S_s-orbit necessity only; not Gödel S5 ontological necessity",
        "checks_passed": (
            direct == sym and len(frame.worlds) == data["orbit_order"]
            and frame.is_equivalence_frame()
            and frame.satisfies_box(history, invariant)
            and not frame.satisfies_box(history, noninvariant)
        ),
    }
