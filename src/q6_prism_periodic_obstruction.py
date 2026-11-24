"""No-go theorem for periodic winding from edge-preserving Q6 maps.

The PRISM product is C64 x C6 x C6 x C6 x C5. A graph homomorphism
maps each Q6 edge to one PRISM edge. Every Q6 graph cycle is an integral
sum of its square plaquette boundaries. A square maps to a closed
four-step walk. Since every target factor has length >= 5, such a walk
has zero winding in every factor. Hence the induced map from the Q6
graph cycle lattice to H1(PRISM) = Z^5 is zero.
"""
from __future__ import annotations
from itertools import product
from typing import Iterable, Optional

PRISM_CYCLE_LENGTHS = (64, 6, 6, 6, 5)
Step = tuple[int, int]


def winding_vector(steps: Iterable[Step]) -> Optional[tuple[int, ...]]:
    """Return integer winding vector for a closed walk, else None.

    A step (i, sign) moves one unit in factor i, sign in {-1,+1}.
    """
    displacement = [0] * len(PRISM_CYCLE_LENGTHS)
    for factor, sign in steps:
        if factor not in range(5) or sign not in (-1, 1):
            raise ValueError("step must be (factor 0..4, sign +/-1)")
        displacement[factor] += sign
    if any(displacement[i] % PRISM_CYCLE_LENGTHS[i] for i in range(5)):
        return None
    return tuple(displacement[i] // PRISM_CYCLE_LENGTHS[i] for i in range(5))


def all_closed_walks(length: int = 4):
    """Enumerate factor-direction/sign words that close in the PRISM."""
    elementary_steps = tuple((i, s) for i in range(5) for s in (-1, 1))
    for walk in product(elementary_steps, repeat=length):
        winding = winding_vector(walk)
        if winding is not None:
            yield walk, winding


def audit() -> dict:
    closed = list(all_closed_walks(4))
    nonzero = [(walk, wind) for walk, wind in closed if any(wind)]
    return {
        "prism_cycle_lengths": list(PRISM_CYCLE_LENGTHS),
        "minimum_factor_cycle_length": min(PRISM_CYCLE_LENGTHS),
        "q6_vertices": 64,
        "q6_edges": 192,
        "q6_square_plaquettes": 240,
        "q6_graph_cycle_rank": 192 - 64 + 1,
        "closed_four_step_target_walk_words": len(closed),
        "four_step_closed_walks_with_nonzero_winding": len(nonzero),
        "maximum_possible_absolute_factor_displacement_in_four_steps": 4,
        "minimum_steps_needed_for_nonzero_winding": min(PRISM_CYCLE_LENGTHS),
        "theorem": "every edge-to-edge graph homomorphism Q6 -> PRISM induces zero map from the Q6 graph cycle lattice to H1(PRISM)=Z^5",
        "scope": "periodic homology of the full cubical PRISM; graph cycle space of the PRISM 1-skeleton is a different, larger object",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2))
