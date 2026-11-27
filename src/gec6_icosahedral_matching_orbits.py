"""Finite A5 orbit audit for matchings of the six icosahedral axes.

Orbit calculations use the permutation projection on unoriented axis labels.
The full signed Q6 action and its compatibility with the GEC-6 contraction
are audited separately; these actions must not be conflated.
"""
from __future__ import annotations

from itertools import combinations, permutations
from collections import Counter
from functools import lru_cache
import numpy as np

from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.rgs_weyl_phase_closure import build_icosahedral_c6

LOWER = (0, 1, 2)
UPPER = (3, 4, 5)
DIAGONAL_MATCHING = ((0, 3), (1, 4), (2, 5))
BLOCK_PARTITION = (0, 1, 2)


@lru_cache(maxsize=1)
def _signed_actions():
    data = icosahedral_signed_axis_representation()
    actions = []
    for matrix in np.asarray(data["signed_permutations"]):
        permutation = tuple(int(np.flatnonzero(matrix[:, j])[0]) for j in range(6))
        signs = tuple(int(matrix[permutation[j], j]) for j in range(6))
        actions.append((permutation, signs))
    if len(set(actions)) != 60:
        raise ArithmeticError("expected 60 distinct signed-axis actions")
    return tuple(actions)


def _canonical_partition(part):
    part = set(part)
    other = set(range(6)) - part
    return tuple(sorted(part if 0 in part else other))


def _all_partitions():
    return tuple((0,) + pair for pair in combinations(range(1, 6), 2))


def _all_matchings(items=(0, 1, 2, 3, 4, 5)):
    items = tuple(items)
    if not items:
        yield ()
        return
    first = items[0]
    for k in range(1, len(items)):
        second = items[k]
        rest = items[1:k] + items[k + 1:]
        for tail in _all_matchings(rest):
            yield tuple(sorted(((first, second),) + tail))


def _image_pair(pair, permutation):
    return tuple(sorted((permutation[pair[0]], permutation[pair[1]])))


def _image_partition(partition, permutation):
    return _canonical_partition(permutation[i] for i in partition)


def _image_matching(matching, permutation):
    return tuple(sorted(_image_pair(pair, permutation) for pair in matching))


def _orbit_decomposition(items, permutations_only, image):
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {image(seed, p) for p in permutations_only}
        orbits.append(tuple(sorted(orbit)))
        unseen.difference_update(orbit)
    return tuple(sorted(orbits, key=lambda orbit: (len(orbit), orbit)))


def _permutation_order(permutation):
    identity = tuple(range(6))
    current = identity
    for order in range(1, 61):
        current = tuple(permutation[current[i]] for i in range(6))
        if current == identity:
            return order
    raise ArithmeticError("permutation order exceeded the A5 group order")


def _cross_matchings(partition):
    lower = tuple(partition)
    upper = tuple(i for i in range(6) if i not in set(lower))
    return tuple(
        tuple(sorted(tuple(sorted(pair)) for pair in zip(lower, target)))
        for target in permutations(upper)
    )


def audit_icosahedral_matching_orbits():
    actions = _signed_actions()
    permutations_only = tuple(action[0] for action in actions)
    pairs = tuple(combinations(range(6), 2))
    partitions = _all_partitions()
    matchings = tuple(_all_matchings())

    if len(set(permutations_only)) != 60:
        raise ArithmeticError("permutation projection is not faithful on the six axes")
    pair_orbits = _orbit_decomposition(pairs, permutations_only, _image_pair)
    partition_orbits = _orbit_decomposition(
        partitions, permutations_only, _image_partition
    )
    matching_orbits = _orbit_decomposition(
        matchings, permutations_only, _image_matching
    )

    c6 = np.asarray(build_icosahedral_c6()["C6"], dtype=float)
    c6_integer = np.rint(c6).astype(int)
    if np.max(np.abs(c6 - c6_integer)) > 1e-8:
        raise ArithmeticError("C6 entries are not integral in this axis convention")

    matching_products = {
        matching: int(np.prod([c6_integer[i, j] for i, j in matching]))
        for matching in matchings
    }
    special = {m for m, product in matching_products.items() if product == 1}
    generic = {m for m, product in matching_products.items() if product == -1}

    lower = set(LOWER)
    block_stabilizer = {
        p for p in permutations_only
        if set(p[i] for i in lower) in (lower, set(UPPER))
    }
    quotient_compatible = {
        p for p, signs in actions if len(set(signs)) == 1
    }
    diagonal = tuple(sorted(tuple(sorted(pair)) for pair in DIAGONAL_MATCHING))
    diagonal_orbit = {
        _image_matching(diagonal, p) for p in permutations_only
    }
    special_cross_by_partition = {
        partition: sorted(set(_cross_matchings(partition)) & special)
        for partition in partitions
    }
    special_orbit = set(matching_orbits[0])
    generic_orbit = set(matching_orbits[1])
    special_stabilizers = []
    special_stabilizer_profiles = []
    for matching in sorted(special_orbit):
        stabilizer = {
            p for p in permutations_only
            if _image_matching(matching, p) == matching
        }
        profile = Counter(_permutation_order(p) for p in stabilizer)
        special_stabilizers.append(frozenset(stabilizer))
        special_stabilizer_profiles.append(
            {str(order): profile[order] for order in sorted(profile)}
        )

    generic_stabilizer_orders = {
        sum(_image_matching(matching, p) == matching for p in permutations_only)
        for matching in generic_orbit
    }

    return {
        "group_order": len(actions),
        "pair_orbit_sizes": list(map(len, pair_orbits)),
        "three_plus_three_partition_orbit_sizes": list(map(len, partition_orbits)),
        "perfect_matching_orbit_sizes": list(map(len, matching_orbits)),
        "special_matching_count_C6_product_plus": len(special),
        "generic_matching_count_C6_product_minus": len(generic),
        "C6_product_classes_equal_A5_matching_orbits": (
            special == special_orbit and generic == generic_orbit
        ),
        "diagonal_matching_in_special_orbit": diagonal in special,
        "diagonal_matching_orbit_size": len(diagonal_orbit),
        "special_matching_stabilizer_orders": sorted({
            len(stabilizer) for stabilizer in special_stabilizers
        }),
        "distinct_special_matching_stabilizers": len(set(special_stabilizers)),
        "special_matching_stabilizer_element_order_profiles": (
            special_stabilizer_profiles
        ),
        "generic_matching_stabilizer_orders": sorted(generic_stabilizer_orders),
        "special_cross_matchings_per_three_plus_three_partition": sorted(
            {len(value) for value in special_cross_by_partition.values()}
        ),
        "special_cross_matchings_for_lower_upper_partition": (
            special_cross_by_partition[BLOCK_PARTITION]
        ),
        "block_partition_stabilizer_order": len(block_stabilizer),
        "GEC6_contraction_kernel_preserving_signed_actions": len(quotient_compatible),
        "quotient_preserving_group_equals_block_partition_stabilizer": (
            quotient_compatible == block_stabilizer
        ),
        "diagonal_matching_orbit_under_quotient_preserving_subgroup": len({
            _image_matching(diagonal, p) for p in quotient_compatible
        }),
        "signed_permutation_sign_products": sorted({
            int(np.prod(signs)) for _, signs in actions
        }),
        "status": (
            "exact orbit enumeration conditional on the numerically reconstructed "
            "60 signed-axis matrices; GEC-6 kernel criterion is exact"
        ),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit_icosahedral_matching_orbits(), indent=2, sort_keys=True))
