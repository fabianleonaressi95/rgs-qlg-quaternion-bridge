"""Exact bookkeeping for the GEC-6 lower/upper trigram decomposition."""
from itertools import permutations, combinations


LOWER = (0, 1, 2)
UPPER = (3, 4, 5)


def pair_sectors():
    lower = tuple((i, j) for i, j in combinations(LOWER, 2))
    mixed = tuple((i, j) for i in LOWER for j in UPPER)
    upper = tuple((i, j) for i, j in combinations(UPPER, 2))
    return lower, mixed, upper


def block_permutations():
    """Return S3 x S3 together with the lower/upper block exchange."""
    out = set()
    for p in permutations(LOWER):
        for q in permutations(UPPER):
            sigma = [None] * 6
            for i, image in zip(LOWER, p):
                sigma[i] = image
            for i, image in zip(UPPER, q):
                sigma[i] = image
            out.add(tuple(sigma))
            swapped = [None] * 6
            for i, image in zip(LOWER, p):
                swapped[i] = image + 3
            for i, image in zip(UPPER, q):
                swapped[i] = image - 3
            out.add(tuple(swapped))
    return tuple(sorted(out))


def wedge_image(pair, sigma):
    """Image of e_i wedge e_j, retaining the orientation sign."""
    i, j = pair
    a, b = sigma[i], sigma[j]
    return (1, (a, b)) if a < b else (-1, (b, a))


def antipodal_vertex(vertex):
    return vertex ^ 63


def audit():
    lower, mixed, upper = pair_sectors()
    all_pairs = set(lower) | set(mixed) | set(upper)
    assert len(lower) == 3 and len(mixed) == 9 and len(upper) == 3
    assert len(all_pairs) == 15
    group = block_permutations()
    assert len(group) == 72

    lower_set, mixed_set, upper_set = set(lower), set(mixed), set(upper)
    for sigma in group:
        for pair in lower:
            assert wedge_image(pair, sigma)[1] in (lower_set | upper_set)
        for pair in mixed:
            assert wedge_image(pair, sigma)[1] in mixed_set

    # Complementing all bits leaves every coordinate-pair direction unchanged.
    antipodal_checks = 0
    for v in range(64):
        av = antipodal_vertex(v)
        assert antipodal_vertex(av) == v
        for i, j in combinations(range(6), 2):
            assert ((v >> i) & 1) ^ ((av >> i) & 1) == 1
            assert ((v >> j) & 1) ^ ((av >> j) & 1) == 1
            antipodal_checks += 1

    return {
        "pair_sector_ranks": [len(lower), len(mixed), len(upper)],
        "total_pair_directions": len(all_pairs),
        "block_preserving_group_order": 36,
        "block_swap_extended_group_order": len(group),
        "antipodal_vertex_checks": 64,
        "antipodal_coordinate_pair_checks": antipodal_checks,
        "status": "exact finite combinatorial audit; representation ranks proved in documentation",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2, sort_keys=True))
