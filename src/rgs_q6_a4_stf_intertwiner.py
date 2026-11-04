"""A4-equivariant linear identification of the Q6 V5 sector with STF(3).

The source is the certified top coexact five-mode A5 module. Restrict it to
a deterministically selected A4 subgroup of the icosahedral rotations and
intertwine it with the tetrahedral action on Sym_0(R^3).

This is a representation intertwiner, not a metric embedding of Q6 or a
physical identification of tensor components.
"""
from __future__ import annotations

from itertools import product
import numpy as np

from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex, induced_cochain_actions


def _group_tables(rotations, tolerance=1e-8):
    n = len(rotations)
    identity = int(np.argmin(np.linalg.norm(rotations - np.eye(3), axis=(1, 2))))
    mul = np.empty((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            d = np.linalg.norm(rotations - rotations[i] @ rotations[j], axis=(1, 2))
            k = int(np.argmin(d))
            if d[k] > tolerance:
                raise ArithmeticError("icosahedral rotations are not closed")
            mul[i, j] = k
    inv = np.empty(n, dtype=int)
    for i in range(n):
        candidates = np.flatnonzero((mul[i] == identity) & (mul[:, i] == identity))
        if len(candidates) != 1:
            raise ArithmeticError("could not identify unique inverse")
        inv[i] = int(candidates[0])
    orders = np.zeros(n, dtype=int)
    for i in range(n):
        x = identity
        for k in range(1, 61):
            x = mul[x, i]
            if x == identity:
                orders[i] = k
                break
        if orders[i] == 0:
            raise ArithmeticError("element order exceeds 60")
    return identity, mul, inv, orders


def _find_a4_subgroup(rotations):
    identity, mul, inv, orders = _group_tables(rotations)
    involutions = [i for i in range(60) if orders[i] == 2]
    order_three = [i for i in range(60) if orders[i] == 3]
    candidates = set()
    for a in involutions:
        for b in involutions:
            if b <= a or mul[a, b] != mul[b, a]:
                continue
            c = mul[a, b]
            if orders[c] != 2:
                continue
            klein = {identity, a, b, c}
            for g in order_three:
                if {mul[mul[g, h], inv[g]] for h in klein} != klein:
                    continue
                subgroup, frontier = set(klein), [g]
                while frontier:
                    x = frontier.pop()
                    if x in subgroup:
                        continue
                    subgroup.add(x)
                    for y in list(subgroup):
                        for z in (mul[x, y], mul[y, x]):
                            if z not in subgroup:
                                frontier.append(z)
                if len(subgroup) == 12:
                    candidates.add(tuple(sorted(subgroup)))
    if not candidates:
        raise ArithmeticError("no A4 subgroup found")
    chosen = min(candidates)
    return {
        "indices": np.asarray(chosen, dtype=int),
        "element_orders": orders,
        "number_of_a4_subgroups_found": len(candidates),
    }


def _stf_basis():
    s2, s6 = np.sqrt(2.0), np.sqrt(6.0)
    return np.asarray([
        np.diag([1., -1., 0.]) / s2,
        np.diag([1., 1., -2.]) / s6,
        np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]]) / s2,
        np.array([[0., 0., 0.], [0., 0., 1.], [0., 1., 0.]]) / s2,
        np.array([[0., 0., 1.], [0., 0., 0.], [1., 0., 0.]]) / s2,
    ])


def _stf_representation(rotation, basis):
    transformed = [rotation @ x @ rotation.T for x in basis]
    return np.asarray([[np.sum(basis[i] * transformed[j]) for j in range(5)]
                       for i in range(5)], dtype=float)


def _constraint_matrix(source_reps, target_reps):
    columns = []
    for i in range(5):
        for j in range(5):
            e = np.zeros((5, 5))
            e[i, j] = 1.0
            columns.append(e)
    return np.vstack([
        np.column_stack([(target @ e - e @ source).ravel() for e in columns])
        for source, target in zip(source_reps, target_reps)
    ])


def build_q6_a4_stf_intertwiner(tolerance=1e-8):
    """Build orthogonal T satisfying rho_STF(g) T = T rho_Q6(g), all g in A4."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    phase = build_q6_d5_phase_space_sector()
    if not phase["all_checks_passed"]:
        raise ArithmeticError("source Q6 five-mode sector failed its audit")
    group = icosahedral_signed_axis_representation()
    rotations = np.asarray(group["rotations"], dtype=float)
    signed = np.asarray(group["signed_permutations"])
    subgroup = _find_a4_subgroup(rotations)
    indices = subgroup["indices"]
    complex_data = build_q6_cochain_complex()
    mode_basis = np.asarray(phase["selected_mode_basis"], dtype=float)
    stf = _stf_basis()
    source_reps, target_reps = [], []
    for index in indices:
        _p0, p1, _p2 = induced_cochain_actions(signed[index], complex_data)
        source_reps.append(mode_basis.T @ p1 @ mode_basis)
        target_reps.append(_stf_representation(rotations[index], stf))

    constraints = _constraint_matrix(source_reps, target_reps)
    _u, singular_values, vh = np.linalg.svd(constraints, full_matrices=True)
    threshold = tolerance * max(1.0, float(singular_values[0]))
    nullity = int(np.count_nonzero(singular_values < threshold))
    if nullity != 3:
        raise ArithmeticError(f"expected Hom_A4 nullity 3, got {nullity}; tail={singular_values[-5:]}")
    null_basis = vh[-nullity:].reshape(nullity, 5, 5)
    candidate = None
    coefficients_used = None
    for coefficients in product(range(-2, 3), repeat=nullity):
        if not any(coefficients):
            continue
        trial = sum(coefficients[k] * null_basis[k] for k in range(nullity))
        if np.linalg.matrix_rank(trial, tol=tolerance) == 5:
            candidate, coefficients_used = trial, coefficients
            break
    if candidate is None:
        raise ArithmeticError("no invertible element found in intertwiner space")
    u, _s, vh_t = np.linalg.svd(candidate, full_matrices=False)
    intertwiner = u @ vh_t
    nz = np.flatnonzero(np.abs(intertwiner.ravel()) > tolerance)
    if len(nz) and intertwiner.ravel()[nz[0]] < 0:
        intertwiner = -intertwiner

    residual = max(float(np.linalg.norm(t @ intertwiner - intertwiner @ s))
                   for s, t in zip(source_reps, target_reps))
    orthogonality = float(np.linalg.norm(intertwiner.T @ intertwiner - np.eye(5)))
    source_orth = max(float(np.linalg.norm(s.T @ s - np.eye(5))) for s in source_reps)
    target_orth = max(float(np.linalg.norm(t.T @ t - np.eye(5))) for t in target_reps)
    orders = subgroup["element_orders"][indices]
    char = {str(o): float(np.mean([np.trace(source_reps[k]) for k in range(12)
                                  if orders[k] == o])) for o in (1, 2, 3)}
    checks = {
        "A4_order_12": len(indices) == 12,
        "source_sector_passed": bool(phase["all_checks_passed"]),
        "Hom_A4_dimension_three": nullity == 3,
        "intertwiner_rank_five": int(np.linalg.matrix_rank(intertwiner)) == 5,
        "equivariance": residual < 20 * tolerance,
        "orthogonal_intertwiner": orthogonality < 20 * tolerance,
        "source_actions_orthogonal": source_orth < 20 * tolerance,
        "target_actions_orthogonal": target_orth < 20 * tolerance,
        "restricted_character": (abs(char["1"] - 5) < 20 * tolerance
                                  and abs(char["2"] - 1) < 20 * tolerance
                                  and abs(char["3"] + 1) < 20 * tolerance),
    }
    return {
        "audit": "Q6_A4_STF3_LINEAR_INTERTWINER",
        "source": "Q6 top coexact lambda=12 five-mode A5 module restricted to A4",
        "target": "Sym_0(R3) under X -> R X R^T",
        "a4_order": len(indices),
        "a4_subgroup_indices": indices,
        "a4_subgroups_found": subgroup["number_of_a4_subgroups_found"],
        "source_character_by_element_order": char,
        "source_representations": np.asarray(source_reps),
        "target_representations": np.asarray(target_reps),
        "stf_basis": stf,
        "intertwiner": intertwiner,
        "intertwiner_nullity": nullity,
        "nullspace_coefficients": coefficients_used,
        "singular_values": singular_values,
        "equivariance_residual": residual,
        "orthogonality_residual": orthogonality,
        "source_orthogonality_residual": source_orth,
        "target_orthogonality_residual": target_orth,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "scope": {
            "established": [
                "A4 restriction has character (5, -1, -1, 1)",
                "STF(3) tetrahedral action has the same A4 representation type",
                "an A4-equivariant linear isomorphism exists",
            ],
            "implementation_audits": [
                "subgroup extraction from the 60 icosahedral rotations",
                "three-dimensional intertwiner nullspace",
                "rank-five orthogonal polar factor",
                "equivariance on all twelve A4 elements",
            ],
            "not_implied": [
                "preferred physical identification of tensor components",
                "global Q6-to-continuum embedding",
                "intrinsic quaternionic structure",
                "Yang-Mills dynamics or a derivation of alpha",
            ],
        },
    }
