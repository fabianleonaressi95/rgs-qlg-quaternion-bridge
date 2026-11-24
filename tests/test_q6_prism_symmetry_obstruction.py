from src.q6_prism_symmetry_obstruction import (
    A6_ORDER,
    TARGET_H1_AUTOMORPHISM_BOUND,
    audit,
    invariant_cycle_boundary_matrix,
    invariant_cycle_space_dimension,
)


def test_a6_order_exceeds_target_h1_symmetry_image_bound():
    assert A6_ORDER == 360
    assert TARGET_H1_AUTOMORPHISM_BOUND == 192
    assert A6_ORDER > TARGET_H1_AUTOMORPHISM_BOUND


def test_invariant_cycle_recurrence_has_only_zero_solution():
    matrix = invariant_cycle_boundary_matrix()
    assert len(matrix) == 7
    assert all(len(row) == 6 for row in matrix)
    # The rows successively force a_0 through a_5 to vanish.
    assert matrix[0] == (-6, 0, 0, 0, 0, 0)
    assert matrix[1] == (1, -5, 0, 0, 0, 0)
    assert matrix[6] == (0, 0, 0, 0, 0, 6)
    assert invariant_cycle_space_dimension() == 0


def test_audit_excludes_equivariant_surjection():
    result = audit()
    assert result["A6_action_on_target_H1_forced_trivial"]
    assert result["equivariant_induced_H1_map"] == "zero"
    assert not result["equivariant_surjection_to_Z5"]
