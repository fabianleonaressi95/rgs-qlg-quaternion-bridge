from scripts.audit_oriented_three_channel_incidence import (
    Incidence, TARGET, audit, automorphisms, charpoly_coefficients_3x3,
    graph_isomorphisms, matrix_from_rows,
)

def target_rows():
    return [
        Incidence("a0", "A", "a1", "A", "+"),
        Incidence("a0", "A", "b0", "B", "-"),
        Incidence("b0", "B", "c0", "C", "+"),
        Incidence("c0", "C", "a2", "A", "-"),
    ]

def test_target_matrix_and_integer_characteristic_polynomial():
    m = matrix_from_rows(target_rows())
    assert m == TARGET
    assert charpoly_coefficients_3x3(m) == (1, -1, 0, -1)

def test_target_type_graph_has_only_identity_automorphism():
    assert len(automorphisms(TARGET)) == 1

def test_target_is_recognized_up_to_relabelling():
    assert graph_isomorphisms(TARGET) == [{"A": "A", "B": "B", "C": "C"}]

def test_relabelled_copy_is_recognized():
    # Observed labels A,B,C represent target C,A,B respectively.
    target_to_observed = {"C": "A", "A": "B", "B": "C"}
    changed = [Incidence(r.parent_id, target_to_observed[r.parent_type],
                         r.child_id, target_to_observed[r.child_type], r.orientation)
               for r in target_rows()]
    assert graph_isomorphisms(matrix_from_rows(changed))

def test_antitwin_class_preserving_action_passes():
    rows = target_rows()
    anti = {x: x + "bar" for x in ("a0", "a1", "a2", "b0", "c0")}
    anti.update({v: k for k, v in list(anti.items())})
    doubled = rows + [Incidence(anti[r.parent_id], r.parent_type,
                                anti[r.child_id], r.child_type, r.orientation)
                      for r in rows]
    result = audit(doubled, anti)["anti_twin"]
    assert result["valid_fixed_point_free_involution"]
    assert result["class_preserving"]
    assert result["incidence_invariant"]

def test_antitwin_class_exchange_is_rejected():
    rows = target_rows()
    anti = {"a0": "b0", "b0": "a0", "a1": "a2", "a2": "a1",
            "c0": "cbar", "cbar": "c0"}
    result = audit(rows, anti)["anti_twin"]
    assert not result["class_preserving"]
