from src.kathara_15_line_graph import (
    PRIMARY_LINES, audit_kathara_15_line_graph, elementary_edges,
)


def test_source_line_inventory_and_center_count():
    assert len(PRIMARY_LINES) == 15
    assert set(PRIMARY_LINES) == set(range(1, 16))
    audit = audit_kathara_15_line_graph()
    assert len(audit.vertices) == 12


def test_oriented_graph_is_connected_planar_and_acyclic():
    audit = audit_kathara_15_line_graph()
    assert audit.connected
    assert audit.planar_straight_line_embedding
    assert audit.acyclic_orientation
    assert len(audit.directed_edges) == 21


def test_face_count_and_euler_formula():
    audit = audit_kathara_15_line_graph()
    assert len(audit.faces) == 11
    assert len(audit.bounded_faces) == 10
    assert audit.euler_characteristic == 2
    assert len(audit.vertices) - len(audit.directed_edges) + len(audit.faces) == 2


def test_expected_bounded_face_boundaries():
    audit = audit_kathara_15_line_graph()
    got = {frozenset(face) for face in audit.bounded_faces}
    expected = {
        frozenset((1, 2, 3)),
        frozenset((1, 2, 4)),
        frozenset((2, 4, 7, 5)),
        frozenset((2, 5, 6, 3)),
        frozenset((5, 7, 10, 8)),
        frozenset((5, 8, 9, 6)),
        frozenset((8, 10, 11)),
        frozenset((8, 11, 9)),
        frozenset((9, 11, 12)),
        frozenset((10, 12, 11)),
    }
    assert got == expected
