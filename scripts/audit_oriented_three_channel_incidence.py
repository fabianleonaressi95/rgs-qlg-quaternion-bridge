#!/usr/bin/env python3
"""Audit an oriented three-channel parent->child incidence table.

CSV schemas:
  children: parent_id,parent_type,child_id,child_type,orientation
  anti-twin: id,antipode

The tool audits supplied data; it never infers or constructs geometry.
"""
from __future__ import annotations
import argparse, csv, itertools, json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

TYPES = ("A", "B", "C")
# Rows = child types; columns = parent types.
TARGET = ((1, 0, 1), (1, 0, 0), (0, 1, 0))

@dataclass(frozen=True)
class Incidence:
    parent_id: str
    parent_type: str
    child_id: str
    child_type: str
    orientation: str

def read_children(path: Path) -> list[Incidence]:
    required = {"parent_id", "parent_type", "child_id", "child_type", "orientation"}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = []
        for line, row in enumerate(reader, 2):
            p, pt = row["parent_id"].strip(), row["parent_type"].strip()
            c, ct, sign = row["child_id"].strip(), row["child_type"].strip(), row["orientation"].strip()
            if not p or not c:
                raise ValueError(f"{path}:{line}: empty element identifier")
            if pt not in TYPES or ct not in TYPES:
                raise ValueError(f"{path}:{line}: types must be one of {TYPES}")
            if sign not in {"+", "-"}:
                raise ValueError(f"{path}:{line}: orientation must be '+' or '-'")
            rows.append(Incidence(p, pt, c, ct, sign))
    if not rows:
        raise ValueError(f"{path}: no incidence rows")
    return rows

def read_antipodes(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"id", "antipode"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        result = {}
        for line, row in enumerate(reader, 2):
            a, b = row["id"].strip(), row["antipode"].strip()
            if not a or not b:
                raise ValueError(f"{path}:{line}: empty identifier")
            if a in result:
                raise ValueError(f"{path}:{line}: duplicate id {a}")
            result[a] = b
    if not result:
        raise ValueError(f"{path}: no antipodal pairs")
    return result

def matrix_from_rows(rows: Iterable[Incidence]) -> tuple[tuple[int, ...], ...]:
    counts = Counter((r.child_type, r.parent_type) for r in rows)
    return tuple(tuple(counts[(child, parent)] for parent in TYPES) for child in TYPES)

def charpoly_coefficients_3x3(m):
    """Return descending coefficients of det(lambda I - M), exactly."""
    a, b, c = m[0]; d, e, f = m[1]; g, h, i = m[2]
    trace = a + e + i
    minors2 = (a*e-b*d) + (a*i-c*g) + (e*i-f*h)
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    return (1, -trace, minors2, -det)

def format_poly(coeffs):
    terms = []
    for coefficient, power in zip(coeffs, (3, 2, 1, 0)):
        if not coefficient:
            continue
        atom = ("" if abs(coefficient) == 1 and power else str(abs(coefficient)))
        atom += "lambda" if power == 1 else (f"lambda^{power}" if power else "")
        if not atom:
            atom = str(abs(coefficient))
        terms.append((("-" if coefficient < 0 else "") if not terms else
                      (" - " if coefficient < 0 else " + ")) + atom)
    return "".join(terms) or "0"

def relabelled_matrix(m, perm):
    out = [[0]*3 for _ in range(3)]
    for ci, child in enumerate(TYPES):
        for pi, parent in enumerate(TYPES):
            out[TYPES.index(perm[child])][TYPES.index(perm[parent])] = m[ci][pi]
    return tuple(tuple(row) for row in out)

def graph_isomorphisms(m):
    return [dict(zip(TYPES, p)) for p in itertools.permutations(TYPES)
            if relabelled_matrix(m, dict(zip(TYPES, p))) == TARGET]

def automorphisms(m):
    return [dict(zip(TYPES, p)) for p in itertools.permutations(TYPES)
            if relabelled_matrix(m, dict(zip(TYPES, p))) == m]

def audit_antitwin(rows, antipodes):
    ids = {r.parent_id for r in rows} | {r.child_id for r in rows}
    missing = sorted(ids - set(antipodes))
    unknown = sorted(set(antipodes.values()) - set(antipodes))
    involution_failures = sorted(x for x in antipodes if antipodes.get(antipodes[x]) != x)
    fixed_points = sorted(x for x, y in antipodes.items() if x == y)
    type_of, conflicts = {}, []
    for r in rows:
        for element, kind in ((r.parent_id, r.parent_type), (r.child_id, r.child_type)):
            if element in type_of and type_of[element] != kind:
                conflicts.append(element)
            type_of[element] = kind
    class_changes = []
    for element, kind in type_of.items():
        partner = antipodes.get(element)
        if partner in type_of and type_of[partner] != kind:
            class_changes.append((element, partner, kind, type_of[partner]))
    edges = Counter((r.parent_id, r.child_id) for r in rows)
    anti_edges = Counter()
    for (p, c), count in edges.items():
        anti_edges[(antipodes.get(p, ""), antipodes.get(c, ""))] += count
    oriented = Counter((r.parent_id, r.child_id, r.orientation) for r in rows)
    anti_same, anti_flip = Counter(), Counter()
    for (p, c, sign), count in oriented.items():
        anti_same[(antipodes.get(p, ""), antipodes.get(c, ""), sign)] += count
        anti_flip[(antipodes.get(p, ""), antipodes.get(c, ""),
                   "+" if sign == "-" else "-")] += count
    same, flipped = oriented == anti_same, oriented == anti_flip
    return {
        "complete_on_incidence_ids": not missing, "missing_ids": missing,
        "unknown_antipode_targets": unknown,
        "involution": not involution_failures, "involution_failures": involution_failures,
        "fixed_points": fixed_points, "element_type_consistent": not conflicts,
        "type_conflicts": sorted(set(conflicts)),
        "class_preserving": not class_changes, "class_changes": class_changes,
        "incidence_invariant": edges == anti_edges,
        "orientation_action": "preserved" if same else "reversed" if flipped else "mixed_or_not_invariant",
        "valid_fixed_point_free_involution": (
            not missing and not unknown and not involution_failures and
            not fixed_points and not conflicts),
    }

def audit(rows, antipodes=None):
    duplicates = [list(k) for k, n in Counter(
        (r.parent_id, r.child_id, r.orientation) for r in rows).items() if n > 1]
    matrix = matrix_from_rows(rows)
    polynomial = charpoly_coefficients_3x3(matrix)
    iso, autos = graph_isomorphisms(matrix), automorphisms(matrix)
    result = {
        "types_order": list(TYPES),
        "matrix_rows_child_columns_parent": [list(row) for row in matrix],
        "target_matrix": [list(row) for row in TARGET],
        "target_graph_isomorphic_up_to_relabelling": bool(iso),
        "type_relabellings": iso,
        "graph_automorphisms": autos,
        "automorphism_group_order": len(autos),
        "graph_automorphism_trivial": len(autos) == 1,
        "characteristic_polynomial_coefficients_descending": list(polynomial),
        "characteristic_polynomial": format_poly(polynomial),
        "target_characteristic_polynomial": "lambda^3 - lambda^2 - 1",
        "characteristic_polynomial_matches_target": polynomial == (1, -1, 0, -1),
        "duplicate_oriented_parent_child_rows": duplicates,
        "orientation_counts": dict(sorted(Counter(r.orientation for r in rows).items())),
        "anti_twin": audit_antitwin(rows, antipodes) if antipodes is not None else None,
    }
    result["abstract_substitution_gate_passed"] = (
        bool(iso) and polynomial == (1, -1, 0, -1) and len(autos) == 1 and not duplicates)
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("children_csv", type=Path)
    parser.add_argument("--anti-twin", type=Path, help="CSV with columns id,antipode")
    parser.add_argument("--output", type=Path, help="write JSON report")
    args = parser.parse_args()
    report = audit(read_children(args.children_csv),
                   read_antipodes(args.anti_twin) if args.anti_twin else None)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0 if report["abstract_substitution_gate_passed"] else 2

if __name__ == "__main__":
    raise SystemExit(main())
