from __future__ import annotations

import json
from pathlib import Path

from src.kathara_15_line_graph import audit_kathara_15_line_graph


def main() -> None:
    a = audit_kathara_15_line_graph()
    report = {
        "vertices": list(a.vertices),
        "primary_lines": {str(k): list(v) for k, v in a.primary_lines.items()},
        "directed_edges": [list(e) for e in a.directed_edges],
        "faces": [list(f) for f in a.faces],
        "bounded_faces": [list(f) for f in a.bounded_faces],
        "euler_characteristic": a.euler_characteristic,
        "connected": a.connected,
        "planar_straight_line_embedding": a.planar_straight_line_embedding,
        "acyclic_orientation": a.acyclic_orientation,
        "scope": "Combinatorial transcription of the published Kathara diagram; not a physical model.",
    }
    out = Path("results/kathara_15_line_graph/audit.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
