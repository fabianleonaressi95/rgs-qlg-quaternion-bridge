#!/usr/bin/env python3
"""Run the exact Q6 Wilson microhistory/modal-bridge audit."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rgs_q6_wilson_microhistory_modal_bridge import (
    audit_q6_wilson_microhistory_modal_bridge,
)


def main() -> None:
    report = audit_q6_wilson_microhistory_modal_bridge()
    report["scope"] = (
        "exact finite combinatorics and local orbit-frame semantics; "
        "not an evaluation of the full Q6-D5 fiber or Gödel's ontological proof"
    )
    out = Path("reports/q6_wilson_microhistory_modal_bridge.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    if not report["checks_passed"]:
        raise SystemExit("Q6 Wilson microhistory/modal bridge audit failed")


if __name__ == "__main__":
    main()
