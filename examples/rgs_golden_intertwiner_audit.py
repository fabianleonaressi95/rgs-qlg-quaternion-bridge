"""Run the RGS/golden cut-and-project affine no-go audit."""
from __future__ import annotations

import json
from pathlib import Path

from src.rgs_psu_intertwiner import export_intertwiner_audit


if __name__ == "__main__":
    report = export_intertwiner_audit(Path("results/rgs_golden_intertwiner.json"))
    print(json.dumps(report, indent=2))
