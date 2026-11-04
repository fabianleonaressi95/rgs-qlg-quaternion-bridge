"""Export the exact RGS–Q6 base/sign-fibre metric audit."""
from pathlib import Path
import json

from src.rgs_psu_fibred_intertwiner import export_fibred_intertwiner_audit


if __name__ == "__main__":
    report = export_fibred_intertwiner_audit(
        Path("results/rgs_psu_fibred_intertwiner")
    )
    print(json.dumps(report, indent=2))
