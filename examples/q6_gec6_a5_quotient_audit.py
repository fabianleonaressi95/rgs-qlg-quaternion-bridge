"""Run the exact GEC-6 square-curl and S6 root-action audit."""
import json
from pathlib import Path

from src.q6_gec6_a5_quotient import audit

result = audit()
out = Path("results/q6_gec6_a5_quotient_audit.json")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, sort_keys=True))
