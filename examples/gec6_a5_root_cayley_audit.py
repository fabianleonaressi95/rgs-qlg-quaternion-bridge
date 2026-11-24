"""Generate the finite GEC-6 A5-root Cayley graph audit."""
import json
from pathlib import Path

from src.gec6_a5_root_cayley_torus import finite_a5_audit

result = finite_a5_audit(3)
out = Path("results/gec6_a5_root_cayley_audit.json")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, sort_keys=True))
