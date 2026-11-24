"""Print the Q3/Q6-to-network finite bridge audit as JSON."""
import json
from pathlib import Path
from src.q3_q6_network_bridge import audit

result = audit()
output = Path("results/q3_q6_network_prism_audit.json")
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
