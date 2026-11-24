"""Write the Q6/PRISM embedding obstruction audit."""
import json
from pathlib import Path
from src.q6_prism_embedding_obstruction import audit

result = audit()
path = Path("results/q6_prism_embedding_audit.json")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
