"""Reproducible audit of the Euclidean cubical Whitney Hodge construction."""
import json
from src.q6_geometric_hodge import audit_cubical_whitney_hodge

if __name__ == "__main__":
    print(json.dumps(audit_cubical_whitney_hodge(), indent=2, sort_keys=True))
