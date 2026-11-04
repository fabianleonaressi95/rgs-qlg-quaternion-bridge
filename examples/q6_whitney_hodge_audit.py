"""Reproducible geometric Whitney-Hodge audit for the finite Q6 complex."""
import json
from src.q6_whitney_hodge import audit_q6_whitney_hodge


if __name__ == "__main__":
    print(json.dumps(audit_q6_whitney_hodge(), indent=2, sort_keys=True))
