"""Run the finite weighted combinatorial Hodge audit for Q6."""
from __future__ import annotations
import json
from src.q6_hodge_laplacian import audit_q6_hodge_laplacians, build_q6_hodge_laplacians


def main() -> None:
    report = audit_q6_hodge_laplacians(build_q6_hodge_laplacians())
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("Q6 Hodge Laplacian audit failed")


if __name__ == "__main__":
    main()
