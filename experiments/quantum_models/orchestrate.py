"""Generate and run the 64-node aperiodic TCP benchmark with Docker Compose.

The node workload only uses ordinary TCP sockets; it does not require
NET_ADMIN or other elevated Linux capabilities.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

TOTAL_NODES = 64
ROOT = Path(__file__).resolve().parent


def build_compose_config(total_nodes: int = TOTAL_NODES) -> dict:
    if total_nodes < 2:
        raise ValueError("total_nodes must be at least 2")

    services = {}
    for node_id in range(total_nodes):
        service_name = f"node_{node_id:02d}"
        services[service_name] = {
            "build": {
                "context": str(ROOT),
                "dockerfile": "Dockerfile",
            },
            "container_name": service_name,
            "environment": {
                "NODE_ID": str(node_id),
                "TOTAL_NODES": str(total_nodes),
            },
            "networks": ["aperiodic_mesh"],
            "init": True,
        }

    return {
        "services": services,
        "networks": {"aperiodic_mesh": {"driver": "bridge"}},
    }


def generate_compose_file(
    output_path: str | Path | None = None,
    total_nodes: int = TOTAL_NODES,
) -> Path:
    """Write Compose-compatible JSON (JSON is a YAML subset) and return its path."""
    config = build_compose_config(total_nodes)
    if output_path is None:
        handle = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".json",
            prefix="rgs-aperiodic-compose-",
            delete=False,
        )
        path = Path(handle.name)
        with handle:
            json.dump(config, handle, indent=2)
            handle.write("\n")
    else:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def run_benchmark(compose_file: str | Path | None = None) -> None:
    """Run Compose and always attempt cleanup, propagating benchmark failures."""
    generated_here = compose_file is None
    compose_path = Path(compose_file) if compose_file is not None else generate_compose_file()
    compose_cmd: Sequence[str] = ("docker", "compose", "--file", str(compose_path))

    try:
        print("[Orchestrator] Avvio container e raccolta metriche di benchmark...")
        result = subprocess.run(
            [*compose_cmd, "up", "--build", "--abort-on-container-exit", "--remove-orphans"],
            capture_output=True,
            text=True,
            check=True,
        )
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("Log di sistema/errori:\n", result.stderr)
        print("[Orchestrator] Benchmark completato.")
    finally:
        print("[Orchestrator] Pulizia risorse Docker...")
        cleanup = subprocess.run(
            [*compose_cmd, "down", "--remove-orphans"],
            capture_output=True,
            text=True,
            check=False,
        )
        if cleanup.returncode != 0:
            print("Avviso: cleanup Docker non riuscito:\n", cleanup.stderr)
        if generated_here:
            compose_path.unlink(missing_ok=True)


if __name__ == "__main__":
    run_benchmark()
