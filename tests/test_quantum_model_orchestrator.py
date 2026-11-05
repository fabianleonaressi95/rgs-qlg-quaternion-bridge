from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, call

from experiments.quantum_models import orchestrate


class OrchestratorTests(unittest.TestCase):
    def test_compose_has_64_nodes_and_unprivileged_networking(self):
        config = orchestrate.build_compose_config()
        self.assertEqual(len(config["services"]), 64)
        self.assertIn("node_00", config["services"])
        self.assertIn("node_63", config["services"])
        for node_id, service in enumerate(config["services"].values()):
            self.assertEqual(service["environment"]["NODE_ID"], str(node_id))
            self.assertEqual(service["environment"]["TOTAL_NODES"], "64")
            self.assertEqual(service["build"]["dockerfile"], "Dockerfile")
            self.assertNotIn("cap_add", service)

    def test_compose_file_is_valid_json_and_uses_requested_size(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "compose.json"
            returned = orchestrate.generate_compose_file(path, total_nodes=5)
            self.assertEqual(returned, path)
            config = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(len(config["services"]), 5)

    def test_total_nodes_must_be_at_least_two(self):
        with self.assertRaises(ValueError):
            orchestrate.build_compose_config(1)

    @patch("experiments.quantum_models.orchestrate.subprocess.run")
    def test_cleanup_runs_even_if_benchmark_fails(self, run):
        run.side_effect = [
            __import__("subprocess").CalledProcessError(1, ["docker"], stderr="failed"),
            __import__("subprocess").CompletedProcess(["docker"], 0, "", ""),
        ]
        with tempfile.TemporaryDirectory() as folder:
            compose = Path(folder) / "compose.json"
            compose.write_text("{}", encoding="utf-8")
            with self.assertRaises(__import__("subprocess").CalledProcessError):
                orchestrate.run_benchmark(compose)
        self.assertEqual(run.call_count, 2)
        self.assertIn("up", run.call_args_list[0].args[0])
        self.assertIn("down", run.call_args_list[1].args[0])


if __name__ == "__main__":
    unittest.main()
