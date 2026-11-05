"""Static contracts for the numbered microtubule coordinate pipelines."""
import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "experiments" / "microtubule_pipeline"


class DensityCoordinateContractTests(unittest.TestCase):
    def test_definition_and_calls_have_matching_arity(self):
        paths = [
            PIPELINE / "mvp_middleware47.py",
            PIPELINE / "mvp_middleware48.py",
            PIPELINE / "mvp_middleware49.py",
            PIPELINE / "mvp_middleware50.py",
            PIPELINE / "qip_qc_phase1_microtubule5.py",
            PIPELINE / "qip_qc_phase1_microtubule6.py",
        ]
        for path in paths:
            with self.subTest(path=path.name):
                tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
                definitions = [
                    node for node in ast.walk(tree)
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and node.name == "density_coordinates"
                ]
                self.assertEqual(len(definitions), 1)
                expected = len(definitions[0].args.posonlyargs) + len(definitions[0].args.args)
                calls = [
                    node for node in ast.walk(tree)
                    if isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "density_coordinates"
                ]
                self.assertTrue(calls, "expected at least one density_coordinates call")
                for call in calls:
                    self.assertEqual(
                        len(call.args), expected,
                        f"{path.name}:{call.lineno} passes {len(call.args)} positional "
                        f"arguments; definition expects {expected}",
                    )


if __name__ == "__main__":
    unittest.main()
