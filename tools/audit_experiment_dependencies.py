#!/usr/bin/env python3
"""Build a conservative static dependency graph for repository experiments."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from collections import defaultdict
from pathlib import Path

TEXT_SUFFIXES = {".py", ".ipynb", ".md", ".rst", ".txt", ".yml", ".yaml", ".json", ".toml", ".sh", ".cfg", ".ini"}
SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "node_modules"}


def iter_files(root):
    for path in root.rglob("*"):
        if path.is_file() and not any(part in SKIP_DIRS for part in path.parts) and path.suffix.lower() in TEXT_SUFFIXES:
            yield path


def module_index(root):
    result = {}
    for path in root.rglob("*.py"):
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        rel = path.relative_to(root).with_suffix("")
        parts = list(rel.parts)
        if parts[-1] == "__init__":
            parts.pop()
        if parts:
            result[".".join(parts)] = path.relative_to(root).as_posix()
    return result


def resolve_relative(module, package, level):
    parts = package.split(".") if package else []
    prefix = parts[:max(0, len(parts) - level + 1)]
    if module:
        prefix.extend(module.split("."))
    return ".".join(prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="reports/EXPERIMENT_DEPENDENCY_GRAPH.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    files = list(iter_files(root))
    rel_files = {p.relative_to(root).as_posix(): p for p in files}
    py_files = [p for p in files if p.suffix.lower() == ".py"]
    modules = module_index(root)
    edges, parse_errors = [], []

    for path in py_files:
        rel = path.relative_to(root).as_posix()
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=rel)
        except (UnicodeDecodeError, SyntaxError) as exc:
            parse_errors.append({"path": rel, "error": str(exc)})
            continue
        package = ".".join(Path(rel).with_suffix("").parts[:-1])
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                target = resolve_relative(node.module or "", package, node.level) if node.level else (node.module or "")
                if target:
                    names.append(target)
                names.extend(f"{target}.{alias.name}" for alias in node.names if alias.name != "*" and target)
            for name in names:
                if name in modules:
                    edges.append({"source": rel, "target": modules[name], "kind": "python_import", "symbol": name, "line": getattr(node, "lineno", None)})

    target_paths = list(rel_files)
    target_basenames = {Path(p).name: p for p in target_paths if len(Path(p).name) >= 12}
    for path in files:
        rel = path.relative_to(root).as_posix()
        try:
            content = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for target in target_paths:
            if target != rel and target in content:
                edges.append({"source": rel, "target": target, "kind": "literal_path_reference"})
        for basename, target in target_basenames.items():
            if target != rel and basename in content and target not in content:
                edges.append({"source": rel, "target": target, "kind": "literal_basename_reference"})

    duplicate_groups = defaultdict(list)
    for rel, path in rel_files.items():
        duplicate_groups[hashlib.sha256(path.read_bytes()).hexdigest()].append(rel)
    duplicates = [{"sha256": sha, "paths": sorted(paths)} for sha, paths in duplicate_groups.items() if len(paths) > 1]
    unique = {(e["source"], e["target"], e["kind"], e.get("symbol"), e.get("line")): e for e in edges}
    edges = sorted(unique.values(), key=lambda e: (e["source"], e["target"], e["kind"]))
    experiment_nodes = sorted(p for p in rel_files if p.startswith("experiments/"))
    report = {
        "schema_version": 1,
        "scope": "static repository scan; no experiment execution",
        "counts": {
            "text_files_scanned": len(files), "python_files_scanned": len(py_files),
            "experiment_files": len(experiment_nodes), "dependency_edges": len(edges),
            "python_import_edges": sum(e["kind"] == "python_import" for e in edges),
            "literal_reference_edges": sum(e["kind"] != "python_import" for e in edges),
            "parse_errors": len(parse_errors), "duplicate_groups": len(duplicates),
        },
        "experiment_nodes": experiment_nodes, "edges": edges,
        "duplicate_groups": duplicates, "parse_errors": parse_errors,
        "limitations": [
            "Computed imports and dynamically constructed paths are not resolved.",
            "Basename matches can be false positives and require human review.",
            "Static edges do not prove runtime use; external references are outside this checkout.",
            "No scientific or numerical validity is inferred.",
        ],
    }
    output = (root / args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report["counts"], indent=2))


if __name__ == "__main__":
    main()
