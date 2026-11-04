#!/usr/bin/env python3
"""Build a conservative static dependency graph for repository experiments.

The report distinguishes Python import edges from literal filename/path mentions.
It does not execute experiments and cannot resolve computed/dynamic paths.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from collections import defaultdict
from pathlib import Path

TEXT_SUFFIXES = {
    ".py", ".ipynb", ".md", ".rst", ".txt", ".yml", ".yaml",
    ".json", ".toml", ".sh", ".cfg", ".ini",
}
SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "node_modules"}


def iter_files(root: Path):
    for path in root.rglob("*"):
        if not path.is_file() or any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix.lower() in TEXT_SUFFIXES:
            yield path


def module_index(root: Path):
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


def resolve_relative(module: str, package: str, level: int):
    package_parts = package.split(".") if package else []
    keep = max(0, len(package_parts) - level + 1)
    prefix = package_parts[:keep]
    if module:
        prefix.extend(module.split("."))
    return ".".join(prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    parser.add_argument("--output", default="reports/EXPERIMENT_DEPENDENCY_GRAPH.json")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    files = list(iter_files(root))
    rel_files = {p.relative_to(root).as_posix(): p for p in files}
    py_files = [p for p in files if p.suffix == ".py"]
    modules = module_index(root)
    edges = []
    parse_errors = []
    file_names = defaultdict(list)
    for rel in rel_files:
        file_names[Path(rel).name].append(rel)

    # Python imports: record only imports resolvable to a repository module.
    for path in py_files:
        rel = path.relative_to(root).as_posix()
        try:
            source = path.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=rel)
        except (UnicodeDecodeError, SyntaxError) as exc:
            parse_errors.append({"path": rel, "error": str(exc)})
            continue
        package = ".".join(Path(rel).with_suffix("").parts[:-1])
        for node in ast.walk(tree):
            candidates = []
            if isinstance(node, ast.Import):
                candidates.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                target = resolve_relative(node.module or "", package, node.level) if node.level else (node.module or "")
                if target:
                    candidates.append(target)
                for alias in node.names:
                    if alias.name != "*":
                        candidates.append(f"{target}.{alias.name}" if target else alias.name)
            for name in candidates:
                if name in modules:
                    edges.append({
                        "source": rel, "target": modules[name],
                        "kind": "python_import", "symbol": name,
                        "line": getattr(node, "lineno", None),
                    })

    # Literal path / filename references in docs, notebooks, workflows and code.
    for path in files:
        rel = path.relative_to(root).as_posix()
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for target in rel_files:
            base = Path(target).name
            if target != rel and (target in text or (len(base) >= 12 and base in text)):
                edges.append({
                    "source": rel, "target": target,
                    "kind": "literal_path_reference" if target in text else "literal_basename_reference",
                })

    # Byte-identical files are reported, never removed.
    duplicate_groups = defaultdict(list)
    for rel, path in rel_files.items():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        duplicate_groups[digest].append(rel)
    duplicates = [
        {"sha256": digest, "paths": sorted(paths)}
        for digest, paths in duplicate_groups.items() if len(paths) > 1
    ]

    unique_edges = {
        (e["source"], e["target"], e["kind"], e.get("symbol"), e.get("line")): e
        for e in edges
    }
    edges = sorted(unique_edges.values(), key=lambda e: (e["source"], e["target"], e["kind"]))
    experiment_nodes = sorted(
        p for p in rel_files if p.startswith("experiments/")
    )
    report = {
        "schema_version": 1,
        "scope": "static repository scan; no experiment execution",
        "counts": {
            "text_files_scanned": len(files),
            "python_files_scanned": len(py_files),
            "experiment_files": len(experiment_nodes),
            "dependency_edges": len(edges),
            "python_import_edges": sum(e["kind"] == "python_import" for e in edges),
            "literal_reference_edges": sum(e["kind"] != "python_import" for e in edges),
            "parse_errors": len(parse_errors),
            "duplicate_groups": len(duplicates),
        },
        "experiment_nodes": experiment_nodes,
        "edges": edges,
        "duplicate_groups": duplicates,
        "parse_errors": parse_errors,
        "limitations": [
            "Computed imports and dynamically constructed paths are not resolved.",
            "Basename matches can be false positives and require human review.",
            "An edge records a textual/static reference, not proof that runtime execution uses it.",
            "No scientific or numerical validity is inferred.",
        ],
    }
    output = (root / args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report["counts"], indent=2))


if __name__ == "__main__":
    main()
