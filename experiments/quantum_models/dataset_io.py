"""Deterministic location resolution for the codon-usage input dataset."""
from __future__ import annotations
import os
from pathlib import Path
from typing import Mapping

def resolve_dataset_path(script_file: str | Path, filename: str = "codon_train.csv", environ: Mapping[str, str] | None = None, cwd: str | Path | None = None) -> Path:
    """Resolve the input CSV without silently depending on the caller CWD.

    CODON_TRAIN_CSV, when set, is authoritative. Otherwise search beside the
    experiment module and then in the current working directory.
    """
    env = os.environ if environ is None else environ
    configured = env.get("CODON_TRAIN_CSV")
    if configured:
        path = Path(configured).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"CODON_TRAIN_CSV points to a missing file: {path}")
        return path.resolve()
    script_dir = Path(script_file).resolve().parent
    work_dir = Path.cwd() if cwd is None else Path(cwd)
    candidates = (script_dir / filename, work_dir / filename)
    for path in candidates:
        if path.is_file():
            return path.resolve()
    searched = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Required dataset {filename!r} was not found. Searched: {searched}. "
        "Place the CSV beside the experiment or set CODON_TRAIN_CSV to its path."
    )
