#!/usr/bin/env python3
"""Stage-1 filter for level counts in a proposed three-channel refinement.

CSV columns: n,A,B,C. Recurrence checks use exact integer arithmetic;
the supergolden frequency target is computed with Decimal arithmetic.
"""
from __future__ import annotations
import argparse, csv, json
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 60
ONE = Decimal(1)

def supergolden() -> Decimal:
    x = Decimal("1.465571231876768026656731225219939108025577568472")
    for _ in range(30):
        x -= (x*x*x - x*x - ONE) / (Decimal(3)*x*x - Decimal(2)*x)
    return +x

def read_counts(path: Path) -> list[dict]:
    required = ("n", "A", "B", "C")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows, seen = [], set()
        for line, row in enumerate(reader, 2):
            try:
                n, a, b, c = (int(row[k]) for k in required)
            except (TypeError, ValueError):
                raise ValueError(f"{path}:{line}: n,A,B,C must be integers")
            if min(n, a, b, c) < 0:
                raise ValueError(f"{path}:{line}: values must be nonnegative")
            if n in seen:
                raise ValueError(f"{path}:{line}: duplicate level n={n}")
            seen.add(n)
            rows.append({"n": n, "A": a, "B": b, "C": c, "total": a+b+c})
    if not rows:
        raise ValueError(f"{path}: no count rows")
    rows.sort(key=lambda r: r["n"])
    if [r["n"] for r in rows] != list(range(rows[0]["n"], rows[-1]["n"]+1)):
        raise ValueError(f"{path}: level labels must be contiguous")
    return rows

def ds(x: Decimal) -> str:
    return format(x, ".16g")

def analyse(rows: list[dict], tolerance: Decimal | None = None) -> dict:
    s = supergolden()
    targets = {"A": s**-2, "B": s**-3, "C": s**-4}
    if abs(sum(targets.values()) - ONE) > Decimal("1e-50"):
        raise ArithmeticError("internal Perron-frequency normalization failed")

    vector_failures, total_failures = [], []
    tested = max(0, len(rows)-3)
    for i in range(tested):
        r0, r2, r3 = rows[i], rows[i+2], rows[i+3]
        for key in ("A", "B", "C"):
            if r3[key] != r2[key] + r0[key]:
                vector_failures.append({"n": r0["n"], "channel": key,
                    "observed": r3[key], "expected": r2[key]+r0[key]})
        if r3["total"] != r2["total"] + r0["total"]:
            total_failures.append({"n": r0["n"], "observed": r3["total"],
                "expected": r2["total"]+r0["total"]})

    frequencies = []
    for row in rows:
        if row["total"] == 0:
            continue
        obs = {k: Decimal(row[k])/Decimal(row["total"]) for k in ("A","B","C")}
        errors = {k: abs(obs[k]-targets[k]) for k in obs}
        frequencies.append({"n": row["n"],
            "fractions": {k: ds(v) for k,v in obs.items()},
            "absolute_errors": {k: ds(v) for k,v in errors.items()},
            "max_absolute_error": ds(max(errors.values()))})

    last = frequencies[-1] if frequencies else None
    freq_pass = None if tolerance is None else (
        last is not None and Decimal(last["max_absolute_error"]) <= tolerance)
    ratios = [{"from_n": a["n"], "to_n": b["n"],
               "ratio": ds(Decimal(b["total"])/Decimal(a["total"]))}
              for a,b in zip(rows, rows[1:]) if a["total"]]

    return {
        "levels": len(rows), "level_range": [rows[0]["n"], rows[-1]["n"]],
        "supergolden_root": ds(s),
        "target_frequencies": {k: ds(v) for k,v in targets.items()},
        "target_frequency_sum": ds(sum(targets.values())),
        "total_recurrence": {"equation": "N[n+3] = N[n+2] + N[n]",
            "tested": tested, "passed": not total_failures, "failures": total_failures},
        "vector_recurrence": {"equation": "v[n+3] = v[n+2] + v[n]",
            "tested": tested, "passed": not vector_failures, "failures": vector_failures},
        "frequency_diagnostics_by_level": frequencies,
        "last_level_frequency_check": {"result": last,
            "tolerance": ds(tolerance) if tolerance is not None else None,
            "passed": freq_pass},
        "total_growth_ratios": ratios,
        "stage1_exact_recurrence_gate_passed": not vector_failures and not total_failures,
    }

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts_csv", type=Path)
    parser.add_argument("--frequency-tolerance", type=Decimal, default=None)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.frequency_tolerance is not None and args.frequency_tolerance < 0:
        parser.error("--frequency-tolerance must be nonnegative")
    report = analyse(read_counts(args.counts_csv), args.frequency_tolerance)
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0 if report["stage1_exact_recurrence_gate_passed"] else 2

if __name__ == "__main__":
    raise SystemExit(main())
