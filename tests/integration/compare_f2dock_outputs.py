#!/usr/bin/env python3
"""Compare F2Dock output tables with tolerance for numeric fields."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def normalized_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(line.split())
    return rows


def is_number(token: str) -> bool:
    return bool(NUMERIC_RE.match(token))


def compare_rows(expected: list[str], actual: list[str], abs_tol: float, rel_tol: float) -> tuple[bool, str]:
    if len(expected) != len(actual):
        return False, f"field count mismatch (expected {len(expected)}, actual {len(actual)})"

    for idx, (exp_tok, act_tok) in enumerate(zip(expected, actual, strict=True)):
        exp_num = is_number(exp_tok)
        act_num = is_number(act_tok)

        if exp_num and act_num:
            exp_val = float(exp_tok)
            act_val = float(act_tok)
            delta = abs(act_val - exp_val)
            limit = abs_tol + rel_tol * max(abs(exp_val), abs(act_val))
            if math.isnan(delta) or delta > limit:
                return (
                    False,
                    f"numeric mismatch at column {idx}: expected {exp_val}, actual {act_val}, "
                    f"delta {delta} > limit {limit}",
                )
        elif exp_tok != act_tok:
            return False, f"token mismatch at column {idx}: expected '{exp_tok}', actual '{act_tok}'"

    return True, ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", required=True)
    parser.add_argument("--actual", required=True)
    parser.add_argument("--rows", type=int, default=100)
    parser.add_argument("--abs-tol", type=float, default=1e-3)
    parser.add_argument("--rel-tol", type=float, default=1e-5)
    args = parser.parse_args()

    expected_path = Path(args.expected)
    actual_path = Path(args.actual)

    if not expected_path.exists():
        print(f"expected output file missing: {expected_path}")
        return 2
    if not actual_path.exists():
        print(f"actual output file missing: {actual_path}")
        return 2

    expected_rows = normalized_rows(expected_path)
    actual_rows = normalized_rows(actual_path)

    if not expected_rows:
        print(f"expected output has no data rows: {expected_path}")
        return 3
    if not actual_rows:
        print(f"actual output has no data rows: {actual_path}")
        return 3

    compare_count = min(args.rows, len(expected_rows), len(actual_rows))
    if compare_count <= 0:
        print("no rows available for comparison")
        return 4

    for row_idx in range(compare_count):
        ok, reason = compare_rows(expected_rows[row_idx], actual_rows[row_idx], args.abs_tol, args.rel_tol)
        if not ok:
            print(f"row {row_idx} mismatch: {reason}")
            print("EXPECTED:")
            print(" ".join(expected_rows[row_idx]))
            print("ACTUAL:")
            print(" ".join(actual_rows[row_idx]))
            return 5

    print(
        f"compared {compare_count} rows successfully "
        f"(expected rows={len(expected_rows)}, actual rows={len(actual_rows)}, "
        f"abs_tol={args.abs_tol}, rel_tol={args.rel_tol})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
