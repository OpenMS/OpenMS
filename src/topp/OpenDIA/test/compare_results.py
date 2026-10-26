#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Justin Sing $
"""Compare direct and tiled OpenDIA exports and all feature/score tables."""
import csv
import math
from pathlib import Path
import sqlite3
import sys


IGNORED_EXPORT = {"run_id", "filename", "run_name", "feature_id"}


def same_value(left, right):
    if left == right:
        return True
    if isinstance(left, (tuple, list)) or isinstance(right, (tuple, list)):
        return (isinstance(left, (tuple, list)) and isinstance(right, (tuple, list))
                and len(left) == len(right)
                and all(same_value(a, b) for a, b in zip(left, right)))
    if left is None or right is None:
        return False
    if isinstance(left, int) and isinstance(right, int):
        return False
    if isinstance(left, str) and isinstance(right, str) and (";" in left or ";" in right):
        return same_value(left.split(";"), right.split(";"))
    try:
        a, b = float(left), float(right)
    except (ValueError, TypeError):
        return False
    return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-8)


def compare_rows(label, columns, left, right):
    if len(left) != len(right) or not left:
        raise ValueError(f"{label}: row count differs or is zero: {len(left)}, {len(right)}")
    for index, (a, b) in enumerate(zip(left, right)):
        for column, x, y in zip(columns, a, b):
            exact = column.upper() == "ID" or column.upper().endswith("_ID") or column in {"peak_group_rank", "decoy"}
            matches = x == y if exact else same_value(x, y)
            if not matches:
                raise ValueError(f"{label} row {index} column {column}: {x!r} != {y!r}")
    print(f"PASS {label}: {len(left)} rows, {len(columns)} columns")


def export(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = [name for name in reader.fieldnames if name not in IGNORED_EXPORT]
        records = list(reader)
    key = lambda row: (row["precursor_id"], row["peak_group_rank"])
    keys = [key(row) for row in records]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate export precursor/rank key")
    return columns, [tuple(row[name] for name in columns) for row in sorted(records, key=key)]


def database(work_dir):
    files = list(work_dir.rglob("workflow.osw"))
    if len(files) != 1:
        raise ValueError(f"Expected one workflow database in {work_dir}, found {len(files)}")
    with sqlite3.connect(f"file:{files[0]}?mode=ro", uri=True) as connection:
        features = connection.execute("SELECT ID, PRECURSOR_ID FROM FEATURE ORDER BY PRECURSOR_ID, EXP_RT, ID").fetchall()
        ranks = {}
        feature_keys = {}
        for feature_id, precursor_id in features:
            ranks[precursor_id] = ranks.get(precursor_id, 0) + 1
            feature_keys[feature_id] = (precursor_id, ranks[precursor_id])
        names = sorted(row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND (name LIKE 'FEATURE%' OR name LIKE 'SCORE%')"))
        result = {}
        for name in names:
            # Table names come from the generated SQLite schema, not command text.
            quoted = '"' + name.replace('"', '""') + '"'
            cursor = connection.execute(f"SELECT * FROM {quoted}")
            source_columns = [column[0] for column in cursor.description]
            columns = [column for column in source_columns if column != "RUN_ID"]
            rows = []
            for values in cursor:
                record = dict(zip(source_columns, values))
                if name == "FEATURE":
                    record["ID"] = feature_keys[record["ID"]]
                if "FEATURE_ID" in record:
                    record["FEATURE_ID"] = feature_keys[record["FEATURE_ID"]]
                rows.append(tuple(record[column] for column in columns))
            result[name] = (columns, sorted(rows, key=repr))
    if not result["FEATURE"][1] or not result["FEATURE_MS2"][1]:
        raise ValueError("Missing feature or MS2 evidence rows")
    return result


def compare(left_out, right_out, left_work, right_work, basename):
    a_columns, a_rows = export(left_out / f"{basename}.results.tsv")
    b_columns, b_rows = export(right_out / f"{basename}.results.tsv")
    if a_columns != b_columns:
        raise ValueError("Export columns differ")
    compare_rows("export", a_columns, a_rows, b_rows)
    left, right = database(left_work), database(right_work)
    if left.keys() != right.keys():
        raise ValueError("Feature/score table sets differ")
    for name in left:
        a_columns, a_rows = left[name]
        b_columns, b_rows = right[name]
        if a_columns != b_columns:
            raise ValueError(f"{name}: columns differ")
        if a_rows or b_rows:
            compare_rows(name, a_columns, a_rows, b_rows)
    print("PASS: identities, evidence, scores and quantities match. Only run/file/feature identities are normalized.")


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        checks = [("1;2", "1.0000001;2", True), ("1;2", "1;3", False),
                  (None, "0", False), ("0", "", False), ("peptideA", "peptideB", False),
                  (((1, 2), 5.0), ((1, 3), 5.0), False), (1, float("inf"), False)]
        if not all(same_value(a, b) == expected for a, b, expected in checks):
            raise ValueError("Comparison self-test failed")
        print("PASS: numeric, identity, missing/zero and aggregate mismatch checks")
    elif len(sys.argv) == 6:
        compare(*(Path(value) for value in sys.argv[1:5]), sys.argv[5])
    else:
        raise SystemExit("usage: compare_results.py LEFT_OUT RIGHT_OUT LEFT_WORK RIGHT_WORK BASENAME")
