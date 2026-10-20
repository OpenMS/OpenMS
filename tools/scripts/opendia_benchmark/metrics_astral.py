#!/usr/bin/env python3
"""Report identification and replicate-quantification metrics for Astral results."""

import argparse
import csv
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def pearson(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return None
    xs, ys = zip(*pairs)
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    denominator = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return sum((x - mx) * (y - my) for x, y in pairs) / denominator if denominator else None


def entity_metrics(values, runs):
    completeness = defaultdict(int)
    cvs = []
    for run_values in values.values():
        present = [value for value in run_values.values() if value > 0 and math.isfinite(value)]
        completeness[len(present)] += 1
        if len(present) == len(runs):
            mean = statistics.fmean(present)
            cvs.append(statistics.pstdev(present) / mean if mean else math.nan)
    correlations = {}
    for left, right in itertools.combinations(runs, 2):
        pairs = [
            (math.log2(v[left]), math.log2(v[right]))
            for v in values.values()
            if v.get(left, 0) > 0 and v.get(right, 0) > 0
        ]
        correlations[f"{left}__{right}"] = {
            "shared": len(pairs),
            "pearson_log2": pearson([p[0] for p in pairs], [p[1] for p in pairs]),
        }
    return {
        "entities": len(values),
        "completeness": {str(key): completeness[key] for key in sorted(completeness)},
        "complete_entities": completeness[len(runs)],
        "median_cv_complete": statistics.median(cvs) if cvs else None,
        "pairwise_log2_correlations": correlations,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--intensity-column", default="Intensity")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    with args.results.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"run_name", "Sequence", "FullPeptideName", "ProteinName", "Charge", args.intensity_column}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing result columns: {sorted(missing)}")
        q_column = next((name for name in ("m_score", "ms2_m_score", "m_score_peptide_run_specific") if name in reader.fieldnames), None)
        if q_column is None:
            raise SystemExit(f"Missing q-value column in {args.results}")
        rows = []
        for row in reader:
            if str(row.get("decoy", "0")).lower() not in {"0", "false", ""}:
                continue
            q_value = number(row.get(q_column))
            if math.isfinite(q_value) and q_value > args.q_threshold:
                continue
            rows.append(row)

    runs = sorted({row["run_name"] for row in rows})
    per_run = {}
    precursor_values = defaultdict(dict)
    protein_values = defaultdict(lambda: defaultdict(float))
    target_rows = entrapment_rows = 0
    target_precursors = set()
    entrapment_precursors = set()

    for row in rows:
        run = row["run_name"]
        precursor = (row["FullPeptideName"], int(row["Charge"]))
        proteins = row["ProteinName"]
        is_entrapment = "_p_target" in proteins
        intensity = number(row[args.intensity_column])
        if is_entrapment:
            entrapment_rows += 1
            entrapment_precursors.add(precursor)
            continue
        target_rows += 1
        target_precursors.add(precursor)
        if math.isfinite(intensity) and intensity > 0:
            precursor_values[precursor][run] = precursor_values[precursor].get(run, 0.0) + intensity
            protein_values[proteins][run] += intensity

    for run in runs:
        selected = [row for row in rows if row["run_name"] == run]
        entrapment = sum("_p_target" in row["ProteinName"] for row in selected)
        target = len(selected) - entrapment
        per_run[run] = {
            "rows": len(selected),
            "target_rows": target,
            "entrapment_rows": entrapment,
            "empirical_fdp": entrapment / max(len(selected), 1),
            "unique_target_precursors": len({(row["FullPeptideName"], row["Charge"]) for row in selected if "_p_target" not in row["ProteinName"]}),
            "unique_target_peptides": len({row["Sequence"] for row in selected if "_p_target" not in row["ProteinName"]}),
            "unique_target_proteins": len({row["ProteinName"] for row in selected if "_p_target" not in row["ProteinName"]}),
        }

    summary = {
        "results": str(args.results.resolve()),
        "q_threshold": args.q_threshold,
        "q_column": q_column,
        "intensity_column": args.intensity_column,
        "runs": runs,
        "identification": {
            "target_rows": target_rows,
            "entrapment_rows": entrapment_rows,
            "empirical_fdp": entrapment_rows / max(target_rows + entrapment_rows, 1),
            "unique_target_precursors": len(target_precursors),
            "unique_entrapment_precursors": len(entrapment_precursors),
            "unique_precursor_empirical_fdp": len(entrapment_precursors) / max(len(target_precursors) + len(entrapment_precursors), 1),
            "per_run": per_run,
        },
        "quantification_repeatability": {
            "note": "Replicate consistency only; this dataset directory contains no quantitative ground truth.",
            "precursor": entity_metrics(precursor_values, runs),
            "protein": entity_metrics(protein_values, runs),
        },
    }
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
