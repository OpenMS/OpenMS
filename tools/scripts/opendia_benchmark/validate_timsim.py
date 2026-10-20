#!/usr/bin/env python3
"""Validate six-run OpenDIA output against packaged TimSim truth."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def pearson(xs, ys):
    pairs = [
        (x, y)
        for x, y in zip(xs, ys)
        if x is not None and y is not None and math.isfinite(x) and math.isfinite(y)
    ]
    if len(pairs) < 2:
        return None
    xs, ys = zip(*pairs)
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    numerator = sum((x - mx) * (y - my) for x, y in pairs)
    denominator = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return numerator / denominator if denominator else None


def mean(values):
    values = [value for value in values if math.isfinite(value)]
    return sum(values) / len(values) if values else None


def log2_ratio(treatment, control):
    return math.log2(treatment / control) if treatment > 0 and control > 0 else None


def best_results(path, q_threshold, q_candidates):
    rows = read_tsv(path)
    if not rows:
        raise SystemExit(f"No rows in {path}")
    required = {"run_name", "Sequence", "Charge", "Intensity", "RT", "EXP_IM"}
    missing = required - rows[0].keys()
    if missing:
        raise SystemExit(f"Missing result columns: {sorted(missing)}")
    q_column = next((name for name in q_candidates if name in rows[0]), None)
    if q_column is None:
        raise SystemExit(f"Missing q-value column in {path}")
    best = {}
    for row in rows:
        if str(row.get("decoy", "0")).lower() not in {"0", "false", ""}:
            continue
        q_value = number(row.get(q_column))
        if math.isfinite(q_value) and q_value > q_threshold:
            continue
        key = (row["run_name"], row["Sequence"].strip().upper(), int(row["Charge"]))
        rank = (q_value if math.isfinite(q_value) else math.inf, -number(row.get("d_score")))
        if key not in best or rank < best[key][0]:
            best[key] = (rank, row)
    return {key: row for key, (_, row) in best.items()}, q_column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("target_results", type=Path)
    parser.add_argument("entrapment_results", type=Path)
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    truth_dir = args.dataset / "truth"
    manifest = json.loads((args.dataset / "metadata/fixture_manifest.json").read_text())
    run_names = set(manifest["runs"])
    realized = read_tsv(truth_dir / "OpenSwathTimSim.realized_truth.tsv")
    entrapment_truth = read_tsv(truth_dir / "OpenSwathTimSim.entrapment_truth.tsv")
    target_results, target_q = best_results(
        args.target_results,
        args.q_threshold,
        ("m_score", "ms2_m_score", "m_score_peptide_run_specific", "m_score_peptide_global"),
    )
    entrapment_results, entrapment_q = best_results(
        args.entrapment_results,
        args.q_threshold,
        ("m_score_peptide_global", "m_score", "ms2_m_score", "m_score_peptide_run_specific"),
    )

    conditions = {}
    expected_keys = set()
    target_keys = set()
    precursor_expected = defaultdict(dict)
    precursor_observed = defaultdict(dict)
    protein_expected = defaultdict(lambda: defaultdict(float))
    protein_observed = defaultdict(lambda: defaultdict(float))
    rt_errors, im_errors = [], []
    per_run = {}

    for row in realized:
        run = row["RunName"]
        if run not in run_names:
            continue
        sequence = row["PeptideSequence"].strip().upper()
        charge = int(row["PrecursorCharge"])
        key = (run, sequence, charge)
        expected_keys.add(key)
        target_keys.add((sequence, charge))
        conditions[run] = row["Condition"]
        realized_intensity = number(row["RealizedInputEvents"])
        precursor = (sequence, charge)
        precursor_expected[precursor][run] = realized_intensity
        protein_expected[row["ProteinId"]][run] += realized_intensity
        observed = target_results.get(key)
        if observed:
            intensity = number(observed["Intensity"])
            precursor_observed[precursor][run] = intensity
            protein_observed[row["ProteinId"]][run] += intensity
            rt_errors.append(abs(number(observed["RT"]) - number(row["RealizedRTApex"])))
            im_errors.append(abs(number(observed["EXP_IM"]) - number(row["RealizedIMApexSQLite"])))

    for run in sorted(run_names):
        expected = sum(key[0] == run for key in expected_keys)
        identified = sum(key in target_results for key in expected_keys if key[0] == run)
        per_run[run] = {"expected": expected, "identified": identified, "recovery": identified / expected}

    expected_precursor_fc = []
    observed_precursor_fc = []
    for precursor in sorted(precursor_expected):
        expected_runs = precursor_expected[precursor]
        observed_runs = precursor_observed[precursor]
        control_runs = [run for run in observed_runs if conditions[run] == "control"]
        treatment_runs = [run for run in observed_runs if conditions[run] == "treatment"]
        control_e = [expected_runs[run] for run in control_runs]
        treatment_e = [expected_runs[run] for run in treatment_runs]
        control_o = [observed_runs[run] for run in control_runs]
        treatment_o = [observed_runs[run] for run in treatment_runs]
        if control_o and treatment_o:
            expected_precursor_fc.append(log2_ratio(sum(treatment_e) / len(treatment_e), sum(control_e) / len(control_e)))
            observed_precursor_fc.append(log2_ratio(sum(treatment_o) / len(treatment_o), sum(control_o) / len(control_o)))

    expected_protein_fc, observed_protein_fc = [], []
    for protein in sorted(protein_expected):
        expected_runs = protein_expected[protein]
        observed_runs = protein_observed[protein]
        control_runs = [run for run in observed_runs if conditions[run] == "control"]
        treatment_runs = [run for run in observed_runs if conditions[run] == "treatment"]
        control_e = [expected_runs[run] for run in control_runs]
        treatment_e = [expected_runs[run] for run in treatment_runs]
        control_o = [observed_runs[run] for run in control_runs]
        treatment_o = [observed_runs[run] for run in treatment_runs]
        if control_o and treatment_o:
            expected_protein_fc.append(log2_ratio(sum(treatment_e) / len(treatment_e), sum(control_e) / len(control_e)))
            observed_protein_fc.append(log2_ratio(sum(treatment_o) / len(treatment_o), sum(control_o) / len(control_o)))

    entrapment_keys = {(row["EntrapmentSequence"].strip().upper(), int(row["PrecursorCharge"])) for row in entrapment_truth}
    classified = []
    for (run, sequence, charge), row in entrapment_results.items():
        key = (sequence, charge)
        if key in target_keys:
            classified.append((run, "target"))
        elif key in entrapment_keys:
            classified.append((run, "entrapment"))
    target_rows = sum(kind == "target" for _, kind in classified)
    entrapment_rows = sum(kind == "entrapment" for _, kind in classified)
    unique_target = len({(sequence, charge) for (_, sequence, charge) in entrapment_results if (sequence, charge) in target_keys})
    unique_entrapment = len({(sequence, charge) for (_, sequence, charge) in entrapment_results if (sequence, charge) in entrapment_keys})

    summary = {
        "q_threshold": args.q_threshold,
        "target_q_column": target_q,
        "entrapment_q_column": entrapment_q,
        "target_identification": {
            "expected_target_run_pairs": len(expected_keys),
            "identified_target_run_pairs": sum(key in target_results for key in expected_keys),
            "overall_recovery": sum(key in target_results for key in expected_keys) / len(expected_keys),
            "per_run": per_run,
            "rt_mae_seconds": mean(rt_errors),
            "im_mae_1_over_k0": mean(im_errors),
        },
        "quantification": {
            "precursor_effect_pairs": len(expected_precursor_fc),
            "precursor_realized_vs_observed_pearson": pearson(expected_precursor_fc, observed_precursor_fc),
            "protein_effect_pairs": len(expected_protein_fc),
            "protein_realized_vs_observed_pearson": pearson(expected_protein_fc, observed_protein_fc),
        },
        "entrapment_export": {
            "target_rows": target_rows,
            "entrapment_rows": entrapment_rows,
            "pooled_empirical_fdp": entrapment_rows / max(target_rows + entrapment_rows, 1),
            "unique_target_precursors": unique_target,
            "unique_entrapment_precursors": unique_entrapment,
            "unique_precursor_empirical_fdp": unique_entrapment / max(unique_target + unique_entrapment, 1),
        },
    }
    text = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
