#!/usr/bin/env python3
"""Validate six-run OpenDIA output against packaged TimSim truth."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from benchmark_metrics import deduplicate, views


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
    if rows and 'm_score' not in rows[0]:
        raise SystemExit(f"Missing m_score column in {path}")
    return deduplicate(rows, min(q_threshold, 0.01)), 'm_score'


def fdr_curve(classified, maximum):
    """Report raw entrapment shares, not search-space-adjusted FDR estimates."""
    points = []
    for threshold in sorted({q for q in (0.0001, 0.0005, 0.001, 0.005, maximum) if q <= maximum}):
        selected = [entrapment for q, entrapment in classified if q <= threshold]
        points.append({
            "q_threshold": threshold,
            "rows": len(selected),
            "entrapment_rows": sum(selected),
            "entrapment_share": sum(selected) / len(selected) if selected else None,
        })
    return points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("target_results", type=Path)
    parser.add_argument("entrapment_results", type=Path)
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.q_threshold) or not 0 <= args.q_threshold <= 1:
        parser.error("--q-threshold must be between 0 and 1")
    if args.q_threshold != 0.01:
        parser.error("--q-threshold must be 0.01 for the fixed 1% benchmark views")

    truth_dir = args.dataset / "truth"
    manifest = json.loads((args.dataset / "metadata/fixture_manifest.json").read_text())
    run_names = set(manifest["runs"])
    realized = read_tsv(truth_dir / "OpenSwathTimSim.realized_truth.tsv")
    entrapment_truth = read_tsv(truth_dir / "OpenSwathTimSim.entrapment_truth.tsv")
    target_results, target_q = best_results(args.target_results, args.q_threshold, ('m_score',))
    entrapment_results, entrapment_q = best_results(args.entrapment_results, args.q_threshold, ('m_score',))

    conditions = {run: ('control' if 'control' in run.lower() else 'treatment') for run in run_names}
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
        per_run[run] = {"expected": expected, "identified": identified, "recovery": identified / expected if expected else None}

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

    native_effects = {}
    native_views = views(read_tsv(args.target_results), sorted(run_names), conditions, args.target_results, read_tsv)
    if native_views['standard']['status'] == 'ok':
        for level, suffix in [('precursor', 'precursor'), ('peptide', 'peptide'), ('protein_label', 'protein')]:
            matrix_path = args.target_results.with_name(args.target_results.name.replace('.results.tsv', f'.{suffix}.matrix.tsv'))
            if not matrix_path.exists() or matrix_path == args.target_results:
                continue
            expected = defaultdict(lambda: defaultdict(float))
            for row in realized:
                entity = row['ProteinId'] if level == 'protein_label' else row['PeptideSequence'].strip().upper()
                if level == 'precursor':
                    entity = f"{entity}/{row['PrecursorCharge']}"
                expected[entity][row['RunName']] += number(row['RealizedInputEvents'])
            effects = []
            for row in read_tsv(matrix_path):
                entity = row['ProteinName'] if level == 'protein_label' else row['FullPeptideName']
                if level == 'precursor':
                    entity = f"{entity}/{row['Charge']}"
                if entity not in expected:
                    continue
                groups = {c: [r for r in sorted(run_names) if conditions[r] == c and math.isfinite(number(row.get(r))) and number(row[r]) > 0] for c in ('control', 'treatment')}
                if not all(groups.values()):
                    continue
                observed_fc = log2_ratio(mean([number(row[r]) for r in groups['treatment']]), mean([number(row[r]) for r in groups['control']]))
                expected_fc = log2_ratio(mean([expected[entity].get(r, 0) for r in groups['treatment']]), mean([expected[entity].get(r, 0) for r in groups['control']]))
                if observed_fc is not None and expected_fc is not None:
                    effects.append({'entity': entity, 'realized_log2_fc': expected_fc, 'observed_log2_fc': observed_fc})
            native_effects[level] = {'source': 'native_matrix', 'n': len(effects), 'effects': effects,
                'pearson': pearson([e['realized_log2_fc'] for e in effects], [e['observed_log2_fc'] for e in effects]),
                'truth_note': 'Realized input events summed by entity, matched observed run support, not nominal fold changes'}

    summary = {
        'schema_version': 2,
        'runs': sorted(run_names),
        'conditions': conditions,
        'views': native_views,
        'entrapment_views': views([row for row in read_tsv(args.entrapment_results)
            if (row['Sequence'].strip().upper(), int(row['Charge'])) in target_keys], sorted(run_names), conditions, args.entrapment_results, read_tsv),
        'fdr_curve_per_run': {run: fdr_curve([(number(row['m_score']), (sequence, charge) in entrapment_keys)
             for (r, sequence, charge), row in entrapment_results.items() if r == run
             and ((sequence, charge) in target_keys or (sequence, charge) in entrapment_keys)], .01) for run in sorted(run_names)},
        "fdr_curve": fdr_curve([
            (number(row[entrapment_q]), (sequence, charge) in entrapment_keys)
            for (_, sequence, charge), row in entrapment_results.items()
            if (sequence, charge) in target_keys or (sequence, charge) in entrapment_keys
        ], args.q_threshold),
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
            'native_standard_effects': native_effects,
            'note': 'Legacy protein intensity is a descriptive sum, not native C++ top-N mean quantification. Use views.standard.quantification for native repeatability.',
            'precursor_effect_values': [{'realized_log2_fc': e, 'observed_log2_fc': o} for e, o in zip(expected_precursor_fc, observed_precursor_fc)],
            'protein_descriptive_sum_effect_values': [{'realized_log2_fc': e, 'observed_log2_fc': o} for e, o in zip(expected_protein_fc, observed_protein_fc)],
            "precursor_effect_pairs": len(expected_precursor_fc),
            "precursor_realized_vs_observed_pearson": pearson(expected_precursor_fc, observed_precursor_fc),
            "protein_effect_pairs": len(expected_protein_fc),
            "protein_realized_vs_observed_pearson": pearson(expected_protein_fc, observed_protein_fc),
        },
        "entrapment_export": {
            "target_rows": target_rows,
            "entrapment_rows": entrapment_rows,
            "raw_entrapment_share": entrapment_rows / max(target_rows + entrapment_rows, 1),
            'interpretation': 'Raw entrapment share, not calibrated FDR',
            "unique_target_precursors": unique_target,
            "unique_entrapment_precursors": unique_entrapment,
            "unique_precursor_raw_entrapment_share": unique_entrapment / max(unique_target + unique_entrapment, 1),
        },
    }
    text = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
