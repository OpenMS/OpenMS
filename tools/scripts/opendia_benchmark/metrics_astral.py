#!/usr/bin/env python3
"""Report Astral IDs and native-matrix replicate metrics, not biological proteins."""

import argparse
import json
import math
from pathlib import Path

from benchmark_metrics import deduplicate, mapping_kind, number, views
from validate_timsim import fdr_curve, read_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--q-threshold', type=float, default=.01)
    parser.add_argument('--intensity-column', default='Intensity')
    parser.add_argument('--runs', nargs='+', help='Expected runs, including runs with no exported rows')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.q_threshold) or not 0 <= args.q_threshold <= 1:
        parser.error('--q-threshold must be between 0 and 1')
    if args.q_threshold != .01:
        parser.error('--q-threshold must be 0.01 for the fixed 1% benchmark views')
    raw = read_tsv(args.results)
    if raw and 'm_score' not in raw[0]:
        raise SystemExit('Missing m_score column, no q-value fallback is permitted')
    runs = set(args.runs or []) | {row['run_name'] for row in raw}
    # Native matrix headers preserve runs even when there are no IDs in them.
    matrix_path = args.results.with_name(args.results.name.replace('.results.tsv', '.precursor.matrix.tsv'))
    if matrix_path != args.results and matrix_path.exists():
        import csv
        with matrix_path.open() as handle:
            columns = next(csv.reader(handle, delimiter='\t'), [])
        runs.update(c for c in columns if c not in {'transition_group_id', 'Sequence', 'FullPeptideName', 'Charge', 'ProteinName', 'GeneName'})
    runs = sorted(runs)
    conditions = {run: 'astral' for run in runs}
    rows = list(deduplicate(raw).values())
    per_run = {}
    for run in runs:
        subset = [r for r in rows if r['run_name'] == run]
        targets = [r for r in subset if mapping_kind(r) == 'target']
        entrapments = [r for r in subset if mapping_kind(r) == 'entrapment']
        per_run[run] = {'rows': len(subset), 'target_rows': len(targets), 'entrapment_rows': len(entrapments),
            'ambiguous_mixed_mapping_rows': sum(mapping_kind(r) == 'ambiguous' for r in subset),
            'raw_entrapment_share': len(entrapments) / (len(targets) + len(entrapments)) if targets or entrapments else None,
            'unique_target_precursors': len(targets),
            'unique_target_peptides': len({r['FullPeptideName'] for r in targets}),
            'unique_target_fasta_groups': len({r['ProteinName'] for r in targets})}
    targets = [r for r in rows if mapping_kind(r) == 'target']
    entrapments = [r for r in rows if mapping_kind(r) == 'entrapment']
    classified = [(number(r['m_score']), mapping_kind(r) == 'entrapment') for r in rows if mapping_kind(r) != 'ambiguous']
    summary = {'schema_version': 2, 'results': str(args.results.resolve()), 'q_threshold': .01,
        'q_column': 'm_score', 'intensity_column': args.intensity_column, 'runs': runs, 'conditions': conditions,
        'protein_label_terminology': 'FASTA groups with peptide pseudo-protein headers, not biological protein identifications or quantification',
        'entrapment_interpretation': 'Raw entrapment share, not calibrated FDR. _p_target marks entrapment. Mixed mappings excluded.',
        'views': views(raw, runs, conditions, args.results, read_tsv, astral=True),
        'fdr_curve': fdr_curve(classified, .01),
        'fdr_curve_per_run': {run: fdr_curve([(number(r['m_score']), mapping_kind(r) == 'entrapment') for r in rows
                                           if r['run_name'] == run and mapping_kind(r) != 'ambiguous'], .01) for run in runs},
        'identification': {'target_rows': len(targets), 'entrapment_rows': len(entrapments),
            'ambiguous_mixed_mapping_rows': sum(mapping_kind(r) == 'ambiguous' for r in rows),
            'raw_entrapment_share': len(entrapments) / (len(targets) + len(entrapments)) if targets or entrapments else None,
            'unique_target_precursors': len({(r['FullPeptideName'], r['Charge']) for r in targets}),
            'unique_entrapment_precursors': len({(r['FullPeptideName'], r['Charge']) for r in entrapments}), 'per_run': per_run},
        'quantification_repeatability': {'note': 'Use views.standard.quantification native matrices. No quantitative ground truth is supplied.'}}
    text = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
