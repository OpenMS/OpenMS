"""Shared, standard-library-only benchmark filtering and replicate statistics."""

import itertools
import math
import re
import statistics
from collections import defaultdict


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def valid_q(value, threshold):
    value = number(value)
    return math.isfinite(value) and 0 <= value <= threshold


def key(row):
    return (row['run_name'], row.get('FullPeptideName') or row['Sequence'], int(row['Charge']))


def deduplicate(rows, threshold=0.01):
    best = {}
    for row in rows:
        if str(row.get('decoy', '0')).lower() not in {'0', 'false', ''}:
            continue
        if not valid_q(row.get('m_score'), threshold):
            continue
        intensity = number(row.get('Intensity'))
        rank = (number(row['m_score']), -intensity if math.isfinite(intensity) else math.inf)
        identity = key(row)
        if identity not in best or rank < best[identity][0]:
            best[identity] = (rank, row)
    return {identity: item[1] for identity, item in best.items()}


def mapping_kind(row):
    labels = [label for label in re.split(r'[;,]', row.get('ProteinName', '')) if label]
    entrapment = ['_p_target' in label for label in labels]
    if not labels or (any(entrapment) and not all(entrapment)):
        return 'ambiguous'
    return 'entrapment' if all(entrapment) else 'target'


def stats(values):
    return {'values': values, 'n': len(values),
            'mean': statistics.fmean(values) if values else None,
            'sample_sd': statistics.stdev(values) if len(values) > 1 else None}


def pearson(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    denominator = math.sqrt(sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))
    return sum((x-mx)*(y-my) for x, y in zip(xs, ys)) / denominator if denominator else None


def entity_metrics(values, runs):
    completeness = defaultdict(int)
    cvs = []
    for run_values in values.values():
        present = [run_values[run] for run in runs if run in run_values and math.isfinite(run_values[run]) and run_values[run] > 0]
        completeness[len(present)] += 1
        if len(present) == len(runs) and len(present) > 1:
            cvs.append(statistics.stdev(present) / statistics.fmean(present))
    correlations = {}
    for left, right in itertools.combinations(runs, 2):
        pairs = [(math.log2(v[left]), math.log2(v[right])) for v in values.values()
                 if math.isfinite(v.get(left, 0)) and math.isfinite(v.get(right, 0)) and v.get(left, 0) > 0 and v.get(right, 0) > 0]
        correlations[f'{left}__{right}'] = {'shared': len(pairs), 'n': len(pairs),
            'pearson_log2': pearson([p[0] for p in pairs], [p[1] for p in pairs])}
    return {'entities': len(values), 'completeness': {str(k): v for k, v in sorted(completeness.items())},
            'complete_entities': completeness[len(runs)], 'cv_n': len(cvs), 'cv_values': cvs,
            'median_cv_complete': statistics.median(cvs) if cvs else None,
            'cv_definition': 'sample SD (ddof=1) / mean, positive complete entities only',
            'pairwise_log2_correlations': correlations}


def views(rows, runs, conditions, results_path, read_tsv, astral=False):
    selected = deduplicate(rows)
    # Unknown global values invalidate the standard view, never silently pass it.
    globals_ = ('m_score_peptide_global', 'm_score_protein_global')
    schema_available = bool(rows) and all(column in rows[0] for column in globals_)
    available = schema_available and all(
        all(math.isfinite(number(row.get(column))) and 0 <= number(row[column]) <= 1 for column in globals_)
        for row in selected.values())
    standard = {k: row for k, row in selected.items() if all(valid_q(row.get(c), .01) for c in globals_)}
    result = {}
    for name, chosen in [('precursor_only', selected), ('standard', standard)]:
        if name == 'standard' and not available:
            result[name] = {'status': 'NA', 'reason': 'Missing or nonfinite global peptide/protein q-values', 'per_run': None, 'count_summary': None, 'quantification': None}
            continue
        per_run = {}
        for run in runs:
            subset = [r for r in chosen.values() if r['run_name'] == run and (not astral or mapping_kind(r) == 'target')]
            per_run[run] = {'precursor': len(subset), 'peptide': len({r.get('FullPeptideName') or r['Sequence'] for r in subset}),
                            'protein_label': len({r['ProteinName'] for r in subset})}
        count_summary = {condition: {level: stats([per_run[r][level] for r in runs if conditions[r] == condition])
                         for level in ('precursor', 'peptide', 'protein_label')} for condition in sorted(set(conditions.values()))}
        result[name] = {'status': 'ok', 'per_run': per_run, 'count_summary': count_summary}
    if result['standard']['status'] == 'ok':
        quant = {}
        for level, suffix in [('precursor', 'precursor'), ('peptide', 'peptide'), ('protein_label', 'protein')]:
            path = results_path.with_name(results_path.name.replace('.results.tsv', f'.{suffix}.matrix.tsv'))
            if path == results_path or not path.exists():
                quant[level] = {'source': None, 'status': 'NA', 'reason': 'Native matrix unavailable'}
                continue
            matrix = read_tsv(path)
            values = {}
            allowed = {r['ProteinName'] for r in standard.values() if not astral or mapping_kind(r) == 'target'}
            allowed_peptides = {r.get('FullPeptideName') or r['Sequence'] for r in standard.values() if not astral or mapping_kind(r) == 'target'}
            allowed_precursors = {(r.get('FullPeptideName') or r['Sequence'], r['Charge']) for r in standard.values() if not astral or mapping_kind(r) == 'target'}
            for index, row in enumerate(matrix):
                if ((level == 'protein_label' and row.get('ProteinName') not in allowed) or
                    (level == 'peptide' and row.get('FullPeptideName') not in allowed_peptides) or
                    (level == 'precursor' and (row.get('FullPeptideName'), row.get('Charge')) not in allowed_precursors)):
                    continue
                identity = row.get('ProteinName') if level == 'protein_label' else row.get('FullPeptideName')
                if level == 'precursor':
                    identity = f"{identity}/{row['Charge']}"
                values[identity or str(index)] = {run: number(row.get(run)) for run in runs}
            quant[level] = {'source': 'native_matrix', 'path': str(path), 'status': 'ok',
                'conditions': {c: entity_metrics(values, [r for r in runs if conditions[r] == c]) for c in sorted(set(conditions.values()))}}
        result['standard']['quantification'] = quant
    return result
