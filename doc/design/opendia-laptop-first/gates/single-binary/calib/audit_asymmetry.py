#!/usr/bin/env python3
"""Descriptive, label-blind candidate/sibling audit. See ASYMMETRY_PLAN.md."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3


def read_tsv(path, required):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f'{path}: missing columns {sorted(required - set(reader.fieldnames or []))}')
        return list(reader)


def number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('nonfinite numeric value')
    return result


def library(path):
    uri = Path(path).resolve().as_uri() + '?mode=ro'
    with sqlite3.connect(uri, uri=True) as connection:
        rows = connection.execute('''SELECT p.ID, p.TRAML_ID, p.DECOY, p.CHARGE,
            pe.MODIFIED_SEQUENCE FROM PRECURSOR p
            JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
            JOIN PEPTIDE pe ON pe.ID=m.PEPTIDE_ID''').fetchall()
        ids = {str(row[0]) for row in connection.execute('SELECT ID FROM PRECURSOR')}
    precursors = {}
    names = set()
    for pid, name, side, charge, sequence in rows:
        pid = str(pid)
        if pid in precursors or name in names or side not in (0, 1) or not sequence:
            raise ValueError('library requires unique IDs/TRAML_IDs and one peptide per precursor')
        if side and not name.startswith('DECOY_'):
            raise ValueError('decoy lacks own-target DECOY_ prefix')
        names.add(name)
        precursors[pid] = dict(name=name, side=side, charge=charge, sequence=sequence,
                               pair=name[6:] if side else name)
    if set(precursors) != ids:
        raise ValueError('unmapped library precursor')
    return precursors


def audit(library_path, features_path=None, reference_path=None, tolerance=3.0):
    precursors = library(library_path)
    groups = collections.defaultdict(set)
    pairs = collections.defaultdict(dict)
    for pid, precursor in precursors.items():
        groups[precursor['side'], precursor['sequence']].add(pid)
        pair = pairs[precursor['pair']]
        if precursor['side'] in pair:
            raise ValueError('ambiguous own-decoy pair')
        pair[precursor['side']] = pid
    windows = (read_tsv(features_path, {'candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value', 'fold'})
               if features_path else [])
    seen = set()
    by_precursor = collections.defaultdict(list)
    by_group = collections.defaultdict(list)
    for row in windows:
        pid = row['precursor_id']
        if pid not in precursors or int(row['is_decoy']) != precursors[pid]['side']:
            raise ValueError('window/library ID or side mismatch')
        if row['candidate_row'] in seen:
            raise ValueError('duplicate candidate_row')
        seen.add(row['candidate_row'])
        for field in ('score', 'q_value', 'rt_apex', 'fold'):
            number(row[field])
        row['rt'] = number(row['rt_apex'])
        by_precursor[pid].append(row)
        p = precursors[pid]
        by_group[p['side'], p['sequence']].append(row)
    reference = {}
    if reference_path:
        for row in read_tsv(reference_path, {'precursor_id', 'is_decoy', 'score'}):
            pid = row['precursor_id']
            if pid not in precursors or int(row['is_decoy']) != precursors[pid]['side']:
                raise ValueError('reference/library mismatch')
            score = number(row['score'])
            if score > -1e29:
                reference[pid] = max(reference.get(pid, -math.inf), score)
    records = []
    for pid, p in sorted(precursors.items()):
        siblings = groups[p['side'], p['sequence']] - {pid}
        own = by_precursor[pid]
        co_counts = []
        for row in own:
            co = {other['precursor_id'] for other in by_group[p['side'], p['sequence']]
                  if other['precursor_id'] != pid and abs(other['rt'] - row['rt']) <= tolerance}
            co_counts.append(len(co))
            # Exported s_n/s_co are the actual first-pass features if present.
            if 's_n' in row and number(row['s_n']) != len(siblings):
                raise ValueError('exported s_n disagrees with library grouping')
            if 's_co' in row and number(row['s_co']) != len(co):
                raise ValueError('exported s_co disagrees with window availability')
        records.append(dict(precursor_id=pid, is_decoy=p['side'], pair_key=p['pair'],
                            library_siblings=len(siblings), candidate_windows=len(own) if features_path else None,
                            candidate_siblings=sum(bool(by_precursor[s]) for s in siblings) if features_path else None,
                            windows_with_coeluting_sibling=sum(c > 0 for c in co_counts) if features_path else None,
                            max_coeluting_siblings=max(co_counts, default=None)))
    by_id = {r['precursor_id']: r for r in records}
    paired = []
    for key, sides in sorted(pairs.items()):
        if set(sides) != {0, 1}:
            continue
        target, decoy = sides[0], sides[1]
        region = 'not_requested'
        if reference_path:
            region = ('missing_reference' if target not in reference or decoy not in reference else
                      'low_reference_both' if max(reference[target], reference[decoy]) <= math.log(2) else
                      'outside_low_reference')
        record = dict(pair_key=key, target_id=target, decoy_id=decoy, reference_region=region)
        for metric in ('library_siblings', 'candidate_windows', 'candidate_siblings', 'windows_with_coeluting_sibling'):
            record['target_' + metric] = by_id[target][metric]
            record['decoy_' + metric] = by_id[decoy][metric]
            t_value, d_value = by_id[target][metric], by_id[decoy][metric]
            record['delta_' + metric] = t_value - d_value if t_value is not None and d_value is not None else None
        paired.append(record)
    return dict(scope='descriptive_only_not_null_exchangeability_test',
                candidate_population='post_pair_support_population_before_final_keep' if features_path else 'unavailable_library_only',
                rt_tolerance_seconds=tolerance, precursors=records, own_decoy_pairs=paired,
                unpaired_precursor_ids=sorted(pid for sides in pairs.values() if len(sides) != 2 for pid in sides.values()),
                gaps=['No pre-population candidate counts or extraction search-space coverage.',
                      'No first-pass sibling scores unless separately retained.',
                      'Low-reference region is a diagnostic proxy, not known null truth.'])


def digest(path):
    result = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', required=True)
    population = parser.add_mutually_exclusive_group(required=True)
    population.add_argument('--features', help='One run .features.tsv, never .precursors.tsv')
    population.add_argument('--library-only', action='store_true', help='Structural sibling counts only, candidate fields null')
    parser.add_argument('--reference', help='Frozen no-GBM own-decoy score table for the same run/seed')
    parser.add_argument('--output', required=True)
    parser.add_argument('--synthetic', action='store_true', help='Use only with fabricated fixtures')
    parser.add_argument('--precheck-commit')
    parser.add_argument('--precheck-pushed', action='store_true')
    args = parser.parse_args()
    if not args.synthetic and not (re.fullmatch('[0-9a-f]{40}', args.precheck_commit or '') and args.precheck_pushed):
        parser.error('real data blocked until precheck committed/pushed, supply full commit and push attestation')
    if args.library_only and args.reference:
        parser.error('--reference requires --features, library-only mode has no score strata')
    result = audit(args.library, args.features, args.reference)
    result['provenance'] = dict(synthetic=args.synthetic, precheck_commit=args.precheck_commit,
                                pushed_attestation=args.precheck_pushed,
                                inputs={str(Path(p).resolve()): digest(p) for p in
                                        (args.library, args.features, args.reference) if p})
    with open(args.output, 'x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    main()
