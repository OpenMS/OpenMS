#!/usr/bin/env python3
"""Summary table (mean over Percolator seeds, with per-seed values) from evaluate.py JSONs.
Usage: summarize.py <eval dir> <out.tsv>"""
import glob, json, os, re, statistics, sys
from collections import defaultdict
rows = defaultdict(dict)
for f in glob.glob(os.path.join(sys.argv[1], '*.json')):
    d = json.load(open(f))
    m = re.match(r'(ko\d)_(\w+?)_s(\d+)_q([\d.]+)$', d['arm'])
    db, arm, seed, q = m.groups()
    c = d['counts']
    rows[(db, arm, q)][int(seed)] = {
        'target_psms': c.get('target', 0),
        'non_variant_psms': c.get('other', 0),
        'variant_psms': c.get('variant', 0),
        'variant_carried': c.get('variant_carried', 0),
        'variant_absent': c.get('variant_absent', 0),
        'variant_unlabelled': c.get('variant_unlabelled', 0),
        'variant_fdp': d['variant_fdp'],
        'knockoff_psms': c.get('knockoff', 0),
        'variants_carried': d['variants']['carried'],
        'variants_absent': d['variants']['absent'],
        'variants_fdp': d['variants']['fdp'],
    }
order = {'off': 0, 'class': 1, 'ko_onetier': 2, 'ko_tiers': 3, 'ko_onetier_fix': 4, 'ko_tiers_fix': 5}
keys = ['target_psms', 'non_variant_psms', 'variant_psms', 'variant_carried', 'variant_absent', 'variant_unlabelled',
        'variant_fdp', 'knockoff_psms', 'variants_carried', 'variants_absent', 'variants_fdp']
with open(sys.argv[2], 'w') as o:
    o.write('database\tarm\tq_cutoff\tseeds\t' + '\t'.join(keys) + '\t' + '\t'.join(k + '_per_seed' for k in keys) + '\n')
    for (db, arm, q) in sorted(rows, key=lambda k: (k[0], float(k[2]), order.get(k[1], 9))):
        r = rows[(db, arm, q)]
        seeds = sorted(r)
        means, per = [], []
        for k in keys:
            vals = [r[s][k] for s in seeds]
            num = [v for v in vals if v is not None]
            means.append(f'{statistics.mean(num):.4g}' if num else 'NA')
            per.append('/'.join('NA' if v is None else (f'{v:.3f}' if isinstance(v, float) else str(v)) for v in vals))
        o.write(f'{db}\t{arm}\t{q}\t{",".join(map(str, seeds))}\t' + '\t'.join(means) + '\t' + '\t'.join(per) + '\n')
print(open(sys.argv[2]).read())
