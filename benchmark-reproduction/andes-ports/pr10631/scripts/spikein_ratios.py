#!/usr/bin/env python3
"""Ratio compression on PXD007683 a05058 (TMT11, yeast spiked into human). Channel groups and levels inferred from
the data and the project description (3-, 2- and 1.5-fold changes): yeast at 3x in 126-127C (ch0-2), 1.5x in
128N-129C (ch3-6), 1x in 130N-131C (ch7-10); human 1:1. Reports, per arm/tool and purity cut-off, the median
observed yeast fold change for each pair (expected 3, 2, 1.5) and the PSMs behind it.
usage: spikein_ratios.py <data_dir> <out_tsv>"""
import sys, csv, statistics, math
d, out = sys.argv[1:3]
G = [range(0, 3), range(3, 7), range(7, 11)]
PAIRS = [((0, 2), 3.0), ((0, 1), 2.0), ((1, 2), 1.5)]
ids = {r['spectrum_reference']: r for r in csv.DictReader(open(f'{d}/search/a05058.ids.tsv'), delimiter='\t') if r['target_decoy'] == 'target'}
def species(acc):
    s = {a.rsplit('_', 1)[-1] for a in acc.split(';') if a}
    return s.pop() if len(s) == 1 else None
w = csv.writer(open(out, 'w', newline=''), delimiter='\t')
w.writerow(['tool', 'arm', 'min_purity', 'yeast_psms', 'fc_3x_expected_3', 'fc_2x_expected_2', 'fc_1.5x_expected_1.5', 'human_psms', 'human_fc_3x_vs_1x_expected_1'])
for tool in ('IA', 'IW'):
    for arm in ('base', 'pr', 'fix'):
        rows = []
        for r in csv.DictReader(open(f'{d}/quant/{arm}/a05058.{tool}.tsv'), delimiter='\t'):
            if tool == 'IA':
                idr = ids.get(r['id_scan_id'] or r['scan_id'])
                if idr is None: continue
                r['accessions'] = idr['accessions']
            elif r['target_decoy'] != 'target':
                continue
            r['_pur'] = float(r['precursor_purity'] or 0)
            rows.append(r)
        for cut in (0.0, 0.5, 0.7, 0.9):
            sel = [r for r in rows if r['_pur'] >= cut]
            res = []
            for sp in ('YEAST', 'HUMAN'):
                g = [[statistics.mean(float(r[f'ch{i}']) for i in grp) for grp in G] for r in sel if species(r['accessions']) == sp]
                g = [x for x in g if min(x) > 0]
                res.append((len(g), [round(2 ** statistics.median(math.log2(x[a] / x[b]) for x in g), 3) for (a, b), _ in PAIRS]))
            (ny, fy), (nh, fh) = res
            w.writerow([tool, arm, cut, ny] + fy + [nh, fh[0]])
