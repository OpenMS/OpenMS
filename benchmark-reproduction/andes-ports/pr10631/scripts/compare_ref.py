#!/usr/bin/env python3
"""Agreement of IsobaricWorkflow precursor purities (single survey scan, -extraction:purity_interpolation false) with the
andes-style reference purity (ref_purity.py), per arm, split by whether the selected ion lies inside the isolation window.
usage: compare_ref.py <data_dir> <out_tsv> run1 run2 ...   (reads quant/<arm>_noint/<run>.IW.tsv)"""
import sys, csv, statistics, os
d, out = sys.argv[1:3]
w = csv.writer(open(out, 'w', newline=''), delimiter='\t')
w.writerow(['run', 'subset', 'arm', 'psms', 'median_abs_diff_to_ref', 'within_0.05_of_ref', 'above_ref_by_0.2', 'purity_ge_0.9', 'ref_purity_ge_0.9'])
for run in sys.argv[3:]:
    sub = 'pxd017803' if os.path.exists(f'{d}/pxd017803/{run}.scans.tsv') else 'arch'
    sc = {r['native_id']: r for r in csv.DictReader(open(f'{d}/{sub}/{run}.scans.tsv'), delimiter='\t')}
    ref = {r['ms2_id']: float(r['ref_purity']) for r in csv.DictReader(open(f'{d}/{sub}/{run}.refpurity.tsv'), delimiter='\t') if r['ref_purity'] != ''}
    arms = {}
    for arm in os.environ.get('ARMS', 'base pr fix').split():
        a = {}
        for r in csv.DictReader(open(f'{d}/quant/{arm}_noint/{run}.IW.tsv'), delimiter='\t'):
            if r['target_decoy'] != 'target': continue
            k = r['scan_id']; k = sc[k]['spectrum_ref'] if sc[k]['ms_level'] == '3' else k
            a[k] = float(r['precursor_purity'] or 0)
        arms[arm] = a
    def inside(s):
        t, x = float(s['iso_target']), float(s['selected_mz']); return t - float(s['iso_lo']) <= x <= t + float(s['iso_hi'])
    for name, cond in (('selected ion inside window', inside), ('selected ion outside window', lambda s: not inside(s))):
        ks = [k for k in next(iter(arms.values())) if all(k in arms[a] for a in arms) and k in ref and cond(sc[k])]
        for arm, a in arms.items():
            diff = [a[k] - ref[k] for k in ks]
            w.writerow([run, name, arm, len(ks), round(statistics.median(abs(x) for x in diff), 4) if ks else '',
                        sum(abs(x) <= 0.05 for x in diff), sum(x > 0.2 for x in diff), sum(a[k] >= 0.9 for k in ks), sum(ref[k] >= 0.9 for k in ks)])
