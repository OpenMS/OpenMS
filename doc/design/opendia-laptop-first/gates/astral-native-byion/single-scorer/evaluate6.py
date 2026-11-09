"""PLAN6.md reading: native CNN + LDA with refine features (R1) and with the trace mask (R2), raws 01 and 02.
Metric as astral-raw02-byion-design-20261109/evaluate.py (score stage2_score); rule check; mechanism AUCs (reported only).
All q-derived counts are decoy-count estimates."""
import os, re, sqlite3, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
sys.path.insert(0, B + '/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled
O = B + '/astral-single-scorer-20261109'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
RULE = {'raw01': (1050, 0.015, 0.070), 'raw02': (1324, 0.015, 0.079)}
BASE = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
def n_at(e, t):
    e = np.asarray(e); fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(fdp <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0
def usage(path):
    s = open(path).read() if os.path.exists(path) else ''
    wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', s); rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', s)
    parts = [float(x) for x in wall.group(1).split(':')] if wall else [float('nan')]
    return sum(v * 60 ** (len(parts) - 1 - i) for i, v in enumerate(parts)) / 60, (int(rss.group(1)) / 1024 ** 2 if rss else float('nan'))
def auc(e, d): return mannwhitneyu(e, d).statistic / (len(e) * len(d))
cmp_txt = open(O + '/CMP6.txt').read() if os.path.exists(O + '/CMP6.txt') else ''
verdict = {}
for arm in ('base', 'R1', 'R2'):
    for raw in ('raw01', 'raw02'):
        d = BASE[raw] if arm == 'base' else f'{O}/{raw}/{arm}'
        if not os.path.exists(d + '/crossfit.tsv'):
            print(f'{raw} {arm}: missing'); verdict.setdefault(arm, []).append(False); continue
        cr = pd.read_csv(d + '/crossfit.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left')
        w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
        curve = [n_at(w.e.values, t) for t in (0.005, 0.01, 0.015, 0.02)]
        best, _ = paired_pooled(cr, 'stage2_score')
        t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
        fdp = {q: t[t.q <= q].e.sum() * (1 + 1 / R) / max((t.q <= q).sum(), 1) for q in (0.01, 0.05)}
        nq = {q: int((t.q <= q).sum()) for q in (0.01, 0.05)}
        minutes, gb = usage(d + '/stderr.log')
        internal = arm != 'R1' or (f'{raw}/R1: crossfit cols 1-7,9 IDENTICAL' in cmp_txt and f'{raw}/R1: features cols 1-16 IDENTICAL' in cmp_txt)
        n1, f1, f5 = RULE[raw]
        ok = curve[1] >= n1 and fdp[0.01] <= f1 and fdp[0.05] <= f5 and gb <= 30 and minutes <= 35 and internal
        if arm != 'base': verdict.setdefault(arm, []).append(ok)
        print(f'{raw} {arm:4s} N at FDP_c 0.5/1/1.5/2%: {curve} | q<=0.01: {nq[0.01]} FDP_c {fdp[0.01]:.3f} | q<=0.05: {nq[0.05]} FDP_c {fdp[0.05]:.3f}'
              f' | {minutes:.1f} min {gb:.1f} GB | internal {internal} | rule {"PASS" if ok else "fail"}', flush=True)
        # Mechanism (reported only): entrapment vs decoy AUC over in-window winners by stage2_score.
        f = pd.read_csv(d + '/features.tsv', sep='\t').merge(cr[['sample_id', 'in_rt_window']], on='sample_id')
        f = f[f.in_rt_window == 1].sort_values(['stage2_score', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').merge(ent, on='precursor_id', how='left')
        e = f[(f.is_decoy == 0) & (f.e == 1)]; dd = f[f.is_decoy == 1]
        cols = [x for x in f.columns if x.startswith('refine_')] + ['transition_support', 'log_matched_scans', 'cnn_z', 'stage2_score']
        print('   AUC(entrapment > decoy winners): ' + ' '.join(f'{x}={auc(e[x].values, dd[x].values):.3f}' for x in cols if x in f), flush=True)
        log = open(d + '/stdout.log').read() if os.path.exists(d + '/stdout.log') else ''
        for line in re.findall(r'Experimental trace refinement:.*', log): print('  ', line)
        for line in re.findall(r'second stage fold=.*', log): print('  ', line[:700])
for arm, oks in verdict.items():
    print(f'{arm}: {"PASS" if all(oks) and len(oks) == 2 else "FAIL"} ({oks})')
