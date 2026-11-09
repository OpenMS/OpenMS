"""PLAN.md reading for the single TOPAZ scorer (S0, S1 on raws 01 and 02). The metric is the native byion reading
(astral-raw02-byion-design-20261109/evaluate.py), with the score column `logit` (no stage 2), plus the rule check and
the reported-only items. All q-derived counts are decoy-count estimates."""
import os, re, sqlite3, sys
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
sys.path.insert(0, B + '/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled
O = B + '/astral-single-scorer-20261109'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BASE = {'raw01': (1029, 0.013, 0.065), 'raw02': (1298, 0.007, 0.074)}
def n_at(e, t):
    e = np.asarray(e); fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(fdp <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0
def usage(path):
    s = open(path).read() if os.path.exists(path) else ''
    wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', s); rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', s)
    parts = [float(x) for x in wall.group(1).split(':')] if wall else [float('nan')]
    minutes = sum(v * 60 ** (len(parts) - 1 - i) for i, v in enumerate(parts)) / 60
    return minutes, (int(rss.group(1)) / 1024 ** 2 if rss else float('nan'))
verdict = {}
for arm in ('S0', 'S1'):
    for raw in ('raw01', 'raw02'):
        d = f'{O}/{raw}/{arm}'
        if not os.path.exists(d + '/crossfit.tsv'):
            print(f'{raw} {arm}: missing'); verdict.setdefault(arm, []).append(False); continue
        cr = pd.read_csv(d + '/crossfit.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left')
        w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('logit', ascending=False)
        curve = [n_at(w.e.values, t) for t in (0.005, 0.01, 0.015, 0.02)]
        best, win = paired_pooled(cr, 'logit')
        t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
        fdp_q = {q: t[t.q <= q].e.sum() * (1 + 1 / R) / max((t.q <= q).sum(), 1) for q in (0.01, 0.05)}
        nq = {q: int((t.q <= q).sum()) for q in (0.01, 0.05)}
        # Null-symmetry R at the decoy-winner median per fold (cnn-topaz-design/PLAN.md).
        rs = []
        for f, g in win.groupby('fold'):
            med = g[g.is_decoy == 1].logit.median(); low = g[g.logit < med]
            rs.append((low.is_decoy.eq(0).sum() / max(low.is_decoy.sum(), 1)) / (g.is_decoy.eq(0).sum() / max(g.is_decoy.sum(), 1)))
        minutes, gb = usage(d + '/stderr.log')
        base = BASE[raw]
        checks = [curve[1] >= round(1.02 * base[0]), fdp_q[0.01] <= 0.015, fdp_q[0.05] <= base[2] + 0.005, minutes <= 35 and gb <= 30]
        verdict.setdefault(arm, []).extend(checks)
        print(f'{raw} {arm}: N@FDP_c 0.5/1/1.5/2% {curve} (base N@1% {base[0]}, need {round(1.02 * base[0])}) | '
              f'q<=0.01 {nq[0.01]} FDP_c {fdp_q[0.01]:.3f} | q<=0.05 {nq[0.05]} FDP_c {fdp_q[0.05]:.3f} (limit {base[2] + 0.005:.3f}) | '
              f'R {" ".join(f"{r:.3f}" for r in rs)} | {minutes:.1f} min {gb:.1f} GB | checks {checks}', flush=True)
        log = open(d + '/stdout.log').read()
        for line in re.findall(r'(?:score iteration=\d+ target_precursors_q01=\d+|Native CNN fold=.*|heuristic_features=.*)', log):
            print('   ', line[:300])
for arm, checks in verdict.items():
    print(f'{arm}: {"PASS" if checks and all(checks) else "FAIL"}')
