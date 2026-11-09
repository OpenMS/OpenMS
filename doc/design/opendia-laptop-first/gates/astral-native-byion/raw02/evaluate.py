"""PLAN.md reading for the raw 01 native byion run: entrapment-matched N at FDP_c (matched_fdp.py metric) and pooled paired
q (the tool's TDC) for search+boundary (reference) and search+boundary+byion, plus FDP_c at q <= 0.01 / 0.05."""
import sys, re, sqlite3
import numpy as np, pandas as pd
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
sys.path.insert(0, D + '/cnn-topaz-design')
from evaluate import paired_pooled
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
LIB = B + '/astral-entrapment-20261107/lib/library.pqp'
c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
def fdp_curve(e):
    e = np.asarray(e); n = np.arange(1, len(e) + 1); fdp = np.cumsum(e) * (1 + 1 / R) / n; out = []
    for t in (0.005, 0.01, 0.015, 0.02):
        ok = np.where(fdp <= t)[0]; out.append(int(ok.max() + 1) if len(ok) else 0)
    return out
O = B + '/astral-raw02-byion-design-20261109'  # raw 02 copy; metric unchanged
for name, path in [('search+boundary', O + '/boundary'), ('search+boundary+byion', O + '/byion')]:
    cr = pd.read_csv(path + '/crossfit.tsv', sep='\t')
    w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
    best, _ = paired_pooled(cr.merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left'), 'stage2_score')
    t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
    qs = ' | '.join(f'q<={q}: {int((t.q <= q).sum())} FDP_c {t[t.q <= q].e.sum() * (1 + 1 / R) / max((t.q <= q).sum(), 1):.3f}' for q in (0.01, 0.05))
    print(f'{name:22s} N at FDP_c 0.5/1/1.5/2%: {fdp_curve(w.e.values)} | {qs}', flush=True)
log = open(O + '/byion/stdout.log').read()
for line in re.findall(r'second stage fold=.*', log):
    print(' ', line[:600])
