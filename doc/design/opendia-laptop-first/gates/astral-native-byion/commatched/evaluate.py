"""k11 reading (PLAN.md): the base evaluate.py metric on the commonness-matched library. Entrapment labels are joined here only.
Usage: evaluate.py raw01|raw02"""
import sys, sqlite3
import numpy as np, pandas as pd
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
sys.path.insert(0, D + '/cnn-topaz-design')
from evaluate import paired_pooled
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
O = B + '/astral-commatched-decoys-20261109'
BASE = {'raw01': (1029, 0.013, 0.065, D + '/stage2-byion-design/byion'), 'raw02': (1298, 0.007, 0.074, B + '/astral-raw02-byion-design-20261109/byion')}
def labels(lib):
    c = sqlite3.connect('file:' + lib + '?mode=ro', uri=True)
    p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
    p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
    ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
    return p, ent
def read(path, lib):
    p, ent = labels(lib); cr = pd.read_csv(path + '/crossfit.tsv', sep='\t')
    w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
    e = w.e.values; fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1)
    n = [int(np.where(fdp <= t)[0].max() + 1) if (fdp <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
    best, _ = paired_pooled(cr.merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left'), 'stage2_score')
    t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
    qf = {q: (int((t.q <= q).sum()), t[t.q <= q].e.sum() * (1 + 1 / R) / max((t.q <= q).sum(), 1)) for q in (0.01, 0.05)}
    return n, qf
run = sys.argv[1]; n1b, q1b, q5b, bpath = BASE[run]
nb, qb = read(bpath, B + '/astral-entrapment-20261107/lib/library.pqp')
nn, qn = read(f'{O}/{run}/crossfit', O + '/library.pqp')
for name, n, q in (('base (recomputed)', nb, qb), ('commonness-matched', nn, qn)):
    print(f'{run} {name:20s} N at FDP_c 0.5/1/1.5/2%: {n} | q<=0.01: {q[0.01][0]} FDP_c {q[0.01][1]:.3f} | q<=0.05: {q[0.05][0]} FDP_c {q[0.05][1]:.3f}')
ok5 = qn[0.05][1] <= 0.06; ok1 = qn[0.01][1] <= 0.015; okn = nn[1] >= 0.95 * nb[1]; drop = qb[0.05][1] - qn[0.05][1]
print(f'{run} conditions: q<=0.05 FDP_c<=0.06 {ok5} | q<=0.01 FDP_c<=0.015 {ok1} | N>=95% base {okn} | 5% drop {drop:.3f}')
