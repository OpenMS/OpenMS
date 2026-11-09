"""Diagnostic (no decision rule; reported only): for N1 and L1, per-fold share of target winners at pooled q<=0.01/0.05,
and FDP_c at q<=0.01/0.05 under fold-wise TDC (each fold its own (D+1)/T) instead of pooled. Tests whether the q
miscalibration is a per-fold score-scale problem (fixable without a second learner) or a within-fold tail problem."""
import sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
sys.path.insert(0, B + '/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
def tdc(best):
    best = best.sort_values('z', ascending=False); d = best.is_decoy.values
    q = np.minimum(1, (np.cumsum(d) + 1) / np.maximum(np.cumsum(1 - d), 1))
    return best.assign(q=np.minimum.accumulate(q[::-1])[::-1])
def fdp(t, qq):
    s = t[t.q <= qq]; return len(s), s.e.sum() * (1 + 1 / R) / max(len(s), 1)
for arm in ('N1', 'L1'):
    for raw in ('raw01', 'raw02'):
        cr = pd.read_csv(f'{raw}/{arm}/crossfit.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left')
        best, w = paired_pooled(cr, 'logit')
        t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
        share = {qq: t[t.q <= qq].fold.value_counts().sort_index().tolist() for qq in (0.01, 0.05)}
        # fold-wise paired TDC on the raw logit within each fold
        rows = []
        for f, g in w.groupby('fold'):
            g = g.assign(z=g.logit, key=g.pair_key).sort_values(['z', 'is_decoy'], ascending=[False, False]).drop_duplicates('key')
            rows.append(tdc(g))
        fw = pd.concat(rows); tf = fw[fw.is_decoy == 0].merge(ent, on='precursor_id')
        print(f'{raw} {arm}: pooled per-fold targets q<=0.01 {share[0.01]} q<=0.05 {share[0.05]} | fold-wise TDC '
              f'q<=0.01 N {fdp(tf, 0.01)[0]} FDP_c {fdp(tf, 0.01)[1]:.3f} | q<=0.05 N {fdp(tf, 0.05)[0]} FDP_c {fdp(tf, 0.05)[1]:.3f}', flush=True)
