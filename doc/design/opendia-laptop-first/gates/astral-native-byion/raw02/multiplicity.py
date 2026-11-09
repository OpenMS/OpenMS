"""k10: candidate multiplicity (declared before running).
Hypothesis: a precursor's winner is the max over its candidate windows. If entrapment (false) targets carry more candidate
windows than decoys, e.g. because their fragments are more common and the search opens more peak groups, their winner
score is a max over more draws. That would inflate the false-target tail relative to decoys at the same null.
Measures (native byion arm, raws 01/02, windows in RT window):
  (a) candidate windows per precursor: mean and quantiles for human targets, entrapment, decoys; AUC of entrapment vs decoy.
  (b) the 5% tail ratio S_e/S_d (entrapment winners above the plain q<=0.05 cut over decoy winners above it, each as a share
      of its class) overall and within strata of candidate count (1, 2-3, 4-7, 8+).
  (c) label-free correction: weight each decoy winner by P_target(count stratum)/P_decoy(count stratum), with targets = all
      target precursors (no labels), then weighted paired TDC.
Reading: SUPPORTED if AUC(entrapment vs decoy count) >= 0.55 on both raws AND S_e/S_d within strata is closer to 1 than
overall (at least half the excess removed) on both. (c) PASSES if weighted q<=0.05 FDP_c <= 0.06 on both raws."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
ARMS = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BINS = [0, 1, 3, 7, 10 ** 9]; LABELS = ['1', '2-3', '4-7', '8+']
for run, arm in ARMS.items():
    f = pd.read_csv(arm + '/features.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score'])
    f = f.merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
    f['e'] = f.e.fillna(0).astype(int)
    s = f.stage2_score.values; dec = f.is_decoy.values == 1; fold = f.fold.values
    # pooled z as the tool
    w0 = f.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    z = s.copy()
    for k in np.unique(fold):
        g = w0[(w0.is_decoy == 1) & (w0.fold == k)].s; med = g.median(); z[fold == k] = (s[fold == k] - med) / (1.4826 * (g - med).abs().median())
    f['z'] = z
    n = f.groupby('precursor_id').size().rename('ncand')
    w = f.sort_values(['z', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').merge(n, on='precursor_id')
    w['cls'] = np.where(w.is_decoy == 1, 'decoy', np.where(w.e == 1, 'entrap', 'human'))
    w['stratum'] = pd.cut(w.ncand, BINS, labels=LABELS)
    print(f'== {run}: winners {len(w)}', flush=True)
    for cl, g in w.groupby('cls'):
        print(f'  {cl:7s} n={len(g):6d} ncand mean {g.ncand.mean():.2f} quantiles {g.ncand.quantile([.25, .5, .75, .95]).tolist()} strata {g.stratum.value_counts(normalize=True).reindex(LABELS).round(3).tolist()}')
    ed = w[w.cls != 'human']
    print(f'  AUC entrapment vs decoy ncand: {roc_auc_score((ed.cls == "entrap").astype(int), ed.ncand):.3f}')
    # plain paired TDC cut at q<=0.05
    d = w.sort_values(['z', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('z', ascending=False)
    dd = d.is_decoy.values == 1; q = np.minimum.accumulate(np.minimum(1, (np.cumsum(dd) + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]
    t5 = d.z.values[(~dd) & (q <= 0.05)].min()
    def ratio(g):
        e = g[g.cls == 'entrap']; dc = g[g.cls == 'decoy']
        if len(e) == 0 or len(dc) == 0: return np.nan, 0, 0
        ke, kd = int((e.z >= t5).sum()), int((dc.z >= t5).sum())
        return (ke / len(e)) / max(kd / len(dc), 1e-12), ke, kd
    r0, ke, kd = ratio(w); print(f'  overall S_e/S_d at q<=0.05 cut: {r0:.2f} (entrap {ke}, decoy {kd})')
    for st in LABELS:
        r, ke, kd = ratio(w[w.stratum == st]); print(f'    stratum {st:4s}: S_e/S_d {r:.2f} (entrap {ke}, decoy {kd})')
    # (c) label-free reweighting of decoy winners by target/decoy stratum shares
    pt = w[w.is_decoy == 0].stratum.value_counts(normalize=True); pd_ = w[w.is_decoy == 1].stratum.value_counts(normalize=True)
    wt = (pt / pd_).to_dict(); d['wd'] = np.where(d.is_decoy == 1, d.stratum.map(wt).astype(float), 0.0)
    Dn = np.cumsum(d.wd.values); qw = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]
    out = []
    for nm, qq_ in (('plain', q), ('count-weighted', qw)):
        for qq in (0.01, 0.05):
            a = d[(~dd) & (qq_ <= qq)]; out.append(f'{nm} q<={qq}: {len(a)} FDP_c {a.e.sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
    print('  weights', {k: round(v, 3) for k, v in wt.items()}); print('  ' + ' | '.join(out), flush=True)
