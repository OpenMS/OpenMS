"""Diagnostic (reported only): where does interference handling gain or lose IDs?
For each raw: accepted target precursors at entrapment FDP_c 1% and 2% (in-window winners by stage2_score), overlap of
base with R1 and R2, entrapment share in gained/lost sets; and, from R1 features (unmasked traces), the refinement
profile of base-accepted IDs (mostly true) versus decoy winners (all false)."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-single-scorer-20261109'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BASE = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
def accepted(d, t):
    cr = pd.read_csv(d + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'precursor_best', 'in_rt_window', 'stage2_score'])
    w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
    fdp = np.cumsum(w.e.values) * (1 + 1 / R) / np.arange(1, len(w) + 1); ok = np.where(fdp <= t)[0]
    return w.iloc[:ok.max() + 1 if len(ok) else 0]
for raw in ('raw01', 'raw02'):
    print(f'== {raw}')
    for t in (0.01, 0.02):
        base = accepted(BASE[raw], t)
        for arm in ('R1', 'R2'):
            a = accepted(f'{O}/{raw}/{arm}', t)
            kb, ka = set(base.precursor_id), set(a.precursor_id)
            gained, lost = a[~a.precursor_id.isin(kb)], base[~base.precursor_id.isin(ka)]
            print(f'  FDP_c {t:.0%} {arm}: base {len(kb)} arm {len(ka)} shared {len(kb & ka)} | gained {len(gained)} (entrapment {gained.e.sum()}) '
                  f'| lost {len(lost)} (entrapment {lost.e.sum()})')
    # Refinement profile from R1 features (traces unmasked), at each precursor's base winning window.
    f = pd.read_csv(f'{O}/{raw}/R1/features.tsv', sep='\t')
    rf = [x for x in f.columns if x.startswith('refine_')]
    base1 = accepted(BASE[raw], 0.01)
    cr = pd.read_csv(BASE[raw] + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'precursor_best', 'in_rt_window'])
    win = cr[(cr.precursor_best == 1) & (cr.in_rt_window == 1)]
    true_like = f[f.sample_id.isin(base1[base1.e == 0].sample_id)]
    decoys = f[f.sample_id.isin(win[win.is_decoy == 1].sample_id)]
    print(f'  refinement at base winners: accepted non-entrapment IDs at 1% (n={len(true_like)}) vs decoy winners (n={len(decoys)})')
    for name, g in (('accepted', true_like), ('decoy', decoys)):
        q = g.refine_quantitative.values; u = g.refine_usable.values
        hist = ' '.join(f'{k}:{(q == k).mean():.2f}' for k in range(7))
        print(f'    {name:8s} quantitative(r>=0.9) share by count {hist} | mean usable(r>=0.75) {u.mean():.2f} '
              f'| zero refine (unrefinable) {(g[rf].abs().sum(axis=1) == 0).mean():.2f} | mean r {g.refine_mean_r.mean():.3f}')
    masked = 6 - true_like.refine_usable.values
    print(f'    accepted IDs that R2 masks >=3 fragments: {(masked >= 3).mean():.2f}; >=1: {(masked >= 1).mean():.2f} '
          '(upper bound: windows with fewer than 6 library fragments count padding as masked)')
