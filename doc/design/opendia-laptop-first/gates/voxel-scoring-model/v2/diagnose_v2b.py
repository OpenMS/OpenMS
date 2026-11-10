"""Voxel v2 diagnosis, part 2 (after the frozen result, no training, no verdict change).
Part 1 showed the lost base hits score above 99 % of the arm's decoy winners. So the N@1 % cut must be set by
entrapment target winners that the arm puts above them. This script measures, per arm and BASE:
  - the top of the target-winner list: entrapment count in the top 500 / 1000 / 1500 target winners,
    and decoy winners above the 1000th target (the decoy-count estimate at that depth);
  - what those top entrapment winners look like (heuristic medians), and how the base ranks them;
  - whether the library input alone separates targets from decoys (a decoy-tell: entrapment shares it).
Entrapment labels are read only for these evaluation numbers."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv',
         'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e

def winners(pop, s):
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    return w.sort_values('s', ascending=False)

for raw in sys.argv[1:] or ['raw01', 'raw02']:
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    H = pd.DataFrame(np.load(f'{O}/{raw}/heuristics.npy'), columns=FEATS)
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    print(f'\n===== {raw}')
    wb = winners(pop, base); tb = wb[wb.is_decoy == 0].reset_index(drop=True)
    base_rank = pd.Series(np.arange(len(tb)), index=tb.precursor_id)
    for arm in ('BASE', 'V0', 'F0', 'S0', 'FS'):
        s = base if arm == 'BASE' else np.load(f'{O}/{raw}/scores_{arm}.npy')
        w = winners(pop, s); t = w[w.is_decoy == 0].reset_index(drop=True); d = w[w.is_decoy == 1]
        t['e'] = t.precursor_id.map(ent).values
        line = f'{arm:4s} entrapment in top 500/1000/1500 targets: ' + '/'.join(str(int(t.e[:k].sum())) for k in (500, 1000, 1500))
        line += f' | decoy winners above target #1000: {(d.s > t.s[999]).sum()}'
        top_e = t[:1500][t.e[:1500] == 1]
        line += f' | top-1500 entrapment: base rank median {np.median(base_rank.loc[top_e.precursor_id]) if len(top_e) else np.nan:.0f}'
        print(line)
        if arm in ('V0', 'FS') and len(top_e):
            hm = H.iloc[top_e.row.values].median(); ht = H.iloc[t[:1000][t.e[:1000] == 0].row.values].median()
            print('     medians      ' + ' '.join(f'{f[:12]:>12s}' for f in FEATS))
            print('     top non-ent  ' + ' '.join(f'{v:12.3f}' for v in ht.values))
            print('     top entrap   ' + ' '.join(f'{v:12.3f}' for v in hm.values))
    # decoy tell from library input alone: logistic regression on sorted library intensities of the winners' precursors
    lib = np.load(f'{O}/{raw}/lib.npy', mmap_mode='r')
    w = winners(pop, base); rows = w.row.values; y = (w.is_decoy.values == 0).astype(int)
    Xl = np.sort(np.asarray(lib[rows]), 1)[:, ::-1]; f = w.fold.values
    p = np.zeros(len(y))
    for k in range(3):
        m = LogisticRegression(max_iter=1000).fit(Xl[f != k], y[f != k]); p[f == k] = m.predict_proba(Xl[f == k])[:, 1]
    e = w.precursor_id.map(ent).values
    print(f'library-only target-vs-decoy AUC {roc_auc_score(y, p):.3f}; entrapment-vs-decoy {roc_auc_score(np.r_[np.ones((e == 1).sum()), np.zeros((y == 0).sum())], np.r_[p[e == 1], p[y == 0]]):.3f}')
