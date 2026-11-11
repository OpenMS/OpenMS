"""Stage-D deep dive, step 4: what information does the D11FWT score leave unused? (descriptive; labels define groups only)
Within score strata, compare raw-evidence measures (dd_evidence.feats, plus per-fragment interference counts) of
  D   stage-D right windows (true by DIA-NN, below the cut)
  DEC decoy winners
  ENT entrapment target winners (false targets, never in DIA-NN's human reference)
Strata: deciles of the D11FWT winner score over D + DEC + ENT restricted to score >= the decoy-winner q90. For each
measure, the stratified AUC (mean over strata of within-stratum AUC, weighted by D count) of D vs DEC and of ENT vs
DEC, oriented so 1 = target-like. A measure with D-vs-DEC far above 0.5 and ENT-vs-DEC near 0.5 is information the
score does not use and that false targets do not share: a usable cause. D-vs-DEC near 0.5 means D is indistinguishable
from decoys on that measure at the same score. ENT-vs-DEC far from 0.5 marks a target/decoy asymmetry (unusable).
Extra measures: n_interf = fragments with central signal whose profile correlates < 0.3 with the sum of the others;
excess = max over fragments of observed central-area share / library share (one fragment far above its library share
indicates interference)."""
import sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4')
from dd_evidence import feats, C0, C1
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
GOOD = {'n_sig': 1, 'n_peak': 1, 'log_area': 1, 'frac_core': 1, 'corr_med': 1, 'top1_share': -1, 'lib_cos': 1, 'rank_agree': 1, 'ppm_abs': -1,
        'ms1_m0': 1, 'ms1_corr': 1, 'iso_ratio': 1, 'n_sig12': 1, 'extra_sig': 1, 'n_interf': -1, 'excess': -1, 'bins_nz': 1, 'n_missing': -1}


def extra(rows, A):
  rows = np.sort(rows); x = np.asarray(A['xic'][rows], np.float64); lib = np.asarray(A['lib'][rows], np.float64)
  area = x[:, :, C0:C1].sum(2); tot = x.sum(1, keepdims=True); oth = tot - x
  xc = x - x.mean(2, keepdims=True); oc = oth - oth.mean(2, keepdims=True)
  r = (xc * oc).sum(2) / np.sqrt((xc * xc).sum(2) * (oc * oc).sum(2) + 1e-12)
  out = pd.DataFrame(index=rows)
  out['n_interf'] = ((area > 0) & (r < 0.3)).sum(1)
  sh = area / np.maximum(area.sum(1, keepdims=True), 1e-12); ls = lib / np.maximum(lib.sum(1, keepdims=True), 1e-12)
  out['excess'] = np.where(ls > 0, sh / np.maximum(ls, 1e-3), 0).max(1)
  return out


def auc(a, b):
  a, b = np.asarray(a, float), np.asarray(b, float); a, b = a[np.isfinite(a)], b[np.isfinite(b)]
  if len(a) == 0 or len(b) == 0: return np.nan
  r = pd.Series(np.r_[a, b]).rank().values; return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


for raw in sys.argv[1:] or ['raw01', 'raw02']:
  A = {k: np.load(f'{O}/{raw}/{k}.npy', mmap_mode='r') for k in ('xic', 'ppm', 'ms1', 'lib', 'xic12', 'orig12')}
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_D11FWT.npy').astype(np.float64)
  w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
  w = w[w.s > -1e20].copy()
  st = pd.read_csv(f'{O}/v3seq/misses_D11FWT_{raw}.tsv', sep='\t').set_index('precursor_id').stage
  w['grp'] = np.where(w.precursor_id.map(st) == 'D right window, below cut', 'D',
                      np.where(w.is_decoy == 1, 'DEC', np.where(w.precursor_id.map(ent) == 1, 'ENT', 'other')))
  q90 = w[w.grp == 'DEC'].s.quantile(0.9)
  sel = w[w.grp.isin(['D', 'DEC', 'ENT']) & (w.s >= q90)].copy()
  # Strata: 10 bins with edges at the D score quantiles, so each bin holds about 10 % of D and the decoys and
  # entrapment winners in the same narrow score range (score-matched comparison).
  edges = np.unique(np.r_[-np.inf, sel[sel.grp == 'D'].s.quantile(np.linspace(0.1, 0.9, 9)).values, np.inf])
  sel['stratum'] = pd.cut(sel.s, edges, labels=False)
  F = feats(sel.row.values, A).join(extra(sel.row.values, A))
  sel = sel.join(F, on='row')
  print(f'\n== {raw}: windows with score >= decoy q90: ' + ', '.join(f'{k} {int((sel.grp == k).sum())}' for k in ('D', 'DEC', 'ENT')))
  print('D / DEC / ENT per stratum:', sel.groupby(['stratum', 'grp']).size().unstack(fill_value=0).to_dict('index'))
  rows = []
  for m, g in GOOD.items():
    ad, ae, wt = [], [], []
    for _, t in sel.groupby('stratum'):
      d, dc, en = g * t[t.grp == 'D'][m], g * t[t.grp == 'DEC'][m], g * t[t.grp == 'ENT'][m]
      if len(d) < 5 or len(dc) < 5: continue
      ad.append(auc(d, dc)); ae.append(auc(en, dc) if len(en) >= 5 else np.nan); wt.append(len(d))
    wt = np.array(wt, float); ad, ae = np.array(ad), np.array(ae); ok = np.isfinite(ae)
    rows.append((m, np.average(ad, weights=wt), np.average(ae[ok], weights=wt[ok]) if ok.any() else np.nan))
  tab = pd.DataFrame(rows, columns=['measure', 'D_vs_DEC', 'ENT_vs_DEC']).set_index('measure')
  tab['usable_gap'] = tab.D_vs_DEC - tab.ENT_vs_DEC
  print('stratified AUC (1 = target-like):'); print(tab.sort_values('usable_gap', ascending=False).round(3).to_string())
  sel.drop(columns=['stratum']).to_csv(f'{O}/v4/dd_band_{raw}.tsv', sep='\t', index=False)
