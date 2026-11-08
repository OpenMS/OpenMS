"""Raw 01 follow-up (descriptive): is there information left in the decision band, and how do our peak boundaries compare to OSW's?"""
import numpy as np, pandas as pd
exec(open('miss01.py').read().split("mD, wD")[0].replace("print(", "(lambda *a, **k: None)("))
mO, wO, missO, wrongO, norO = funnel('OSW', o)


def auc(a, b):
    a, b = a.dropna().values, b.dropna().values
    r = pd.Series(np.r_[a, b]).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


w = win[win.in_rt_window == 1]
hi = acc.stage2_score.min()
q05 = w[(w.is_decoy == 0) & (w.q <= 0.05)].stage2_score.min()
band = w[(w.stage2_score >= q05) & (w.stage2_score < hi)]
bt, bd = band[band.is_decoy == 0], band[band.is_decoy == 1]
print(f'decision band q 0.01..0.05: stage2 {q05:.1f}..{hi:.1f}, targets {len(bt)}, decoys {len(bd)}')
print('within-band target-vs-decoy AUC (incremental information over stage2_score):')
for k in cols + ['ms2_abs_ppm', 'mz', 'charge']:
    print(f'  {k:22s} {auc(bt[k], bd[k]):.2f}   target med {bt[k].median():.2f}  decoy med {bd[k].median():.2f}')
both = (bt.fragment_coelution > 0.6) & (bt.library_correlation > 0.6)
bothd = (bd.fragment_coelution > 0.6) & (bd.library_correlation > 0.6)
print(f'band rows with coelution>0.6 and libcorr>0.6: targets {int(both.sum())}/{len(bt)}, decoys {int(bothd.sum())}/{len(bd)}')

print('\npeak boundaries on OSW-accepted precursors, right winner (seconds):')
oo = pd.read_csv(OSW, sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'leftWidth', 'rightWidth'])
oo = oo[oo.run_name.str.endswith('_01')]
oo['mkey'] = oo.FullPeptideName + '/' + oo.Charge.astype(str)
oo = oo.merge(lib[['mkey', 'precursor_id']], on='mkey').assign(osw_width=lambda t: t.rightWidth - t.leftWidth)
j = wO[wO.right].merge(oo[['precursor_id', 'osw_width']], on='precursor_id')
print('  ours', j.peak_width_seconds.quantile([.25, .5, .75]).round(1).to_dict(), '| OSW', j.osw_width.quantile([.25, .5, .75]).round(1).to_dict(),
      '| ratio median', round((j.peak_width_seconds / j.osw_width).median(), 2))
print('  matched_scans per second of our width, median', round((j.matched_scans / j.peak_width_seconds).median(), 2))
