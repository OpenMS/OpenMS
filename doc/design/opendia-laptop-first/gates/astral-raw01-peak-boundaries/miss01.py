"""Raw 01 miss analysis (descriptive, diagnostic only; nothing fitted or selected; entrapment labels not used).
References: DIA-NN 1% exact library matches and same-library OpenSWATH+Percolator 1% (run 01).
Our result: frozen grid2d CNN + stage-2 search set (stage2-search-design/search)."""
import numpy as np, pandas as pd, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
LIB = B + '/astral-entrapment-20261107/lib/library.pqp'
OSW = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv'
pd.set_option('display.width', 250)

lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet')
lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
dn = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
dn = dn[dn.Run.str.endswith('_01')][['precursor_id', 'RT']].assign(ref_rt=lambda x: x.RT * 60)[['precursor_id', 'ref_rt']]
o = pd.read_csv(OSW, sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'Intensity', 'd_score'])
o = o[o.run_name.str.endswith('_01')]
o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
o = o.merge(lib[['mkey', 'precursor_id']], on='mkey').rename(columns={'RT': 'ref_rt', 'Intensity': 'osw_intensity'})
print('references run 01: DIA-NN exact', len(dn), '| OSW+Percolator in library', len(o), 'both', len(set(dn.precursor_id) & set(o.precursor_id)))

c = pd.read_csv(D + '/stage2-search-design/search/crossfit.tsv', sep='\t')
c['q'] = pd.to_numeric(c.q_value, errors='coerce')
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t').rename(columns={'candidate_row': 'sample_id'}).drop(columns=['precursor_id', 'is_decoy'])
x = pd.read_csv(D + '/stage2-design/features.tsv', sep='\t', usecols=['sample_id', 'fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope'])
f2 = pd.read_csv(D + '/stage2-search-design/search/features.tsv', sep='\t', usecols=['sample_id', 'cnn_z'])
c = c.merge(s, on='sample_id').merge(x, on='sample_id').merge(f2, on='sample_id')
con = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
pr = pd.read_sql('select ID precursor_id, CHARGE charge, PRECURSOR_MZ mz from PRECURSOR', con); con.close()
c = c.merge(pr, on='precursor_id')
ncand = c.groupby('precursor_id').size().rename('n_candidates')
c = c.merge(ncand, on='precursor_id')
c['log_apex'] = np.log10(c.apex_intensity.clip(lower=1))
win = c[(c.precursor_best == 1)]
acc = win[(win.q <= 0.01) & (win.is_decoy == 0)]
dec = win[win.is_decoy == 1]
top_dec = dec.nlargest(int(0.01 * len(dec)), 'stage2_score')
print('our accepted at q<=0.01', len(acc), '| decoy winners', len(dec))

cols = ['cnn_z', 'stage2_score', 'transition_support', 'matched_scans', 'peak_width_seconds', 'abs_drt', 'log_apex', 'n_candidates',
        'fragment_coelution', 'library_correlation', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope', 'ms2_abs_ppm']
c['abs_drt'] = c.apex_rt_error_seconds.abs()
win = c[c.precursor_best == 1]; acc = win[(win.q <= 0.01) & (win.is_decoy == 0)]
dec = win[win.is_decoy == 1]; top_dec = dec.nlargest(int(0.01 * len(dec)), 'stage2_score')


def funnel(name, ref):
    m = c[c.is_decoy == 0].merge(ref[['precursor_id', 'ref_rt']], on='precursor_id')
    m['right'] = (m.rt_apex - m.ref_rt).abs() <= 10
    g = m.groupby('precursor_id')
    w = m[m.precursor_best == 1]
    rows = [('reference precursors', len(ref)),
            ('any candidate window', g.ngroups),
            ('a candidate within 10 s of ref RT', int(g.right.any().sum())),
            ('winner at ref RT', int(w.right.sum())),
            ('accepted q<=0.01', int((w.q <= 0.01).sum())),
            ('accepted q<=0.01 and at ref RT', int(((w.q <= 0.01) & w.right).sum())),
            ('accepted q<=0.05', int((w.q <= 0.05).sum()))]
    print(f'\n== funnel vs {name}')
    for k, v in rows:
        print(f'  {k:40s} {v:6d}  ({v / len(ref):.0%})')
    right_miss = w[w.right & ~(w.q <= 0.01)]
    wrong = w[~w.right]
    # A missed right-peak precursor whose right candidate exists but lost the winner slot
    lost = m[m.right & (m.precursor_best == 0)].precursor_id.unique()
    lost = np.setdiff1d(lost, w[w.right].precursor_id)
    no_right = np.setdiff1d(ref.precursor_id, m[m.right].precursor_id)
    print(f'  loss: no candidate at ref RT {len(no_right)}, right candidate lost winner slot {len(lost)}, right winner not accepted {len(right_miss)}')
    return m, w, right_miss, wrong, no_right


mD, wD, missD, wrongD, norD = funnel('DIA-NN 1% exact', dn)
mO, wO, missO, wrongO, norO = funnel('OSW+Percolator 1%', o)

print('\n== feature medians (winners)')
groups = {'accepted (ours)': acc, 'OSW-acc, right peak, missed': missO, 'DIA-NN, right peak, missed': missD,
          'OSW-acc, wrong winner': wrongO, 'top 1% decoy winners': top_dec, 'all decoy winners': dec}
print(pd.DataFrame({k: g[cols].median() for k, g in groups.items()}).round(2).T.assign(n=[len(g) for g in groups.values()]).to_string())

print('\n== missed right-peak (OSW ref): q-value distribution of the winner')
print(missO.q.describe(percentiles=[.1, .25, .5, .75, .9]).round(3).to_dict())
print('stage2 z of 1% cut (lowest accepted target):', round(acc.stage2_score.min(), 2), '| missed OSW right-peak median', round(missO.stage2_score.median(), 2))

print('\n== where they differ: share of group by transition_support (of 6)')
for k in ['accepted (ours)', 'OSW-acc, right peak, missed', 'top 1% decoy winners']:
    print(f'  {k:30s}', groups[k].transition_support.value_counts(normalize=True).sort_index().round(2).to_dict())

print('\n== charge breakdown: OSW-accepted vs ours accepted')
oc = o.merge(pr, on='precursor_id').charge.value_counts().sort_index().to_dict()
ac = acc.charge.value_counts().sort_index().to_dict()
mc = missO.charge.value_counts().sort_index().to_dict()
print('  OSW acc', oc, '\n  ours acc', ac, '\n  OSW right-peak missed', mc)

print('\n== OSW intensity tertiles: recall of OSW-accepted precursors by us (q<=0.01)')
oo = o.merge(wO[['precursor_id', 'q', 'right']], on='precursor_id', how='left')
oo['ours'] = (oo.q <= 0.01)
oo['tert'] = pd.qcut(oo.osw_intensity.rank(method='first'), 3, labels=['low', 'mid', 'high'])
print(oo.groupby('tert', observed=True).ours.mean().round(2).to_dict(), '| overall', round(oo.ours.mean(), 2))
oo['dbin'] = pd.qcut(oo.d_score.rank(method='first'), 4, labels=['q1(weak)', 'q2', 'q3', 'q4(strong)'])
print('by OSW d_score quartile', oo.groupby('dbin', observed=True).ours.mean().round(2).to_dict())

print('\n== wrong winners (OSW ref): |our winner RT - ref RT| quantiles (s), and was right candidate present?')
print((wrongO.rt_apex - wrongO.ref_rt).abs().quantile([.25, .5, .75, .9]).round(0).to_dict())
rc = mO[mO.right].groupby('precursor_id').stage2_score.max()
ww = wrongO.set_index('precursor_id')
both = ww.join(rc.rename('right_score'), how='inner')
print('  right candidate present for', len(both), 'of', len(ww), '| winner-minus-right stage2 median', round((both.stage2_score - both.right_score).median(), 2),
      '| right candidate in RT window share', round(mO[mO.right & mO.precursor_id.isin(both.index)].groupby('precursor_id').in_rt_window.max().mean(), 2))

print('\n== no candidate at ref RT (OSW ref): what is the closest candidate?')
nr = mO[mO.precursor_id.isin(norO)].assign(d=lambda t: (t.rt_apex - t.ref_rt).abs()).groupby('precursor_id').d.min()
print('  with any candidate', len(nr), 'of', len(norO), '| closest |dRT| quantiles', nr.quantile([.25, .5, .75]).round(0).to_dict())
lr = o[o.precursor_id.isin(norO)].merge(c[['precursor_id']].drop_duplicates(), on='precursor_id', how='left', indicator=True)
print('  ref RT - library-calibrated window? share of these refs whose lib dRT>75 s not computable here; n no candidate at all', int((lr._merge == 'left_only').sum()))

print('\n== separability inside the decision band: OSW right-peak missed vs top 1% decoys (AUC per feature)')
def auc(a, b):
    a, b = a.dropna().values, b.dropna().values
    r = pd.Series(np.r_[a, b]).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
for k in cols:
    print(f'  {k:22s} {auc(missO[k], top_dec[k]):.2f}')
