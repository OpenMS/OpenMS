"""Follow-up on the 626 right-peak-but-rejected DIA-NN IDs (evaluation only): how far below the cut, what signal do they
have, and what evidence DIA-NN reported for them. Signal: OpenDIA's apex/area from search.features.tsv at the winning
window; DIA-NN: its Evidence, CScore, Ms1 area, fragment count used. Also MS1: does DIA-NN's MS1 profile correlation help?"""
import numpy as np, pandas as pd, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/gap-20261114/armE_raw01'
exec(open(B + '/gap-20261114/funnel.py').read().split("# Why stage 1 loses")[0].split('print(\'\\nFunnel')[0])
ref_all = pd.read_parquet(B + '/diann-pairlib-20261113/report_raw01.parquet')
ref_all['key'] = ref_all['Modified.Sequence'] + '/' + ref_all['Precursor.Charge'].astype(str)
cols = [c for c in ('Evidence', 'CScore', 'Ms1.Area', 'Ms1.Profile.Corr', 'Precursor.Quantity', 'Q.Value', 'PEP') if c in ref_all]
ra = ref_all.drop_duplicates('key').set_index('key')[cols]
human = human.set_index('precursor_id')
near_is_winner_ = near_is_winner
rej = sorted(near_is_winner_ - accepted); acc = sorted(near_is_winner_ & accepted)
sf = pd.read_csv(f'{D}/search.features.tsv', sep='\t', usecols=['candidate_row', 'precursor_id', 'apex_intensity', 'area', 'coelution', 'peaks_in_window'])
sf = sf.set_index('candidate_row')
wr = w.set_index('precursor_id')
def describe(ids, name):
    sid = wr.loc[ids, 'sample_id'].values
    f = sf.loc[sid]
    keys = human.loc[ids, 'key'].values
    r = ra.loc[keys]
    print(f'{name:9s} n {len(ids):4d} | OpenDIA score p50 {wr.loc[ids, "s"].median():6.2f} | apex int p50 {f.apex_intensity.median():9.0f} area p50 {f.area.median():10.0f} '
          f'peaks_in_window p50 {f.peaks_in_window.median():.0f} | DIA-NN ' + ' '.join(f'{c} {r[c].median():.3g}' for c in cols))
describe(acc, 'accepted')
describe(rej, 'rejected')
# Intensity strata: acceptance rate by DIA-NN quantity decile
both = sorted(near_is_winner_)
q = ra.loc[human.loc[both, 'key'].values, 'Precursor.Quantity'].values
dec_ = pd.qcut(q, 5, labels=False)
acc_flag = np.isin(both, acc)
print('acceptance of right-peak winners by DIA-NN quantity quintile (low->high):', [round(acc_flag[dec_ == i].mean(), 3) for i in range(5)])
# Score distance to cut
s_rej = wr.loc[rej, 's'].values
print('rejected score quantiles (cut %.2f):' % cut, np.quantile(s_rej, [.1, .25, .5, .75, .9]).round(2))
# Rank of rejected among all target winners: how many would we need to accept
rank = tw.reset_index(drop=True).precursor_id.reset_index().set_index('precursor_id')['index']
print('rejected rank among target winners p50 / p90:', int(np.median(rank.loc[rej])), int(np.quantile(rank.loc[rej], .9)), 'of', len(tw))
