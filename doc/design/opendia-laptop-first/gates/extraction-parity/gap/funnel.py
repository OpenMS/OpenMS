"""Where OpenDIA loses DIA-NN's IDs (evaluation only; DIA-NN is a reference, nothing is fitted to it).
Run: arm E (C++ predict_decoys library), raw01, D11FDWT score. Reference: DIA-NN 2.2 on the same targets with its own
decoys (diann-pairlib-20261113, Q.Value <= 0.01, human targets). Population windows (population.tsv) are OpenDIA's
candidate peaks per precursor (up to 5) for pairs with transition_support >= 3.
Funnel per reference precursor:
  1. in the scored population at all (else: search/population loss);
  2. has a candidate window with apex within 1 FWHM-ish (6 s) of DIA-NN's RT (else: no candidate at DIA-NN's peak);
  3. that window is OpenDIA's precursor winner (else: wrong-peak choice);
  4. the winner is accepted at entrapment FDP_c 1 % (the rank cut N@FDP_c 1 %; else: scoring loss).
Then compares the 20 scorer inputs of right-peak-but-rejected against accepted and against decoy winners near the cut.
Usage: funnel.py (paths below)."""
import json, os, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
RAW = os.environ.get('RAW', 'raw01')  # raw02: DIA-NN reference is diann-decoys-20261114/own/raw02 (same targets, own decoys)
D = B + f'/gap-20261114/armE_{RAW}'
REF = {'raw01': B + '/diann-pairlib-20261113/report_raw01.parquet', 'raw02': B + '/diann-decoys-20261114/own/raw02/report.parquet'}[RAW]
TOL = 6.0
NAMES = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y', 'f12_support',
         'f12_coelution', 'f12_libcos', 'f12_extra_support', 'f12_extra_coelution', 'f12_extra_share', 'd_max', 'd6_max', 'd_lib']
c = sqlite3.connect(f'file:{D}/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, p.DECOY decoy, p.CHARGE z, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pre['key'] = pre.ms + '/' + pre.z.astype(str)
pre['ent'] = pre.precursor_id.map(ent)
tg = pre[pre.decoy == 0].set_index('key')
k = 1 + 1 / ((tg.ent == 1).sum() / (tg.ent == 0).sum())

ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str)
ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key')
ref = ref[ref.key.isin(tg.index)].copy()
ref['precursor_id'] = tg.loc[ref.key, 'precursor_id'].values
ref['ent'] = tg.loc[ref.key, 'ent'].values
ref['rt_s'] = ref.RT * 60
human = ref[ref.ent == 0]

pop = pd.read_csv(f'{D}/population.tsv', sep='\t')
X = np.c_[np.load(f'{D}/heuristics.npy'), np.load(f'{D}/feat12.npy'), np.load(f'{D}/featD.npy')]
s = np.load(f'{D}/scores_D11FDWT.npy').astype(float)
pop['s'] = s; pop['row'] = np.arange(len(pop))
# Winners and the N@FDP_c 1 % cut, as evaluate_v2.py (standardised by decoy winners per fold).
w = pop.sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
    w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
tw = w[w.is_decoy == 0].sort_values('s', ascending=False).copy()
tw['ent'] = tw.precursor_id.map(ent)
fdp = np.cumsum(tw.ent.values) * k / np.arange(1, len(tw) + 1)
n1 = int(np.where(fdp <= 0.01)[0].max() + 1)
cut = tw.s.iloc[n1 - 1]
accepted = set(tw.precursor_id.iloc[:n1])
print(f'OpenDIA arm E {RAW}: target winners {len(tw)}, N@FDP_c 1 % {n1} (score cut {cut:.2f})')
print(f'DIA-NN q<=0.01 in library: {len(ref)} ({int(ref.ent.sum())} entrapment); human {len(human)}')

pp = pop[pop.precursor_id.isin(human.precursor_id)].merge(human[['precursor_id', 'rt_s']], on='precursor_id')
pp['d'] = (pp.rt_apex - pp.rt_s).abs()
in_pop = set(pp.precursor_id)
near = pp[pp.d <= TOL]
has_near = set(near.precursor_id)
win_rows = w.set_index('precursor_id').row
near_is_winner = set(near[near.row.values == win_rows.loc[near.precursor_id].values].precursor_id)
acc_right = near_is_winner & accepted
acc_wrong = (in_pop & accepted) - near_is_winner
H = set(human.precursor_id)
print('\nFunnel (human DIA-NN IDs):')
print(f'  1. in OpenDIA scored population            {len(in_pop):5d} / {len(H)}   lost: {len(H - in_pop)}')
print(f'  2. a candidate peak within {TOL:.0f} s of DIA-NN RT {len(has_near):5d}           lost: {len(in_pop - has_near)}')
print(f'  3. that peak is the winner                  {len(near_is_winner):5d}           lost: {len(has_near - near_is_winner)}')
print(f'  4. winner accepted at FDP_c 1 %             {len(acc_right):5d}           lost: {len(near_is_winner - acc_right)}')
print(f'  (accepted with a different peak than DIA-NN: {len(acc_wrong)})')
acc_all = set(tw.precursor_id.iloc[:n1])
print(f'  OpenDIA accepted not in DIA-NN q<=0.01: {len(acc_all - set(ref.precursor_id))}')

if not os.path.exists(f'{D}/search.features.tsv'): raise SystemExit
# Why stage 1 loses: are the missing ones absent from the search table or below the population support rule?
sf = pd.read_csv(f'{D}/search.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex', 'transition_support'])
miss1 = H - in_pop
sm = sf[sf.precursor_id.isin(miss1)].merge(human[['precursor_id', 'rt_s']], on='precursor_id')
print(f'\nStage-1 losses ({len(miss1)}): with any search candidate {sm.precursor_id.nunique()}, '
      f'candidate within {TOL:.0f} s {sm[(sm.rt_apex - sm.rt_s).abs() <= TOL].precursor_id.nunique()}, '
      f'max transition_support of those near: {sm[(sm.rt_apex - sm.rt_s).abs() <= TOL].groupby("precursor_id").transition_support.max().value_counts().sort_index().to_dict()}')
cal = json.load(open(f'{D}/cal.json'))['rt']
lrt = tg.set_index('precursor_id').lrt
pred = np.interp(lrt.loc[list(miss1)].values, cal['knots_library_rt'], cal['knots_run_rt_s'])
err = np.abs(pred - human.set_index('precursor_id').loc[list(miss1), 'rt_s'].values)
print(f'  their |DIA-NN RT - map| p50 {np.median(err):.1f} s, outside the {cal["tolerance_half_width_s"]:.1f} s window: {np.mean(err > cal["tolerance_half_width_s"]):.3f}')

# Scoring loss: features of right-peak-rejected vs accepted vs decoys near the cut.
rej = near_is_winner - accepted
rej_rows = w.set_index('precursor_id').loc[list(rej)]
acc_rows = w.set_index('precursor_id').loc[list(acc_right)]
dec = w[(w.is_decoy == 1) & (w.s >= cut - 1.0)]
print(f'\nScoring loss: {len(rej)} right-peak winners rejected. Their score: median {rej_rows.s.median():.2f}, '
      f'within 1 of the cut {np.mean(rej_rows.s >= cut - 1):.2f}; decoy winners within 1 of the cut {len(dec)}')
print(f'  {"feature":30s} {"rejected":>9s} {"accepted":>9s} {"decoys@cut":>10s}')
for j, nm in enumerate(NAMES):
    print(f'  {nm:30s} {np.median(X[rej_rows.row.values, j]):9.3f} {np.median(X[acc_rows.row.values, j]):9.3f} {np.median(X[dec.row.values, j]):10.3f}')
# Rejected by peak strength: matched scans and support
print('  rejected by transition_support:', pd.Series(X[rej_rows.row.values, 0]).value_counts().sort_index().to_dict())
print('  accepted by transition_support:', pd.Series(X[acc_rows.row.values, 0]).value_counts().sort_index().to_dict())
# Wrong-peak losses: where is the winner?
wp = has_near - near_is_winner
wr = w.set_index('precursor_id').loc[list(wp)].merge(human.set_index('precursor_id')[['rt_s']], left_index=True, right_index=True)
print(f'\nWrong-peak: {len(wp)}; winner |RT - DIA-NN| p50 {np.median(np.abs(wr.rt_apex - wr.rt_s)):.1f} s; winner accepted {len(wp & accepted)}')
