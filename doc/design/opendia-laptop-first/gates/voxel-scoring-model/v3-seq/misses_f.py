"""D11 misses, part f (descriptive): upstream groups. A (not in population): best transition_support of the pair in the full
search table. B (no window within 10 s): distance from the DIA-NN RT to the nearest candidate, and its RT against the library."""
import json, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
P = {'raw01': B + '/astral-design-raw01-20261107/pipeline', 'raw02': B + '/astral-entrapment-raw02-20261108'}
CAL = {'raw01': B + '/astral-design-raw01-20261107/window300k', 'raw02': B + '/astral-entrapment-raw02-20261108/cal'}
for raw in ('raw01', 'raw02'):
    st = pd.read_csv(f'{O}/v3seq/misses_D11_{raw}.tsv', sep='\t')
    cand = pd.read_csv(P[raw] + '/search.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex', 'transition_support', 'apex_rt_error_seconds'])
    sig = float(json.load(open(CAL[raw] + '/cal.json'))['rt']['sigma_s'])
    m = cand.merge(st, on='precursor_id'); m['d'] = (m.rt_apex - m.rt_s).abs()
    A = m[m.stage == 'A not in population']
    print(f'\n{raw} A ({(st.stage == "A not in population").sum()}): with any search candidate {A.precursor_id.nunique()}; '
          f'best transition_support among them: {A.groupby("precursor_id").transition_support.max().value_counts().sort_index().to_dict()}; '
          f'transition_support of the candidate nearest the DIA-NN RT: {A.sort_values("d").drop_duplicates("precursor_id").transition_support.value_counts().sort_index().to_dict()}')
    Bm = m[m.stage == 'B no window within 10 s']; nb = Bm.sort_values('d').drop_duplicates('precursor_id')
    print(f'{raw} B ({(st.stage == "B no window within 10 s").sum()}): nearest candidate to the DIA-NN RT, distance s quantiles 25/50/75/90: '
          + ' '.join(f'{v:.0f}' for v in nb.d.quantile([.25, .5, .75, .9])))
    # Where is the DIA-NN RT relative to the calibrated library RT? apex_rt_error_seconds = apex - predicted; predicted = apex - err.
    nb = nb.assign(pred=nb.rt_apex - nb.apex_rt_error_seconds); nb['dia_err'] = (nb.rt_s - nb.pred).abs()
    print(f'{raw} B: |DIA-NN RT - calibrated library RT| quantiles 25/50/75/90 s: ' + ' '.join(f'{v:.0f}' for v in nb.dia_err.quantile([.25, .5, .75, .9]))
          + f' (sigma {sig:.1f}); beyond 3 sigma: {int((nb.dia_err > 3 * sig).sum())}')
    E = m[m.stage == 'E accepted'].sort_values('d').drop_duplicates('precursor_id'); E = E.assign(pred=E.rt_apex - E.apex_rt_error_seconds)
    print(f'{raw} E (accepted): |DIA-NN RT - calibrated library RT| median {(E.rt_s - E.pred).abs().median():.0f} s')
