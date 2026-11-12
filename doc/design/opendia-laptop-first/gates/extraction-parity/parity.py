"""Trace-by-trace extraction parity: OpenDIA xic native-scan traces vs DIA-NN 2.2.0 --xic export, raw01.

Pairs each DIA-NN fragment trace (pr, annotation) with our trace of the same library fragment (matched on product m/z)
and each DIA-NN point with our scan of the same RT. Reports coverage, zero/nonzero agreement, intensity ratio and
correlation, split by DIA-NN q-value group. Labels are not used.
"""
import sqlite3, sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
PY = sys.argv[1] if len(sys.argv) > 1 else 'TRACE_PARITY_raw01.txt'
out = open(X + '/' + PY, 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
tr = pd.read_sql("""select tp.PRECURSOR_ID precursor_id, t.PRODUCT_MZ product_mz, t.ANNOTATION annotation
  from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tp on tp.TRANSITION_ID=t.ID""", c)
tr = tr[tr.precursor_id.isin(idm.precursor_id)]

ours = pd.read_csv(X + '/ours/traces.tsv', sep='\t')
# slot -> library annotation: slots are the candidate's fragments in (calibrated) m/z order, the library order is the same.
tr = tr.sort_values(['precursor_id', 'product_mz']); tr['slot'] = tr.groupby('precursor_id').cumcount()
ours = ours.merge(tr[['precursor_id', 'slot', 'annotation', 'product_mz']], on=['precursor_id', 'slot'], how='left')
ppm_shift = ((ours.fragment_mz - ours.product_mz) / ours.product_mz * 1e6)
say('slot map: library annotation found for', f'{ours.annotation.notna().mean():.4f}', 'of our rows; calibrated-minus-library m/z ppm',
    ppm_shift.describe()[['min', '50%', 'max']].round(3).to_dict())
ours = ours.merge(idm[['precursor_id', 'Precursor.Id', 'Q.Value']], on='precursor_id')
ours['annotation'] = ours.annotation.str.replace(r'^\D*', lambda m: m.group(0), regex=True)

dn = pd.read_parquet(X + '/diann/report_raw01_xic/raw01.xic.parquet')
dn = dn[~dn.feature.isin(['index'])]
dn['rt_s'] = dn.rt.astype(float) * 60.0
dms1 = dn[dn.feature == 'ms1']; dn = dn[dn.feature != 'ms1']
rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop'])
grp = pd.cut(rep.set_index('Precursor.Id')['Q.Value'], [-1, 0.01, 0.05, 1], labels=['q<=1%', '1-5%', '>5%'])

# --- coverage of the window ---
o_rng = ours.groupby('Precursor.Id').rt.agg(['min', 'max'])
d_rng = dn.groupby('pr').rt_s.agg(['min', 'max'])
rr = d_rng.join(o_rng, lsuffix='_dn', rsuffix='_us', how='left')
apex = rep.set_index('Precursor.Id').RT * 60
cov = ((rr.min_us <= rr.min_dn + 2) & (rr.max_us >= rr.max_dn - 2))
inside = ((rr.min_us <= apex.reindex(rr.index) - 10) & (rr.max_us >= apex.reindex(rr.index) + 10))
say(f'\nwindow: DIA-NN precursors {len(rr)}; we wrote traces for {rr.min_us.notna().sum()}')
say(f'  our RT window contains DIA-NN apex +-10 s: {inside.mean():.4f}; contains DIA-NN full +-60 s export: {cov.mean():.4f}')
say('  by group, apex +-10 s inside our window:', inside.groupby(grp.reindex(rr.index).values, observed=True).mean().round(4).to_dict())

# --- fragment pairing ---
dfr = dn[['pr', 'feature']].drop_duplicates(); ofr = ours[['Precursor.Id', 'annotation']].drop_duplicates()
ofr.columns = ['pr', 'feature']
both = dfr.merge(ofr, on=['pr', 'feature'])
say(f'\nfragments: DIA-NN exported {len(dfr)}, ours (library) {len(ofr)}, paired {len(both)} '
    f'({len(both) / len(dfr):.4f} of DIA-NN traces)')

# --- point pairing on RT (same scans: the RT of DIA-NN's point equals our scan RT up to float rounding) ---
dn2 = dn.merge(both, on=['pr', 'feature'])
us2 = ours.rename(columns={'Precursor.Id': 'pr', 'annotation': 'feature'}).merge(both, on=['pr', 'feature'])
dn2 = dn2.sort_values('rt_s'); us2 = us2.sort_values('rt')
m = pd.merge_asof(dn2[['pr', 'feature', 'rt_s', 'value']], us2[['pr', 'feature', 'rt', 'intensity', 'observed_mz', 'fragment_mz']],
                  left_on='rt_s', right_on='rt', by=['pr', 'feature'], direction='nearest', tolerance=0.5)
m['dt'] = (m.rt - m.rt_s)
say(f'\npoints: DIA-NN points on paired traces {len(m)}; with an our-scan within 0.5 s {m.rt.notna().mean():.4f}')
say('  |RT difference| s quantiles', m.dt.abs().quantile([.5, .9, .99, 1]).round(4).to_dict())
mm = m[m.rt.notna()].copy()
mm['g'] = grp.reindex(mm.pr).values
dz = mm.value > 0; uz = mm.intensity > 0
say(f'  nonzero: DIA-NN {dz.mean():.4f}, ours {uz.mean():.4f}')
say(f'  DIA-NN>0 & ours=0 {(dz & ~uz).mean():.4f} | DIA-NN=0 & ours>0 {(~dz & uz).mean():.4f} | both>0 {(dz & uz).mean():.4f}')
b = mm[dz & uz]
ratio = b.intensity / b.value
say('  ours/DIA-NN intensity ratio where both > 0, quantiles', ratio.quantile([.01, .1, .5, .9, .99]).round(4).to_dict())
say(f'  exactly equal (rel diff < 1e-4): {(np.abs(ratio - 1) < 1e-4).mean():.4f}')
say('  by group: share DIA-NN>0 & ours=0, share ours>0 & DIA-NN=0, median ratio')
for g, s in mm.groupby('g', observed=True):
    d1 = s.value > 0; u1 = s.intensity > 0; bb = s[d1 & u1]
    say(f'    {g}: n={len(s)} miss={(d1 & ~u1).sum() / max(1, d1.sum()):.4f} extra={(~d1 & u1).sum() / max(1, u1.sum()):.4f} '
        f'ratio={float((bb.intensity / bb.value).median()):.4f}')

# --- where the DIA-NN points we miss sit (near the DIA-NN apex or in the tails) and their size ---
miss = mm[dz & ~uz].copy()
miss['from_apex'] = (miss.rt_s - apex.reindex(miss.pr).values).abs()
say('\nmissed DIA-NN points: |RT - DIA-NN apex| s quantiles', miss.from_apex.quantile([.1, .5, .9]).round(1).to_dict())
rel = miss.value / mm[dz].groupby(['pr', 'feature']).value.transform('max').reindex(miss.index)
say('  missed point / trace max quantiles', rel.quantile([.5, .9, .99]).round(4).to_dict())
near = mm[(mm.rt_s - apex.reindex(mm.pr).values).abs() <= 5]
dn_ = near.value > 0; un_ = near.intensity > 0
say(f'  within 5 s of DIA-NN apex: DIA-NN>0 & ours=0 share of DIA-NN>0 {(dn_ & ~un_).sum() / dn_.sum():.4f}; '
    f'ours>0 & DIA-NN=0 share of ours>0 {(~dn_ & un_).sum() / un_.sum():.4f}')

# --- per-trace correlation over the shared scans ---
def corr(s):
    a = s.value.values; u = s.intensity.values
    if a.std() == 0 or u.std() == 0: return np.nan
    return np.corrcoef(a, u)[0, 1]
cr = mm.groupby(['pr', 'feature']).apply(corr, include_groups=False)
say('\nper-trace Pearson r (ours vs DIA-NN, same scans):', cr.quantile([.05, .1, .25, .5]).round(4).to_dict(),
    f'| r>0.99 {float((cr > 0.99).mean()):.4f} | r<0.9 {float((cr < 0.9).mean()):.4f}')
mm.to_parquet(X + '/parity_points_raw01.parquet')
out.close()
