"""Our accepted set vs OSW's on raw 01, located on the OSW-extracted traces: for OSW-accepted precursors, split by
whether we accept them, and compare trace-level quantities at the OSW peak (from the XIC, baseline-free):
summed fragment apex, number of fragments with signal, median pairwise correlation inside OSW boundaries, and the
fraction of total fragment signal in +-97 s that lies inside OSW boundaries (an interference/background index)."""
import os, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic01'
frag = np.load(S + '/frag.npy', mmap_mode='r'); rt = np.load(S + '/rt.npy'); npts = np.load(S + '/npts.npy'); meta = pd.read_parquet(S + '/meta.parquet')
row = pd.Series(np.arange(len(meta)), index=meta.precursor_id)
lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet'); lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
o = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'leftWidth', 'rightWidth'])
o = o[o.run_name.str.endswith('_01')]; o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
o = o.merge(lib[['mkey', 'precursor_id']], on='mkey').drop_duplicates('precursor_id')
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'precursor_best', 'q_value'])
acc = set(cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (pd.to_numeric(cr.q_value, errors='coerce') <= 0.01)].precursor_id)
o['ours'] = o.precursor_id.isin(acc)
rows = []
for r in o.itertuples():
    i = row.get(r.precursor_id)
    if i is None: continue
    n = npts[i]; t = rt[i, :n]; X = np.asarray(frag[i, :, :n], np.float64)
    inside = (t >= r.leftWidth) & (t <= r.rightWidth)
    if inside.sum() < 3: continue
    seg = X[:, inside]; s = seg.sum(0)
    cc = np.corrcoef(seg + 1e-9 * np.random.rand(*seg.shape))
    rows.append(dict(ours=r.ours, apex=np.log10(s.max() + 1), nsig=int((seg.max(1) > 0).sum()),
                     corr=np.nanmedian(cc[np.triu_indices(6, 1)]), frac_in=seg.sum() / max(X.sum(), 1),
                     width=r.rightWidth - r.leftWidth, frag_min_share=np.sort(seg.sum(1) / max(seg.sum(), 1))[0]))
d = pd.DataFrame(rows)
print('OSW-accepted precursors with XIC', len(d), '| also accepted by us', int(d.ours.sum()))
print(d.groupby('ours').median().round(3).T)
d['apex_bin'] = pd.qcut(d.apex, 5, labels=False)
print('\nour recall of OSW IDs by OSW-peak apex intensity quintile:', d.groupby('apex_bin').ours.mean().round(2).to_dict())
d['corr_bin'] = pd.qcut(d['corr'], 5, labels=False)
print('our recall by fragment-correlation quintile (inside OSW bounds):', d.groupby('corr_bin').ours.mean().round(2).to_dict())
d['in_bin'] = pd.qcut(d.frac_in, 5, labels=False)
print('our recall by in-peak signal fraction quintile (low = more background/interference):', d.groupby('in_bin').ours.mean().round(2).to_dict())
