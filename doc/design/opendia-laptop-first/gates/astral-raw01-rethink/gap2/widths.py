"""Peak width / boundary facts on the same precursors (raw 01): DIA-NN-confirmed in our library AND accepted by OSW.
Ours: search group (rt_start..rt_end), the stage-2 boundary picker (boundary_log_width/fwhm, ms2 bins).
OSW: leftWidth..rightWidth. DIA-NN: RT.Start..RT.Stop and FWHM."""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
DN = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet'
lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet'); lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
r = pd.read_parquet(DN, columns=['Run', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop', 'FWHM'])
r = r[r.Run.str.endswith('_01') & (r['Q.Value'] <= 0.01)]
r['key'] = r['Stripped.Sequence'] + '/' + r['Precursor.Charge'].astype(str)
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')][['precursor_id', 'key', 'mkey' if 'mkey' in ref else 'modseq', 'charge']]
ref['mkey'] = ref.modseq + '/' + ref.charge.astype(str)
r = r.drop_duplicates('key').merge(ref[['precursor_id', 'key', 'mkey']], on='key')
o = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'leftWidth', 'rightWidth'])
o = o[o.run_name.str.endswith('_01')]; o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
m = r.merge(o, on='mkey', suffixes=('', '_osw'))
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'precursor_id', 'rt_start', 'rt_apex', 'rt_end', 'matched_scans']).rename(columns={'candidate_row': 'sample_id'})
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t', usecols=['sample_id', 'boundary_log_width', 'boundary_log_fwhm'])
s = s.merge(f, on='sample_id', how='left')
m['dn_rt'] = m.RT * 60
c = s.merge(m[['precursor_id', 'dn_rt']], on='precursor_id')
c = c[(c.rt_apex - c.dn_rt).abs() <= 10].sort_values('matched_scans', ascending=False).drop_duplicates('precursor_id')
m = m.merge(c.drop(columns='dn_rt'), on='precursor_id')
print('precursors confirmed by DIA-NN and OSW, with our candidate at the DIA-NN apex:', len(m))
q = lambda x: f'median {np.median(x):5.1f}  IQR {np.quantile(x, .25):5.1f}-{np.quantile(x, .75):5.1f}'
print('DIA-NN  base width  (RT.Stop-RT.Start) s:', q((m['RT.Stop'] - m['RT.Start']) * 60))
print('DIA-NN  FWHM                         s:', q(m.FWHM * 60))
print('OSW     base width  (right-left)     s:', q(m.rightWidth - m.leftWidth))
print('ours    search group (end-start)     s:', q(m.rt_end - m.rt_start))
print('ours    boundary picker base width   s:', q(np.exp(m.boundary_log_width)))
print('ours    boundary picker FWHM         s:', q(np.exp(m.boundary_log_fwhm)))
print('ours    matched scans in group        :', q(m.matched_scans))
# apex agreement
print('apex offset ours - DIA-NN s:', q(m.rt_apex - m.dn_rt), '| OSW - DIA-NN', q(m.RT_osw - m.dn_rt))
# overlap of our group with DIA-NN span
lo = np.maximum(m.rt_start, m['RT.Start'] * 60); hi = np.minimum(m.rt_end, m['RT.Stop'] * 60)
print('fraction of our group inside DIA-NN span:', q(np.clip(hi - lo, 0, None) / (m.rt_end - m.rt_start).clip(lower=1e-6)))
m.to_parquet('gap2/common.parquet')
