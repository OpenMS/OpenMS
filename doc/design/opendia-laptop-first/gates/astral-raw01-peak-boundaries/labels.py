"""Step A (diagnostic, raw 01): quality of the anchor peak-boundary labels from OpenDIACalibrate -out_anchor_peaks.
Compares picked widths with the stage-1 FWHM, our search candidate widths and same-library OpenSWATH boundaries
(the latter only as a reference for scale; never used as labels)."""
import numpy as np, pandas as pd, json, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
O = D + '/anchor-peaks'
OSW = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv'

p = pd.read_csv(O + '/anchor_peaks.tsv', sep='\t')
cal = json.load(open(O + '/cal.json'))
print('anchors', len(p), '| found', int(p.found.sum()), '| clean', int(p.clean.sum()), f'({p.clean.mean():.0%})')
print('boundary reasons left', p.left_reason.value_counts().to_dict(), 'right', p.right_reason.value_counts().to_dict())
print('agreeing fragments', p.agreeing_fragments.value_counts().sort_index().to_dict())
c = p[p.clean == 1]
q = [.1, .25, .5, .75, .9]
print('clean FWHM s', c.fwhm_s.quantile(q).round(2).to_dict(), '| stage-1 FWHM (cal.json)', round(cal['stage1']['peak_fwhm_s'], 2))
print('clean base width s', c.width_s.quantile(q).round(2).to_dict())
print('clean scans in peak', c.scans.quantile(q).to_dict())
print('apex shift vs calibration apex (s)', (c.apex_rt - c.calibration_rt).abs().quantile([.5, .9, .99]).round(2).to_dict())
print('asymmetry (right-apex)/(apex-left) median', round(((c.right_rt - c.apex_rt) / (c.apex_rt - c.left_rt).replace(0, np.nan)).median(), 2))

# Search candidates for the same precursors near the same apex.
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex', 'peak_width_seconds'])
c = c.assign(precursor_id=c.precursor_id.astype(int))
j = c.merge(s, on='precursor_id')
j = j[(j.rt_apex - j.apex_rt).abs() <= 3].sort_values('peak_width_seconds').drop_duplicates('precursor_id')
print(f'\nsame peak in search output: {len(j)} | search width median {j.peak_width_seconds.median():.1f} s vs picked base width {j.width_s.median():.1f} s')

# OpenSWATH boundaries for the same precursors (scale reference only).
lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet')
lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
o = pd.read_csv(OSW, sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'leftWidth', 'rightWidth'])
o = o[o.run_name.str.endswith('_01')]
o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
o = o.merge(lib[['mkey', 'precursor_id']], on='mkey')
k = c.merge(o, on='precursor_id')
k = k[(k.RT - k.apex_rt).abs() <= 5]
if len(k):
    print(f'same peak in OSW output: {len(k)} | OSW width median {(k.rightWidth - k.leftWidth).median():.1f} s vs picked base width {k.width_s.median():.1f} s, FWHM {k.fwhm_s.median():.1f} s')
    print('  |left diff| s', (k.leftWidth - k.left_rt).abs().quantile([.5, .9]).round(1).to_dict(), '|right diff| s', (k.rightWidth - k.right_rt).abs().quantile([.5, .9]).round(1).to_dict())
