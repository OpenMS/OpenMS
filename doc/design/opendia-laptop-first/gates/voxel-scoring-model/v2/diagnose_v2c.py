"""Voxel v2 diagnosis part 3: what drives the voxel logit among DECOY winners (no labels beyond is_decoy).
Spearman of the arm logit with each heuristic and with the window's log total fragment intensity."""
import numpy as np, pandas as pd, sys
from scipy.stats import spearmanr
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); H = np.load(f'{O}/{raw}/heuristics.npy'); xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r')
    for arm in ('V0', 'FS'):
        s = np.load(f'{O}/{raw}/scores_{arm}.npy')
        w = pop.assign(s=s, row=np.arange(len(pop))).sort_values('s', ascending=False).drop_duplicates('precursor_id')
        d = w[w.is_decoy == 1].row.values
        inten = np.log1p(np.asarray(xic[np.sort(d)]).sum((1, 2)))[np.argsort(np.argsort(d))]
        r = {f: spearmanr(s[d], H[d, i])[0] for i, f in enumerate(FEATS)}; r['log_intensity'] = spearmanr(s[d], inten)[0]
        print(raw, arm, 'decoy winners n', len(d), ' '.join(f'{k}={v:+.2f}' for k, v in r.items()))
