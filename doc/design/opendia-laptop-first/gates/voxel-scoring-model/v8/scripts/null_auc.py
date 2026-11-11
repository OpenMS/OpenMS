"""Null-window symmetry: AUC(target > decoy) of a score over windows far from the library RT (|RT error| > 3 sigma).
No peptide signal is expected there, so a sequence or library tell shows up as AUC != 0.5. Window level and per-precursor
best null window. Usage: null_auc.py <arm> [...] (BASE = stage2_score). No entrapment labels."""
import sys, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/rtlin-20261111'  # rtlin rebuild root
BASEF = {r: O + '/' + r + '/byion/features.tsv' for r in ('raw01', 'raw02')}
CAL = {r: O + '/' + r + '/cal/cal.json' for r in ('raw01', 'raw02')}
def auc(p, n): return mannwhitneyu(p, n).statistic / (len(p) * len(n))
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); H = np.load(f'{O}/{raw}/heuristics.npy')
    sig = json.load(open(CAL[raw]))['rt']['sigma_s']; null = H[:, 3] > 3 * sig
    for arm in sys.argv[1:]:
        s = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values if arm == 'BASE' else np.load(f'{O}/{raw}/scores_{arm}.npy')
        d = pop[null].assign(s=s[null]); t, c = d[d.is_decoy == 0].s.values, d[d.is_decoy == 1].s.values
        w = d.sort_values('s', ascending=False).drop_duplicates('precursor_id'); tw, cw = w[w.is_decoy == 0].s.values, w[w.is_decoy == 1].s.values
        print(f'{raw} {arm:5s} null windows T {len(t)} D {len(c)}  AUC {auc(t, c):.4f} | null winners T {len(tw)} D {len(cw)}  AUC {auc(tw, cw):.4f}', flush=True)
