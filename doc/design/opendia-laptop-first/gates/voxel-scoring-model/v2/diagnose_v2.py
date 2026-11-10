"""Voxel v2 diagnosis (after the frozen result, no new training, does not change any verdict).
Question: why do the voxel arms lose IDs against the base? Splits the loss into
  (a) window choice: which candidate window each precursor's winner is, and
  (b) ranking: how the arm orders winners against decoy winners.
Swap test: rank the base's winner windows by the arm score (arm ranking, base window choice), and the arm's winner
windows by the base score (base ranking, arm window choice). Also: score saturation, the decoys that the arm puts on
top, and the arm score of base hits that are calibration anchors versus non-anchors.
Entrapment labels are read here only for the evaluation numbers, as in evaluate_v2.py."""
import sqlite3, sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv',
         'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
ANCH = {'raw01': B + '/astral-design-raw01-20261107/window300k/anchors.tsv', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/anchors.tsv'}
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e

def winners(pop, s):
    """Per-precursor winner row (max score), then per-fold decoy-winner median/MAD standardisation (as evaluate_v2)."""
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    return std(w)

def std(w):
    w = w.copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    return w

def accepted(w, t=0.01):
    """Target winners accepted at entrapment FDP_c <= t (the evaluate_v2 N@t set) and the score cut."""
    tw = w[w.is_decoy == 0].sort_values('s', ascending=False)
    e = tw.precursor_id.map(ent).values
    f = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(f <= t)[0]
    n = int(ok.max() + 1) if len(ok) else 0
    return tw.iloc[:n], (tw.s.values[n - 1] if n else np.inf)

def regather(pop, s, win_rows):
    """Score fixed winner windows (one row per precursor) with score s, standardised the same way."""
    w = pop.iloc[win_rows].assign(s=s[win_rows], row=win_rows)
    return std(w)

for raw in sys.argv[1:] or ['raw01', 'raw02']:
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    H = pd.DataFrame(np.load(f'{O}/{raw}/heuristics.npy'), columns=FEATS)
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    anchors = set(pd.read_csv(ANCH[raw], sep='\t').query('used == 1').precursor_id)
    wb = winners(pop, base); ab, cut_b = accepted(wb)
    print(f'\n===== {raw}: base N@1% {len(ab)}  (anchors among them {ab.precursor_id.isin(anchors).sum()})')
    for arm in ('V0', 'F0', 'S0', 'FS'):
        s = np.load(f'{O}/{raw}/scores_{arm}.npy')
        wa = winners(pop, s); aa, cut_a = accepted(wa)
        lost = ab[~ab.precursor_id.isin(aa.precursor_id)]; gain = aa[~aa.precursor_id.isin(ab.precursor_id)]
        # (a) window choice: same winner window as base?
        wa_i = wa.set_index('precursor_id'); wb_i = wb.set_index('precursor_id')
        same_lost = (wa_i.loc[lost.precursor_id, 'row'].values == wb_i.loc[lost.precursor_id, 'row'].values).mean()
        dt_lost = np.abs(wa_i.loc[lost.precursor_id, 'rt_apex'].values - wb_i.loc[lost.precursor_id, 'rt_apex'].values)
        # swap tests
        n_arm_on_base_win = len(accepted(regather(pop, s, wb.row.values))[0])
        n_base_on_arm_win = len(accepted(regather(pop, base, wa.row.values))[0])
        # saturation: share of all logits in the top 1 % of the logit range; ties at the top of target winners
        raw_top = np.sort(s)[::-1]
        sat = (s >= raw_top[0] - 0.01 * (raw_top[0] - raw_top[-1])).mean()
        # decoys above the arm cut, their heuristics against base-accepted targets
        dec_hi = wa[(wa.is_decoy == 1) & (wa.s >= cut_a)]
        dec_hi_b = wb[(wb.is_decoy == 1) & (wb.s >= cut_b)]
        # arm score (standardised, arm winner) of base hits: anchors vs non-anchors, as the percentile among arm decoy winners
        dws = np.sort(wa[wa.is_decoy == 1].s.values)
        pct = lambda pids: 100 * np.searchsorted(dws, wa_i.loc[pids, 's'].values) / len(dws)
        an = ab.precursor_id[ab.precursor_id.isin(anchors)]; na = ab.precursor_id[~ab.precursor_id.isin(anchors)]
        print(f'--- {arm}: N@1% {len(aa)}  lost {len(lost)}  gained {len(gain)}  (gained entrapment {int(gain.precursor_id.map(ent).sum())})')
        print(f'    window choice: lost base hits with the SAME winner window {same_lost:.2f}; median |dRT| of moved winners '
              f'{np.median(dt_lost[dt_lost > 0]) if (dt_lost > 0).any() else 0:.1f} s')
        print(f'    swap: arm score on BASE windows N@1% {n_arm_on_base_win};  base score on ARM windows N@1% {n_base_on_arm_win}')
        print(f'    raw logit range {raw_top[-1]:.1f}..{raw_top[0]:.1f}; share of windows in top 1 % of range {sat:.4f}; '
              f'target winners tied at max {(np.isclose(wa[wa.is_decoy==0].s, wa.s.max())).sum()}')
        print(f'    base hits, arm-score percentile among arm decoy winners: anchors median {np.median(pct(an)):.1f} (n {len(an)}), '
              f'non-anchors median {np.median(pct(na)):.1f}, non-anchors below the 99th {(pct(na) < 99).mean():.2f}')
        hb = H.iloc[dec_hi.row.values].median(); ht = H.iloc[ab.row.values].median(); hdb = H.iloc[dec_hi_b.row.values].median() if len(dec_hi_b) else ht * np.nan
        hl = H.iloc[wa_i.loc[lost.precursor_id, 'row'].values].median()
        print(f'    medians        {"":10s} ' + ' '.join(f'{f[:12]:>12s}' for f in FEATS))
        for name, h in (('base hits', ht), ('arm windows of lost', hl), (f'arm top decoys ({len(dec_hi)})', hb), (f'base top decoys ({len(dec_hi_b)})', hdb)):
            print(f'    {name:26s} ' + ' '.join(f'{v:12.3f}' for v in h.values))
