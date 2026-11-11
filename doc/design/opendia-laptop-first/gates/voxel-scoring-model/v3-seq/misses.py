"""D11 miss analysis (descriptive, no training, no change to any pipeline). Reference: DIA-NN 2.3.1 1 % precursors that
exactly match a library precursor (astral-diann-gap-20261108/ref_exact.parquet). DIA-NN is a reference only.
Funnel per raw: library precursor -> in the scoring population (pair max transition_support >= 3) -> a candidate window
within 10 s of the DIA-NN RT -> the D11 winner is that window -> the D11 winner is accepted at FDP_c 1 % (same ranking as
evaluate_v2.py: per-fold decoy median/MAD, precursor winners). Then the features of each loss group against accepted
targets and the top decoy winners. Entrapment labels are evaluation only."""
import sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239; TOL = 10.0; ARM = sys.argv[1] if len(sys.argv) > 1 else 'D11'
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
SHOW = ['transition_support', 'abs_library_rt_error', 'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement',
        'byion_b', 'byion_y', 'log_int', 'n_frag']
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pre = pd.read_sql("select p.ID precursor_id, p.CHARGE charge, p.DECOY decoy from PRECURSOR p", c).set_index('precursor_id')
ref_all = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
out = []
def P(*a): s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

for raw, tag in (('raw01', '_01'), ('raw02', '_02')):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    H = pd.DataFrame(np.load(f'{O}/{raw}/heuristics.npy'), columns=FEATS)
    xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r')
    s = np.load(f'{O}/{raw}/scores_{ARM}.npy')
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True])
    w['wrank'] = w.groupby('precursor_id').cumcount()
    # Per-fold decoy median/MAD normalisation of precursor winners (evaluate_v2.py).
    win = w[w.wrank == 0].copy()
    for f_, g in win[win.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
    t = win[win.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values).sort_values('s', ascending=False)
    f = np.cumsum(t.e.values) * (1 + 1 / R) / np.arange(1, len(t) + 1); ok = np.where(f <= 0.01)[0]; n1 = int(ok.max() + 1) if len(ok) else 0
    cut = t.s.values[n1 - 1]; acc = set(t.precursor_id.values[:n1])
    dwin = win[win.is_decoy == 1].sort_values('s', ascending=False)
    # Reference: DIA-NN 1 %, exact library match, this run, human (entrapment excluded so 'missed' means a likely-true ID).
    ref = ref_all[ref_all.Run.str.endswith(tag) & ~ref_all.entrap].drop_duplicates('precursor_id').copy(); ref['rt_s'] = ref.RT * 60
    ref['e'] = ref.precursor_id.map(ent)
    P(f'\n=== {raw} {ARM}: accepted at FDP_c 1 % = {n1} (score cut {cut:.2f}); DIA-NN reference (human, exact match) = {len(ref)}')
    inpop = ref[ref.precursor_id.isin(set(pop.precursor_id))]
    cw = w[w.is_decoy == 0].merge(ref[['precursor_id', 'rt_s']], on='precursor_id'); cw['right'] = (cw.rt_apex - cw.rt_s).abs() <= TOL
    has_right = set(cw[cw.right].precursor_id)
    rw = cw[cw.right].sort_values('wrank').drop_duplicates('precursor_id')  # best-ranked window at the DIA-NN RT
    wincw = cw[cw.wrank == 0].set_index('precursor_id')
    right_win = set(wincw[wincw.right].index)
    ref['stage'] = np.select(
        [~ref.precursor_id.isin(set(pop.precursor_id)), ~ref.precursor_id.isin(has_right), ~ref.precursor_id.isin(right_win),
         ~ref.precursor_id.isin(acc)],
        ['A not in population', 'B no window within 10 s', 'C wrong window won', 'D right window, below cut'], 'E accepted')
    P(ref.stage.value_counts().sort_index().to_string())
    acc_ref = ref[ref.stage == 'E accepted']; P(f'accepted and in DIA-NN: {len(acc_ref)} of {n1} accepted ({len(acc_ref) / max(n1, 1):.0%})')
    # Accepted by D11 but not in the DIA-NN reference.
    P(f'accepted, not DIA-NN-confirmed: {n1 - len(acc_ref)} (entrapment among them: {int(t.e.values[:n1].sum())} of all accepted)')

    # A: why not in the population. The population keeps pairs with max transition_support >= 3.
    # D: distance of the right window's score to the cut, and what drives it.
    def feats(rows):
        d = H.iloc[rows].copy(); x = np.asarray(xic[np.sort(rows)])[np.argsort(np.argsort(rows))]
        d['log_int'] = np.log10(x.sum((1, 2)) + 1); d['n_frag'] = (x.max(2) > 0).sum(1); return d
    norm = dict(zip(win.precursor_id, win.s))
    rw = rw.set_index('precursor_id')
    groups = {}
    groups['E accepted (right window)'] = rw.loc[rw.index.intersection(acc_ref.precursor_id), 'row'].values
    dD = ref[ref.stage == 'D right window, below cut'].precursor_id
    groups['D right window, below cut'] = rw.loc[dD, 'row'].values
    cC = ref[ref.stage == 'C wrong window won'].precursor_id
    groups['C right window (lost to a wrong one)'] = rw.loc[cC, 'row'].values
    groups['C wrong winner'] = w[(w.wrank == 0) & w.precursor_id.isin(cC)].row.values
    groups['top 1 % decoy winners'] = dwin.row.values[:max(len(dwin) // 100, 1)]
    groups['decoy winners above cut'] = dwin[dwin.s >= cut].row.values
    P('median features (window rows):')
    P(f'{"group":38s} {"n":>5s} ' + ' '.join(f'{k[:11]:>11s}' for k in SHOW))
    for name, rows in groups.items():
        if len(rows) == 0: continue
        m = feats(np.asarray(rows)).median()
        P(f'{name:38s} {len(rows):5d} ' + ' '.join(f'{m[k]:11.3f}' for k in SHOW))
    # D: how far below the cut (normalised winner score; the right window is the winner for D).
    dsc = pd.Series([norm[p] for p in dD], dtype=float)
    P(f'D normalised score quantiles (cut {cut:.2f}): ' + ', '.join(f'q{int(q * 100)} {dsc.quantile(q):.2f}' for q in (0.1, 0.25, 0.5, 0.75, 0.9)))
    P(f'D within 2x of the cut: {int((dsc >= cut / 2).sum())}, within 0.8x: {int((dsc >= 0.8 * cut).sum())}')
    # By charge, and DIA-NN q (an easy ID for DIA-NN vs a borderline one).
    ref['charge'] = ref.precursor_id.map(pre.charge)
    P('by charge (stage counts):'); P(pd.crosstab(ref.charge, ref.stage).to_string())
    ref['diann_q'] = pd.cut(ref['Q.Value'] if 'Q.Value' in ref else pd.Series(np.nan, index=ref.index), [0, 1e-3, 3e-3, 1e-2], include_lowest=True)
    # C: how far the wrong winner is from the DIA-NN RT.
    cwin = wincw.loc[wincw.index.intersection(cC)]
    P(f'C wrong winner |RT - DIA-NN RT| s: median {(cwin.rt_apex - cwin.rt_s).abs().median():.1f}, q90 {(cwin.rt_apex - cwin.rt_s).abs().quantile(.9):.1f}')
    ref[['precursor_id', 'key', 'charge', 'rt_s', 'stage']].to_csv(f'{O}/v3seq/misses_{ARM}_{raw}.tsv', sep='\t', index=False)
open(f'{O}/v3seq/MISSES_{ARM}.txt', 'w').write('\n'.join(out) + '\n')
