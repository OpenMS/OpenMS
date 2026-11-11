"""Where the DIA-NN 1 % precursors (human, exact library match; reference only) are lost for D11F / D11FW,
at the FDP_c 1 % and 2 % cuts, with the D-stage score gap and the 12-fragment library cosine of the right window."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
ref_all = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
for raw, tag in (('raw01', '_01'), ('raw02', '_02')):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); F = np.load(f'{O}/{raw}/feat12.npy'); H = np.load(f'{O}/{raw}/heuristics.npy')
    ref = ref_all[ref_all.Run.str.endswith(tag) & ~ref_all.entrap].drop_duplicates('precursor_id').copy(); ref['rt_s'] = ref.RT * 60
    for arm in ('D11F', 'D11FW'):
        s = np.load(f'{O}/{raw}/scores_{arm}.npy')
        w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True])
        win = w.drop_duplicates('precursor_id').copy()
        for f_, g in win[win.is_decoy == 1].groupby('fold'):
            med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
        t = win[win.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values).sort_values('s', ascending=False)
        f = np.cumsum(t.e.values) * (1 + 1 / R) / np.arange(1, len(t) + 1)
        cuts = {}
        for q in (0.01, 0.02):
            n = int(np.where(f <= q)[0].max() + 1); cuts[q] = (n, t.s.values[n - 1], set(t.precursor_id.values[:n]))
        cw = pop.assign(row=np.arange(len(pop)), s=s)[pop.is_decoy == 0].merge(ref[['precursor_id', 'rt_s']], on='precursor_id')
        cw['right'] = (cw.rt_apex - cw.rt_s).abs() <= 10
        has_right = set(cw[cw.right].precursor_id); winrow = dict(zip(win.precursor_id, win.row))
        right_rows = set(cw[cw.right].row)
        right_win = {p for p in has_right if winrow.get(p) in right_rows}
        inpop = set(pop.precursor_id)
        line = []
        for q, (n, cut, acc) in cuts.items():
            st = np.select([~ref.precursor_id.isin(inpop), ~ref.precursor_id.isin(has_right), ~ref.precursor_id.isin(right_win), ~ref.precursor_id.isin(acc)], ['A', 'B', 'C', 'D'], 'E')
            v = pd.Series(st).value_counts().reindex(list('ABCDE')).fillna(0).astype(int)
            line.append(f'@{int(q*100)}% N={n}: ' + ' '.join(f'{k}={v[k]}' for k in 'ABCDE'))
        print(raw, arm, ' | '.join(line))
        # D-stage detail at 1 %
        n, cut, acc = cuts[0.01]
        dset = ref[ref.precursor_id.isin(right_win) & ~ref.precursor_id.isin(acc)].precursor_id
        ns = win.set_index('precursor_id').s
        d = ns.loc[dset]; print(f'   D at 1 %: n {len(d)}, score/cut median {np.median(d / cut):.2f}; above 0.5x cut {int((d > 0.5 * cut).sum())}; '
              f'accepted at 2 % {int(dset.isin(cuts[0.02][2]).sum())}; transition_support median {np.median(H[winrow_arr := np.array([winrow[p] for p in dset]), 0]):.0f}, '
              f'f12_libcos median {np.median(F[winrow_arr, 2]):.2f} (accepted E median {np.median(F[[winrow[p] for p in acc if p in set(ref.precursor_id)], 2]):.2f})')
