"""D11 misses, part d (descriptive): entrapment target winners near the FDP_c 1 % cut, and the score shape near the cut."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); H = np.load(f'{O}/{raw}/heuristics.npy')
    s = np.load(f'{O}/{raw}/scores_D11.npy')
    st = pd.read_csv(f'{O}/v3seq/misses_D11_{raw}.tsv', sep='\t').set_index('precursor_id').stage
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    t = w[w.is_decoy == 0].assign(e=lambda x: x.precursor_id.map(ent).values).sort_values('s', ascending=False).reset_index(drop=True)
    t['N'] = np.arange(1, len(t) + 1); t['n_e'] = t['e'].cumsum(); t['fdp'] = t['n_e'] * (1 + 1 / R) / t['N']
    t['stage'] = t.precursor_id.map(st).fillna('not in DIA-NN ref')
    print(f'\n=== {raw}: first 15 entrapment target winners (N = rank, FDP_c at that rank)')
    for _, r in t[t.e == 1].head(15).iterrows():
        h = H[int(r.row)]
        print(f'   N {int(r.N):5d} fdp {r.fdp:.4f} fold {int(r.fold)} s {r.s:8.1f} | ts {h[0]:.0f} rt_err {h[3]:5.1f} coel {h[6]:.2f} libcorr {h[7]:.2f} apex {h[8]:.2f} b {h[9]:.0f} y {h[10]:.0f}')
    n1 = int(np.where(t.fdp <= 0.01)[0].max() + 1)
    seg = t.iloc[n1:n1 + 1000]
    print(f'  ranks {n1 + 1}..{n1 + 1000} after the cut: entrapment {int(seg.e.sum())}, DIA-NN-confirmed {int((seg.stage == "D right window, below cut").sum())}, '
          f'not in DIA-NN ref {int((seg.stage == "not in DIA-NN ref").sum())}')
    print('  N at FDP_c 1/1.5/2/3/5 %:', [int(np.where(t.fdp <= x)[0].max() + 1) for x in (0.01, 0.015, 0.02, 0.03, 0.05)])
    print('  DIA-NN-confirmed (stage D/E) among the top N at those cuts:', [int(t.iloc[:int(np.where(t.fdp <= x)[0].max() + 1)].stage.isin(['E accepted', 'D right window, below cut']).sum()) for x in (0.01, 0.015, 0.02, 0.03, 0.05)])
