"""Voxel v2 diagnosis part 5: intensity and fragment concentration of top false winners against base hits (V0, FS)."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'; R = 148983 / 151239
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv', 'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
def winners(pop, s):
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    return w.sort_values('s', ascending=False)
def desc(xic, lib, rows):
    x = np.asarray(xic[np.sort(rows)]); l = np.asarray(lib[np.sort(rows)]); a = x.sum(2)  # [n,6] per-fragment area
    tot = a.sum(1); top = np.sort(a, 1)[:, ::-1]
    nfr = (a > 0.05 * tot[:, None]).sum(1)
    # correlation of per-fragment apex-centred traces with their median trace (coelution in the voxel)
    xc = x[:, :, 12:20]; z = xc - xc.mean(2, keepdims=True); z /= np.sqrt((z ** 2).sum(2, keepdims=True)) + 1e-12
    m = z.mean(1, keepdims=True); m /= np.sqrt((m ** 2).sum(2, keepdims=True)) + 1e-12; co = (z * m).sum(2)
    lc = (np.sqrt(a) * np.sqrt(l)).sum(1) / (np.sqrt(a.sum(1)) * np.sqrt(l.sum(1)) + 1e-12)
    return f'MS2 sum {np.median(tot):8.2g}  top1 share {np.median(top[:, 0] / tot):.2f}  frags>5% {np.median(nfr):.0f}  frags r>0.8 {np.median((co > 0.8).sum(1)):.0f}  lib cos {np.median(lc):.2f}'
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r'); lib = np.load(f'{O}/{raw}/lib.npy', mmap_mode='r')
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    wb = winners(pop, base); tb = wb[wb.is_decoy == 0].copy(); tb['e'] = tb.precursor_id.map(ent).values
    f = np.cumsum(tb.e.values) * (1 + 1 / R) / np.arange(1, len(tb) + 1); n1 = int(np.where(f <= 0.01)[0].max() + 1)
    print(f'\n== {raw}\n  base human hits@1%     ' + desc(xic, lib, tb.iloc[:n1][tb.e.iloc[:n1] == 0].row.values))
    for arm in ('V0', 'FS'):
        w = winners(pop, np.load(f'{O}/{raw}/scores_{arm}.npy')); t = w[w.is_decoy == 0].reset_index(drop=True); t['e'] = t.precursor_id.map(ent).values
        d = w[w.is_decoy == 1]; cut = t.s.iloc[1499]
        print(f'  {arm} entrap top-1500     ' + desc(xic, lib, t.iloc[:1500][t.e.iloc[:1500] == 1].row.values))
        print(f'  {arm} decoys above cut    ' + desc(xic, lib, d[d.s >= cut].row.values))
        print(f'  {arm} all decoy winners   ' + desc(xic, lib, d.sample(5000, random_state=1).row.values))
        print(f'  {arm} rate above cut: entrapment {(t.iloc[:1500].e == 1).sum() / (t.e == 1).sum() * 1e3:.2f}/1000, decoys {(d.s >= cut).sum() / len(d) * 1e3:.2f}/1000')
