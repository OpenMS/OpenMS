"""In-data test of the decoy mass-shift hypothesis (CALIB.md reading 2). Design diagnostic using cached scores.

Each decoy has a precursor mass shift against its target: +-12, +-28 or +-44 Da, from the pseudo-reverse residue swap.
If the shift is what makes decoys under-represent false targets in the tail, decoy tail survival should depend on the
shift: for example, lower for decoys moved to sparser m/z regions (positive shifts), or for larger shifts. Tail survival is
the share of decoy winners above the q<=0.02 and q<=0.05 cuts, shown against the entrapment winners' share. If survival is flat
across shift groups, the shift hypothesis is weakened.
"""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
d = pd.read_sql("""select p.ID pid, p.PRECURSOR_MZ mz, p.CHARGE z, p.DECOY dec, group_concat(pr.PROTEIN_ACCESSION,';') acc
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
d['key'] = d.acc.str.replace('DECOY_', '')
j = d[d.dec == 1].merge(d[d.dec == 0][['key', 'z', 'mz']], on=['key', 'z'], suffixes=('', '_t'))
j['dM'] = ((j.mz - j.mz_t) * j.z).round(0)
shift = j.drop_duplicates('pid').set_index('pid').dM
z = np.load(f'{S.S}/calib_scores.npz')
print(f'raw {RUN}')
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    sh = shift.reindex(S.prec[dw]).values
    for qq in (0.02, 0.05):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        se = (sc[ew] >= t).mean()
        parts = []
        for lab, msk in (('neg', sh < 0), ('pos', sh > 0), ('|12|', np.abs(sh) == 12), ('|28|', np.abs(sh) == 28), ('|44|', np.abs(sh) == 44)):
            k = int((sc[dw[msk]] >= t).sum()); n = int(msk.sum())
            parts.append(f'{lab} {k}/{n} ratio {se / max(1e-12, k / max(1, n)):.2f}')
        print(f'  {name.upper()} cut q<={qq}: S_e {se * 1e4:.2f}e-4 | S_e/S_d by decoy shift: ' + ' | '.join(parts))
