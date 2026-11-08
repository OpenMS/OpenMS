"""Are tail entrapment hits I/L-equivalent to a human target, or identical to a human target under I=L plus a K/Q or N/D
swap? Also: how many were reported by DIA-NN on the same run (a reference only)? Design diagnostic using cached scores."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
d = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq, group_concat(pr.PROTEIN_ACCESSION,';') acc
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c).set_index('pid')
d['ent'] = d.acc.str.split(';').apply(lambda a: all('_p_target' in x for x in a))
il = lambda s: s.replace('I', 'L'); il2 = lambda s: il(s).replace('Q', 'K').replace('N', 'D')
hum = d[~d.ent]; H1 = set(hum.seq.map(il)); H2 = set(hum.seq.map(il2))
e = d[d.ent]
print(f'raw {RUN}: library entrapment precursors {len(e)}; I/L-equal to human {e.seq.map(il).isin(H1).mean():.4f}; +Q/K,N/D {e.seq.map(il2).isin(H2).mean():.4f}')
dn = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet', columns=['Run', 'Modified.Sequence', 'Precursor.Charge', 'Q.Value'])
dn = dn[dn.Run.str.endswith('_' + RUN) & (dn['Q.Value'] <= 0.01)]
dnk = set(zip(dn['Modified.Sequence'].str.replace('(UniMod:4)', '', regex=False).str.replace('(UniMod:35)', '(UniMod:35)', regex=False), dn['Precursor.Charge']))
z = np.load(f'{S.S}/calib_scores.npz')
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    for qq in (0.01, 0.05):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        hp = S.prec[w[(~S.dec[w]) & S.entr[w] & (sc[w] >= t)]]
        x = d.loc[[p for p in hp if p in d.index]]
        indn = sum((s.replace('(UniMod:4)', ''), zz) in dnk for s, zz in zip(x.mseq, x.z))
        print(f'  {name.upper()} q<={qq}: entrapment hits {len(x)}  I/L-equal {x.seq.map(il).isin(H1).sum()}  +QK/ND {x.seq.map(il2).isin(H2).sum()}  in DIA-NN 1% {indn}')
