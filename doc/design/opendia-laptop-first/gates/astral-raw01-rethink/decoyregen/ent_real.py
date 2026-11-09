"""Entrapment-not-false check (CALIB section 7, candidate 2). Are the entrapment winners in the score tail real signal?
For family O, MLP seeds 0/1/2, take entrapment winners above the q<=0.05 cut ('tail') and all entrapment winners
('bulk'). For each: is the same precursor (modified sequence + charge) a DIA-NN ID at Q<=0.01 in this run, and is our
apex within 10 s of DIA-NN's RT? Also: DIA-NN's global entrapment ID rate.
Reading, declared before running: candidate 2 is supported on a run if removing the tail entrapment winners that DIA-NN
also reports within 10 s brings the 5% S_e/S_d band to include 1 for >=2 of 3 seeds. (This removes them from the tail
count only, keeping the denominators, so it is a bound on how much of the excess such winners carry.)
Note: a DIA-NN entrapment ID at the same RT does not prove the peptide is present; DIA-NN's own FDP_c is about 0.007,
so both tools may pick the same interfering signal. The check measures concentration, not truth."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
import sqlite3
RUN = os.environ['OPENDIA_RUN']; RUNNO = RUN[-2:]
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
dn = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                     columns=['Run', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Protein.Ids'])
dn = dn[dn.Run.str.endswith('_' + RUNNO) & (dn['Q.Value'] <= 0.01)]
dn['key'] = dn['Modified.Sequence'] + '/' + dn['Precursor.Charge'].astype(str)
dent = dn['Protein.Ids'].str.split(';').map(lambda a: all('_p_target' in x for x in a))
print(f'run {RUN}: DIA-NN 1% IDs {len(dn)}, entrapment {int(dent.sum())} ({dent.mean():.4f})', flush=True)
DRT = dict(zip(dn.key, dn.RT * 60))
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
pp = pd.read_sql('select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID', c).drop_duplicates('pid')
KEY = dict(zip(pp.pid, pp.ms + '/' + pp.z.astype(str)))
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
S.m = meta[keep].reset_index(drop=True); apex = S.m.apex.values
def band(ke, ne, sd):
    return f'{ke / ne / sd:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / sd:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / sd:.2f}]'
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    hit = np.array([KEY.get(S.prec[i]) in DRT for i in ew])
    near = np.array([(KEY.get(S.prec[i]) in DRT) and abs(DRT[KEY[S.prec[i]]] - apex[i]) <= 10 for i in ew])
    tail = sc[ew] >= t; sd = (sc[dw] >= t).mean()
    print(f'seed{seed}: entrapment tail {int(tail.sum())} | DIA-NN ID: tail {hit[tail].mean():.3f} bulk {hit.mean():.4f} | same RT (10 s): tail {near[tail].mean():.3f} ({int(near[tail].sum())}) bulk {near.mean():.4f}'
          f' | 5% S_e/S_d {band(int(tail.sum()), len(ew), sd)} -> without same-RT DIA-NN entrapment {band(int((tail & ~near).sum()), len(ew), sd)}', flush=True)
