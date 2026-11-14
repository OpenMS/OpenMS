"""Post-hoc asymmetry check for block I (evaluation only). The design doc's risk: decoys share signal less often than
targets, so a shared-signal penalty could hit false targets and decoys unequally and bias TDC. Check: among precursor
winners above the q <= 0.10 and q <= 0.05 cuts of an arm (paired TDC as diagnose.py), rho = E (1 + 1/r) / D split by
tertiles of the winner's i_explained and i_nclaim (tertile edges from all decoy winners, label-free), for the control
D11FDWT and D11FDIWT. Also: share of decoy vs target winners with any claimant contribution. rho near 1 in every
stratum means the decoys still sample false targets where the block acts. Usage: asym_featI.py <root>."""
import sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import beta as betad
O = sys.argv[1]
c = sqlite3.connect(f'file:{O}/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
tgt = pd.read_sql('select ID precursor_id from PRECURSOR where DECOY = 0', c).precursor_id
R = (ent.loc[tgt] == 1).sum() / (ent.loc[tgt] == 0).sum(); k = 1 + 1 / R


def rho(E, D):
    if D == 0: return f'E {E:4d} D {D:4d} rho   inf'
    lo = betad.ppf(0.05, E, D + 1) if E > 0 else 0.0; hi = betad.ppf(0.95, E + 1, D)
    f = lambda p: k * p / (1 - p)
    return f'E {E:4d} D {D:4d} rho {k * E / D:5.2f} [{f(lo):4.2f}, {f(hi):5.2f}]'


for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); FI = np.load(f'{O}/{raw}/featI.npy')
    for arm in ('D11FDWT', 'D11FDIWT'):
        s = np.load(f'{O}/{raw}/scores_{arm}.npy').astype(float)
        w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
        for f_, g in w[w.is_decoy == 1].groupby('fold'):
            med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
        best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
        T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
        best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
        w['ie'] = FI[w.row.values, 0]; w['nc'] = FI[w.row.values, 3]
        dw = w[w.is_decoy == 1]
        print(f'==== {raw} {arm}: winners with a claimant contribution: targets {np.mean(w[w.is_decoy == 0].ie > 0):.3f} decoys {np.mean(dw.ie > 0):.3f}')
        for col, name in (('ie', 'i_explained'), ('nc', 'i_nclaim')):
            edges = np.unique(np.quantile(dw[col], [1 / 3, 2 / 3]))
            for qc in (0.05, 0.10):
                cut = best[best.q <= qc].s.min()
                a = w[w.s >= cut].copy(); a['e'] = a.precursor_id.map(ent)
                a['t'] = np.searchsorted(edges, a[col].values, side='right')
                line = f'  {name:11s} q<={qc:.2f}'
                for t_ in range(len(edges) + 1):
                    g = a[a.t == t_]
                    line += f' | tertile {t_}: N {int((g.is_decoy == 0).sum()):4d} ' + rho(int(((g.is_decoy == 0) & (g.e == 1)).sum()), int((g.is_decoy == 1).sum()))
                print(line, flush=True)
