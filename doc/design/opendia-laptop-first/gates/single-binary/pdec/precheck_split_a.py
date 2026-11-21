"""PRECHECK_SPLIT part A (diagnostic, 21 Nov 2026): does a peak-splitting arm cut real peaks?
Compares the best window per precursor (precursors.tsv) of the def run with an arm run, against DIA-NN's peak extent.
Usage: precheck_split_a.py <arm name> <arm run root>   (arm runs are <root>/<R>_<L>_<arm>/; def runs from ceil0-20261121)"""
import sys, glob, sqlite3, numpy as np, pandas as pd
from scipy.stats import spearmanr
ARM, ROOT = sys.argv[1:3]
C = '/nfs/research/juan/singj/opendia'; W0 = f'{C}/ceil0-20261121'
LIB = {'r0': f'{C}/libsize-20261121/r0/lib.pqp', 'r7': f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'}
DN = {('D', 'r0'): f'{C}/libsize-20261121/r0/diann/D/report.parquet', ('H', 'r0'): f'{C}/libsize-20261121/r0/diann/H/report.parquet',
      ('D', 'r7'): f'{C}/entlib-diann5-20261119/D/report.parquet', ('H', 'r7'): f'{C}/entlib-diann5-20261119/H/report.parquet'}
labs = {}
def labels(L):
    if L not in labs:
        c = sqlite3.connect(f'file:{LIB[L]}?mode=ro', uri=True)
        labs[L] = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
          min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
          from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
          join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
    return labs[L]
def best(rd):
    f = [x for x in glob.glob(f'{rd}/*.precursors.tsv') if not x.endswith('/precursors.tsv')][0]
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'rt_start', 'rt_end', 'area', 'score', 'q_value'])
    t = t[(t.is_decoy == 0) & (t.score > -1e29)].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    return t.set_index('precursor_id')
out = []
for R in ('D', 'H'):
    for L in ('r0', 'r7'):
        tg = labels(L); ent = tg.set_index('pid').e
        d = pd.read_parquet(DN[(R, L)], columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop', 'Precursor.Quantity'])
        d = d[d['Q.Value'] <= 0.01]; t2 = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
        k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
        d['pid'] = t2.pid.reindex(k).values; d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
        d = d[ent.reindex(d.pid).values == 0].sort_values('Q.Value').drop_duplicates('pid').set_index('pid')
        d['rt'] = d.RT * 60; d['s0'] = d['RT.Start'] * 60; d['s1'] = d['RT.Stop'] * 60
        A = best(f'{W0}/{R}_{L}'); B = best(f'{ROOT}/{R}_{L}_{ARM}')
        accA = set(A.index[(A.q_value <= 0.01) & (A.index.map(ent).fillna(1) == 0)])
        accB = set(B.index[(B.q_value <= 0.01) & (B.index.map(ent).fillna(1) == 0)])
        dn = set(d.index)
        def cov(T, ids):
            x = T.reindex(ids); y = d.reindex(ids)
            ov = (np.minimum(x.rt_end, y.s1) - np.maximum(x.rt_start, y.s0)).clip(lower=0)
            return ov / (y.s1 - y.s0).clip(lower=1e-6)
        def split_flag(ids):
            a = A.reindex(ids); b = B.reindex(ids)
            wa = a.rt_end - a.rt_start; wb = b.rt_end - b.rt_start
            return ((a.rt_apex - b.rt_apex).abs() <= 3) & (wb < 0.7 * wa) & (b.rt_start >= a.rt_start - 1) & (b.rt_end <= a.rt_end + 1)
        common = [p for p in (accA & accB & dn) if abs(A.at[p, 'rt_apex'] - d.at[p, 'rt']) <= 6]
        sp = split_flag(common)
        ca, cb = cov(A, common), cov(B, common)
        la, lb = np.log(A.area.reindex(common).clip(lower=1)), np.log(B.area.reindex(common).clip(lower=1))
        lq = np.log(d['Precursor.Quantity'].reindex(common).clip(lower=1))
        rho_a, rho_b = spearmanr(la, lq)[0], spearmanr(lb, lq)[0]
        spl = [p for p, f in zip(common, sp) if f]
        lost = list((accA & dn) - accB); gained = list((accB & dn) - accA)
        lost_split = split_flag(lost).mean() if lost else np.nan
        g = pd.DataFrame({'b': B.rt_apex.reindex(gained).values, 'a': A.rt_apex.reindex(gained).values, 'dn': d.rt.reindex(gained).values})
        true_rec = ((g.b - g.dn).abs() <= 6) & ~((g.a - g.dn).abs() <= 6)
        row = dict(run=R, lvl=L, arm=ARM, common=len(common), split_share=sp.mean(),
                   cov_def_split=ca[sp.values].median() if len(spl) else np.nan, cov_arm_split=cb[sp.values].median() if len(spl) else np.nan,
                   cov_def_all=ca.median(), cov_arm_all=cb.median(),
                   area_ratio_split=(B.area.reindex(spl) / A.area.reindex(spl)).median() if spl else np.nan,
                   rho_def=rho_a, rho_arm=rho_b, d_rho=rho_b - rho_a,
                   lost=len(lost), lost_split_share=lost_split, gained=len(gained), gained_true_share=true_rec.mean() if len(g) else np.nan,
                   width_def=(A.rt_end - A.rt_start).reindex(common).median(), width_arm=(B.rt_end - B.rt_start).reindex(common).median(),
                   width_dn=(d.s1 - d.s0).reindex(common).median())
        row['ruleA'] = bool(row['split_share'] <= 0.05 and row['d_rho'] >= -0.01 and (row['lost_split_share'] <= 0.25) and row['gained_true_share'] >= 0.60)
        out.append(row); print(row, flush=True)
pd.set_option('display.width', 250)
print(pd.DataFrame(out).round(4).to_string(index=False))
