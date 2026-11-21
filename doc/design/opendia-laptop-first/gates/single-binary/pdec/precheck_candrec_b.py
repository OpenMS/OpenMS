"""PRECHECK_CANDREC part B evaluation (diagnostic, 21 Nov 2026).
For each run (D, H) and level (r0, r7), compares arms def / V8M8 / M10:
  recovery  share of the part-A C IDs (candrec_a_<R>_<L>.tsv from the def run) with a target window within 6 s of DIA-NN RT
  HeLa IDs at own q<=0.01 (from <run>.precursors.tsv), N@FDP_c 1 % and bounds at r7, wall time and peak RSS from the log,
  windows written, and how many C IDs are accepted at q<=0.01 in the arm.
Usage: precheck_candrec_b.py <ceil0 dir> <candrec dir>"""
import sys, os, glob, re, sqlite3, numpy as np, pandas as pd
W0, W = sys.argv[1:3]
C = '/nfs/research/juan/singj/opendia'
LIB = {'r0': f'{C}/libsize-20261121/r0/lib.pqp', 'r7': f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'}
DN = {('D', 'r0'): f'{C}/libsize-20261121/r0/diann/D/report.parquet', ('H', 'r0'): f'{C}/libsize-20261121/r0/diann/H/report.parquet',
      ('D', 'r7'): f'{C}/entlib-diann5-20261119/D/report.parquet', ('H', 'r7'): f'{C}/entlib-diann5-20261119/H/report.parquet'}
R7 = 7.0539; k = 1 + 1 / R7
def wall(log):
    t = open(log).read(); m = re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([0-9:.]+)', t)
    rss = re.findall(r'Maximum resident set size \(kbytes\): (\d+)', t)
    p = [float(x) for x in m[-1].split(':')]; s = p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) > 2 else 0)
    return s, int(rss[-1]) / 1e6
labs = {}
def labels(L):
    if L in labs: return labs[L]
    c = sqlite3.connect(f'file:{LIB[L]}?mode=ro', uri=True)
    tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
      min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
      join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
    labs[L] = tg; return tg
rows = []
for R in ('D', 'H'):
    for L in ('r0', 'r7'):
        tg = labels(L); ent = tg.set_index('pid').e
        A = pd.read_csv(f'{W0}/{R}_{L}/candrec_a_{R}_{L}.tsv', sep='\t', index_col=0)
        d = pd.read_parquet(DN[(R, L)], columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT'])
        d = d[d['Q.Value'] <= 0.01]; t2 = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
        kk = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
        d['pid'] = t2.pid.reindex(kk).values; d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
        d = d.sort_values('Q.Value').drop_duplicates('pid').set_index('pid'); dn_rt = d.RT * 60
        for arm in ('def', 'V8M8', 'M10'):
            rd = f'{W0}/{R}_{L}' if arm == 'def' else f'{W}/{R}_{L}_{arm}'
            log = f'{W0}/{R}_{L}.log' if arm == 'def' else f'{W}/{R}_{L}_{arm}.log'
            if not os.path.exists(log) or 'exit 0' not in open(log).read(): print(f'{R} {L} {arm}: missing'); continue
            pt = [f for f in glob.glob(f'{rd}/*.precursors.tsv') if not f.endswith('/precursors.tsv')][0]
            t = pd.read_csv(pt, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
            t = t[(t.is_decoy == 0) & (t.score > -1e29)].sort_values('score', ascending=False).drop_duplicates('precursor_id')
            t['e'] = t.precursor_id.map(ent).fillna(0).values
            acc = t[t.q_value <= 0.01]; N, E = len(acc), int(acc.e.sum())
            row = dict(run=R, lvl=L, arm=arm, hela=N - E)
            if L == 'r7':
                e = t.sort_values(['q_value', 'score'], ascending=[True, False]).e.values
                fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
                row.update(lower=E / N, combined=E * k / N, n_fdp1=int(ok.max() + 1))
            ft = glob.glob(f'{rd}/*.features.tsv')[0]
            Fw = pd.read_csv(ft, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex'])
            row['windows'] = len(Fw)
            Ft = Fw[(Fw.is_decoy == 0) & Fw.precursor_id.isin(A.index)].copy()
            Ft['drt'] = (Ft.rt_apex - dn_rt.reindex(Ft.precursor_id).values).abs()
            rec = Ft.groupby('precursor_id').drt.min().reindex(A.index) <= 6
            row['recovered'] = rec.mean()
            for cat in ('near (6-15 s)', 'picker (<5 windows)', 'cap (5 windows)'):
                m = A.cat == cat
                row[f'rec_{cat.split()[0]}'] = rec[m].mean() if m.any() else np.nan
            row['C_accepted'] = int(A.index.isin(acc.precursor_id[acc.e == 0]).sum())
            dnh = set(d.index); row['dn_cov'] = len(set(acc.precursor_id[acc.e == 0]) & dnh) / len(dnh)
            row['wall_s'], row['rss_gb'] = wall(log)
            rows.append(row); print(row, flush=True)
df = pd.DataFrame(rows)
for (R, L), g in df.groupby(['run', 'lvl']):
    b = g[g.arm == 'def'].iloc[0]
    df.loc[g.index, 'd_hela'] = g.hela / b.hela - 1; df.loc[g.index, 'wall_x'] = g.wall_s / b.wall_s
    if L == 'r7': df.loc[g.index, 'd_fdp1'] = g.n_fdp1 / b.n_fdp1 - 1; df.loc[g.index, 'd_lower_pp'] = (g.lower - b.lower) * 100
pd.set_option('display.width', 250)
print(df.round(4).to_string(index=False))
