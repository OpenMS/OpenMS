"""LIBSIZE (diagnostic, 21 Nov 2026): how much does the entrapment ratio r cost OpenDIA (current default) against DIA-NN?

Nested subsets of the r~7 library: every HeLa (non-entrapment) target is kept; an entrapment target peptide is kept when
crc32(unmodified sequence) % 10000 < f x 10000, so r1 is a subset of r3, which is a subset of r7. Each kept target keeps its
paired DECOY_ precursor (OpenDIA). DIA-NN gets the matching targets-only TSV (it makes its own decoys).
Levels: r0 (f = 0, no entrapment: N only, no FDP estimate), r1 (f = 1/7.0539), r3 (f = 3/7.0539); r7 = the existing runs.

Usage:  libsize.py build <out_root>
        libsize.py eval  <out_root>
"""
import os, sys, json, zlib, sqlite3, numpy as np, pandas as pd
C = '/nfs/research/juan/singj/opendia'
SRC = f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'
TSV = f'{C}/entlib-diann5-20261119/targets_r5.tsv'
R7 = 7.0539
LEVELS = {'r0': 0.0, 'r1': 1.0 / R7, 'r3': 3.0 / R7}
LABEL_SQL = """select p.ID pid, p.TRAML_ID traml, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, pe.UNMODIFIED_SEQUENCE us,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID"""

def labels(path):
    c = sqlite3.connect(f'file:{path}?mode=ro', uri=True); q = pd.read_sql(LABEL_SQL, c); c.close(); return q

def build(root):
    q = labels(SRC)
    tg = q[q.d == 0].copy()
    tg['h'] = [zlib.crc32(s.encode()) % 10000 for s in tg.us]
    print(f'source: precursors {len(q)}, targets {len(tg)}, entrapment targets {int(tg.e.sum())}', flush=True)
    meta = {}
    for name, f in LEVELS.items():
        out = f'{root}/{name}'; os.makedirs(out, exist_ok=True)
        keep_t = tg[(tg.e == 0) | (tg.h < f * 10000)]
        keep_traml = set(keep_t.traml)
        dec = q[(q.d == 1) & q.traml.str.replace('^DECOY_', '', regex=True).isin(keep_traml)]
        ids = list(keep_t.pid) + list(dec.pid)
        n_e, n_h = int(keep_t.e.sum()), int((keep_t.e == 0).sum())
        meta[name] = dict(targets=len(keep_t), hela=n_h, entrapment=n_e, decoys=len(dec), r=n_e / n_h)
        print(f'{name}: targets {len(keep_t)} (HeLa {n_h}, entrapment {n_e}, r {n_e / n_h:.4f}), decoys {len(dec)}', flush=True)
        # PQP subset
        dst = f'{out}/lib.pqp'
        if os.path.exists(dst): os.remove(dst)
        c = sqlite3.connect(dst, uri=True); c.execute(f"attach database 'file:{SRC}?mode=ro' as s")
        tables = c.execute("select name, sql from s.sqlite_master where type='table'").fetchall()
        for _, sql in tables: c.execute(sql)
        c.execute('create temp table sel(id integer primary key)'); c.executemany('insert into sel values (?)', [(int(i),) for i in ids])
        c.execute('insert into PRECURSOR select * from s.PRECURSOR where ID in (select id from sel)')
        c.execute('insert into TRANSITION_PRECURSOR_MAPPING select * from s.TRANSITION_PRECURSOR_MAPPING where PRECURSOR_ID in (select id from sel)')
        c.execute('insert into TRANSITION select * from s.TRANSITION where ID in (select TRANSITION_ID from TRANSITION_PRECURSOR_MAPPING)')
        c.execute('insert into PRECURSOR_PEPTIDE_MAPPING select * from s.PRECURSOR_PEPTIDE_MAPPING where PRECURSOR_ID in (select id from sel)')
        c.execute('insert into PEPTIDE select * from s.PEPTIDE where ID in (select PEPTIDE_ID from PRECURSOR_PEPTIDE_MAPPING)')
        c.execute('insert into TRANSITION_PEPTIDE_MAPPING select * from s.TRANSITION_PEPTIDE_MAPPING where TRANSITION_ID in (select ID from TRANSITION)')
        c.execute('insert into PEPTIDE_PROTEIN_MAPPING select * from s.PEPTIDE_PROTEIN_MAPPING where PEPTIDE_ID in (select ID from PEPTIDE)')
        c.execute('insert into PROTEIN select * from s.PROTEIN where ID in (select PROTEIN_ID from PEPTIDE_PROTEIN_MAPPING)')
        c.execute('insert into PEPTIDE_GENE_MAPPING select * from s.PEPTIDE_GENE_MAPPING where PEPTIDE_ID in (select ID from PEPTIDE)')
        c.execute('insert into GENE select * from s.GENE where ID in (select GENE_ID from PEPTIDE_GENE_MAPPING)')
        for t in ('VERSION', 'COMPOUND', 'PRECURSOR_COMPOUND_MAPPING'): c.execute(f'insert into {t} select * from s.{t}')
        for (sql,) in c.execute("select sql from s.sqlite_master where type='index' and sql is not null").fetchall(): c.execute(sql)
        c.commit(); c.close()
        # DIA-NN targets TSV subset (same targets)
        keys = set(keep_t.ms + '/' + keep_t.z.astype(str))
        first = True; n_rows = 0
        with open(f'{out}/targets.tsv', 'w') as fo:
            for ch in pd.read_csv(TSV, sep='\t', dtype=str, keep_default_na=False, chunksize=2_000_000):
                k = ch.ModifiedPeptideSequence + '/' + ch.PrecursorCharge
                sub = ch[k.isin(keys)]; n_rows += len(sub)
                sub.to_csv(fo, sep='\t', index=False, header=first); first = False
        meta[name]['tsv_rows'] = n_rows
        print(f'{name}: PQP and TSV ({n_rows} transitions) written', flush=True)
    json.dump(meta, open(f'{root}/levels.json', 'w'), indent=1)

def eval_(root):
    from scipy.stats import chi2
    q = labels(SRC); tg = q[q.d == 0]
    ent = tg.set_index('pid').e
    keyed = (tg.ms + '/' + tg.z.astype(str)).values; lab = pd.Series(tg.e.values, index=keyed); lab = lab[~lab.index.duplicated()]
    pid_of = pd.Series(tg.pid.values, index=keyed); pid_of = pid_of[~pid_of.index.duplicated()]
    meta = json.load(open(f'{root}/levels.json'))
    meta['r7'] = dict(r=R7)
    runs = {}
    for run in ('D', 'H'):
        for lvl in ('r0', 'r1', 'r3', 'r7'):
            r = meta[lvl]['r']; k = 1 + 1 / r if r > 0 else np.nan
            # OpenDIA
            od_dirs = ([f'{root}/{lvl}/opendia/{run}_s{s}' for s in (0, 1)] if lvl != 'r7'
                       else [f'{C}/entlib-rtr-20261121/runs/{run}_rtr_s{s}' for s in (0, 1)])
            for s, dd in enumerate(od_dirs):
                f = [x for x in os.listdir(dd) if x.endswith('.precursors.tsv') and x != 'precursors.tsv']
                t = pd.read_csv(f'{dd}/{f[0]}', sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
                t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
                tt = t[t.is_decoy == 0].copy(); tt['e'] = tt.precursor_id.map(ent).fillna(0).values
                acc = tt[tt.q_value <= 0.01]; N, E = len(acc), int(acc.e.sum())
                row = dict(run=run, lvl=lvl, r=r, tool='OpenDIA', seed=s, N_q01=N, hela_q01=N - E, E=E,
                           thr=float(acc.score.min()) if N else np.nan)
                if r > 0:
                    e = tt.sort_values(['q_value', 'score'], ascending=[True, False]).e.values
                    fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
                    row.update(lower=E / N, combined=E * k / N, N_fdp1=int(ok.max() + 1) if len(ok) else 0)
                row['ids'] = set(acc[acc.e == 0].precursor_id)
                runs[(run, lvl, 'OpenDIA', s)] = row
            # DIA-NN
            rep = f'{root}/{lvl}/diann/{run}/report.parquet' if lvl != 'r7' else f'{C}/entlib-diann5-20261119/{run}/report.parquet'
            d = pd.read_parquet(rep, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'PEP'])
            d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
            d = d.sort_values(['Q.Value', 'PEP']).drop_duplicates('key'); d = d[d.key.isin(lab.index)].copy()
            d['e'] = lab.reindex(d.key).values; d['pid'] = pid_of.reindex(d.key).values
            acc = d[d['Q.Value'] <= 0.01]; N, E = len(acc), int(acc.e.sum())
            row = dict(run=run, lvl=lvl, r=r, tool='DIA-NN', seed=0, N_q01=N, hela_q01=N - E, E=E, thr=np.nan)
            if r > 0:
                fdp = np.cumsum(d.e.values) * k / np.arange(1, len(d) + 1); ok = np.where(fdp <= 0.01)[0]
                row.update(lower=E / N, combined=E * k / N, N_fdp1=int(ok.max() + 1) if len(ok) else 0)
            row['ids'] = set(acc[acc.e == 0].pid)
            runs[(run, lvl, 'DIA-NN', 0)] = row
    rows = [{k: v for k, v in x.items() if k != 'ids'} for x in runs.values()]
    df = pd.DataFrame(rows)
    pd.set_option('display.width', 220)
    print(df.round(4).to_string(index=False))
    print('\n== per tool and level (OpenDIA: mean of 2 seeds): HeLa IDs at own q<=0.01, N at true FDP 1 %, share of r0 HeLa IDs')
    g = df.groupby(['run', 'tool', 'lvl']).agg(r=('r', 'first'), hela_q01=('hela_q01', 'mean'), N_fdp1=('N_fdp1', 'mean'),
                                               lower=('lower', 'mean'), combined=('combined', 'mean'), thr=('thr', 'mean')).reset_index()
    for (run, tool), x in g.groupby(['run', 'tool']):
        base = x[x.lvl == 'r0'].hela_q01.values[0]
        x = x.assign(share_of_r0=x.hela_q01 / base)
        print(f'\n{run} {tool}\n' + x.round(4).to_string(index=False))
    print('\n== OpenDIA / DIA-NN, HeLa IDs at own q<=0.01 and N at true FDP 1 %')
    for run in ('D', 'H'):
        for lvl in ('r0', 'r1', 'r3', 'r7'):
            o = g[(g.run == run) & (g.tool == 'OpenDIA') & (g.lvl == lvl)].iloc[0]; n = g[(g.run == run) & (g.tool == 'DIA-NN') & (g.lvl == lvl)].iloc[0]
            s = f'  {run} {lvl}: own-q HeLa {o.hela_q01:.0f} / {n.hela_q01:.0f} = {o.hela_q01 / n.hela_q01:.3f}'
            if lvl != 'r0': s += f';  N@FDP_c 1 % {o.N_fdp1:.0f} / {n.N_fdp1:.0f} = {o.N_fdp1 / n.N_fdp1:.3f}'
            print(s)
    print('\n== retention: of the r0 HeLa IDs (seed 0), share still identified at each level, and overlap with DIA-NN r0 IDs')
    for run in ('D', 'H'):
        for tool in ('OpenDIA', 'DIA-NN'):
            b = runs[(run, 'r0', tool, 0)]['ids']
            print(f'  {run} {tool}: ' + '  '.join(f'{lvl} {len(b & runs[(run, lvl, tool, 0)]["ids"]) / max(len(b), 1):.3f}' for lvl in ('r1', 'r3', 'r7')))
        o, n = runs[(run, 'r0', 'OpenDIA', 0)]['ids'], runs[(run, 'r0', 'DIA-NN', 0)]['ids']
        print(f'  {run} r0: OpenDIA IDs in DIA-NN {len(o & n) / len(o):.3f}; DIA-NN IDs in OpenDIA {len(o & n) / len(n):.3f}')

if __name__ == '__main__':
    {'build': build, 'eval': eval_}[sys.argv[1]](sys.argv[2])
