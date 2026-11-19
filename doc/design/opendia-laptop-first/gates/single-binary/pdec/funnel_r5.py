import sqlite3, sys, numpy as np, pandas as pd
L = '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg['key'] = tg.ms + '/' + tg.z.astype(str); lab = tg.drop_duplicates('key').set_index('key')
for R in ('D', 'H'):
    d = pd.read_parquet(f'{R}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
    d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
    d = d[d['Q.Value'] <= 0.01].drop_duplicates('key'); d = d[lab.e.reindex(d.key).values == 0]
    d['pid'] = lab.pid.reindex(d.key).values
    o = pd.read_csv(f'/home/singjc/.jcode/scratch/dia2/entlib10/cut/{R}_cur_s0.tsv.gz', sep='\t').set_index('precursor_id')
    o = o[o.score > -1e29]
    d['in'] = d.pid.isin(o.index); d['oq'] = o.q_value.reindex(d.pid).values; d['ort'] = o.rt_apex.reindex(d.pid).values
    d['drt'] = (d.ort - d.RT * 60).abs()
    n = len(d); acc = d.oq <= 0.01
    print(f'== {R}: DIA-NN non-entrapment q<=0.01 {n}')
    print(f'  no OpenDIA window at all (not a search candidate / dropped by population) {(~d["in"]).mean():.3f}')
    m = d[d['in'] & ~acc]
    print(f'  scored but not accepted {len(m) / n:.3f}; of these best window within 6 s of DIA-NN RT {(m.drt <= 6).mean():.3f}, q<=0.05 {(m.oq <= 0.05).mean():.3f}, q<=0.2 {(m.oq <= 0.2).mean():.3f}')
    print(f'  accepted {acc.mean():.3f} (right RT {(d[acc].drt <= 6).mean():.3f})')
    d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1 low', 'Q2', 'Q3', 'Q4 high'])
    print('  accepted by DIA-NN quantity quartile:', d.groupby('qq', observed=True).apply(lambda g: round((g.oq <= 0.01).mean(), 3)).to_dict())
