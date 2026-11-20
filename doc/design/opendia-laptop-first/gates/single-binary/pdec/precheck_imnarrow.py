"""PRECHECK_IMNARROW: per-quartile TPR of DIA-NN targets' right window (tool score) at fixed decoy FPR, default vs a
narrower mobility arm. r~1.56 D s0. Usage: precheck_imnarrow.py <arm> [<arm> ...] (dirs under gap7/ev with features).
Labels for evaluation only."""
import glob, sys, sqlite3, numpy as np, pandas as pd
G = '/home/singjc/.jcode/scratch/dia2/gap7/ev'
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
lab = tg[tg.e == 0].assign(k=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('k').set_index('k')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]; k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[k.isin(lab.index).values]; d['pid'] = lab.pid.reindex(k[k.isin(lab.index)]).values; d = d.drop_duplicates('pid').set_index('pid')
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
cols = ['precursor_id', 'is_decoy', 'rt_apex', 'score', 'x_im_centroid_offset', 'x_im_central_share', 'f12_libcos', 'boundary_coelution']
for arm in sys.argv[1:]:
    F = pd.read_csv(glob.glob(f'{G}/{arm}/*.features.tsv')[0], sep='\t', usecols=cols)
    s = np.where(F.score > -1e29, F.score, -np.inf); F['s'] = s
    neg = F[F.is_decoy == 1].groupby('precursor_id').s.max(); neg = neg[np.isfinite(neg)].values
    T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)]
    T = T[np.abs(T.rt_apex - d.RT.reindex(T.precursor_id).values * 60) <= 6]
    R = T.sort_values('s', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    qq = d.qq.reindex(R.index)
    line = []
    for fpr in (1e-4, 1e-3, 1e-2):
        thr = np.quantile(neg, 1 - fpr)
        # share of ALL DIA-NN IDs (missing right window counts as a miss)
        line.append(f'FPR {fpr:g}: ' + ' '.join(f'{q} {(R.s[qq == q] > thr).sum() / (d.qq == q).sum():.3f}' for q in ['Q1', 'Q2', 'Q3', 'Q4'])
                    + f' all {(R.s > thr).sum() / len(d):.3f}')
    print(f'{arm:9s} windows {len(F):8d} | ' + ' | '.join(line))
    Rq = R[qq.isin(['Q1', 'Q2']).values]
    print(f'          Q1-Q2 right windows: median f12_libcos {Rq.f12_libcos.median():.3f}, boundary_coelution {Rq.boundary_coelution.median():.3f}, '
          f'im_central_share {Rq.x_im_central_share.median():.3f}')
