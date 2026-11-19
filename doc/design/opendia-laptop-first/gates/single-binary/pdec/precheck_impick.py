"""PRECHECK_IMPICK (diagnostic): would RT/IM (2D) peak picking choose the right peak more often?
Data: r~1.56 HeLa D seed 0, current default (evidence block), -features_out (ev/D_new_s0): every candidate window with
x_im_centroid_offset (|mobility centroid - library/calibrated 1/K0| in bins of 0.0075), x_im_central_share, x_im_coelution.
Truth for evaluation only: DIA-NN 2.2 RT of non-entrapment q <= 0.01 precursors (r~7 D report, same raw file, same
peptides by modified sequence + charge). A window is 'right' if its apex is within 6 s of DIA-NN's RT.
Questions: (1) for DIA-NN IDs, is the right peak among the candidates at all? (2) when the kept window is wrong but a
right candidate exists, is the right one more mobility-consistent? (3) oracle rules on the existing candidates."""
import glob, sqlite3, numpy as np, pandas as pd
W = '/home/singjc/.jcode/scratch/dia2/gap7/ev/D_new_s0'
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
F = pd.read_csv(glob.glob(f'{W}/*.features.tsv')[0], sep='\t')
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg['key'] = tg.ms + '/' + tg.z.astype(str); lab = tg[tg.e == 0].drop_duplicates('key').set_index('key')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'IM', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]
d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[d.key.isin(lab.index)].drop_duplicates('key'); d['pid'] = lab.pid.reindex(d.key).values
d = d.set_index('pid')
T = F[F.is_decoy == 0].copy()
T = T[T.precursor_id.isin(d.index)]
T['drt'] = np.abs(T.rt_apex - (d.RT * 60).reindex(T.precursor_id).values)
T['right'] = T.drt <= 6
T['kept'] = T.score > -1e29
print(f'DIA-NN non-entrapment q<=0.01 IDs mapped into this library: {len(d)}; with an OpenDIA candidate: {T.precursor_id.nunique()}')
g = T.groupby('precursor_id')
s = pd.DataFrame({'n': g.size(), 'any_right': g.right.any(),
                  'kept_right': g.apply(lambda x: bool((x.kept & x.right).any())),
                  'q': g.apply(lambda x: x.loc[x.kept, 'q_value'].min() if x.kept.any() else 1.0)})
s['acc'] = s.q <= 0.01
s['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].reindex(s.index).clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print(f'candidates per precursor median {s.n.median():.0f}')
cat = np.select([s.acc & s.kept_right, s.acc & ~s.kept_right, ~s.acc & s.kept_right, ~s.acc & s.any_right & ~s.kept_right],
                ['accepted, right peak', 'accepted, other peak', 'missed: right peak kept, scored low',
                 'missed: right peak a candidate but not kept'], 'missed: right peak not a candidate')
s['cat'] = cat
print('\nOutcome of DIA-NN IDs (share), overall and by DIA-NN quantity quartile:')
tab = pd.crosstab(s.cat, s.qq, normalize='columns').round(3); tab['all'] = s.cat.value_counts(normalize=True).round(3)
print(tab.to_string())

# (2) Choice errors: kept wrong, right candidate exists. Compare the kept window with the best right candidate.
ce = s[s.any_right & ~s.kept_right].index
K = T[T.precursor_id.isin(ce) & T.kept].set_index('precursor_id')
R = T[T.precursor_id.isin(ce) & T.right].sort_values('drt').drop_duplicates('precursor_id').set_index('precursor_id')
K, R = K.loc[R.index.intersection(K.index)], R.loc[R.index.intersection(K.index)]
print(f'\nChoice errors (kept window wrong, right peak among candidates): {len(R)}')
for col, better in (('x_im_centroid_offset', 'lower'), ('x_im_central_share', 'higher'), ('x_im_coelution', 'higher'),
                    ('abs_library_rt_error', 'lower'), ('f12_libcos', 'higher'), ('x_ms2_log_area', 'higher')):
    r, k = R[col].values, K[col].values
    win = np.mean(r < k) if better == 'lower' else np.mean(r > k)
    print(f'  {col:22s} right p50 {np.median(r):7.3f}  kept-wrong p50 {np.median(k):7.3f}  right is {better} in {win:.2f} of pairs (ties {np.mean(r == k):.2f})')
# Kept-wrong windows that look like a mobility-offset interferer: centroid > 2 bins (0.015 1/K0) off centre.
print(f'  kept-wrong windows with mobility centroid > 2 bins off centre: {np.mean(K.x_im_centroid_offset > 2):.2f}; right windows: {np.mean(R.x_im_centroid_offset > 2):.2f}')

# (3) Oracle-free rules on the existing candidates: how often does each rule pick a right window (among precursors that
# have one), compared with the tool's kept window?
H = T[T.precursor_id.isin(s[s.any_right].index)]
def pick(col, asc):
    return H.sort_values(col, ascending=asc).drop_duplicates('precursor_id').set_index('precursor_id').right.mean()
print(f'\nAmong DIA-NN IDs with a right candidate ({s.any_right.sum()}): share where the rule picks a right window')
print(f'  tool (kept window)               {s[s.any_right].kept_right.mean():.3f}')
print(f'  smallest mobility centroid offset {pick("x_im_centroid_offset", True):.3f}')
print(f'  largest central mobility share    {pick("x_im_central_share", False):.3f}')
print(f'  smallest |library RT error|        {pick("abs_library_rt_error", True):.3f}')
print(f'  largest MS2 area                   {pick("x_ms2_log_area", False):.3f}')
# Mobility of the kept-wrong peak relative to DIA-NN: offset in the evidence (bins) for wrong vs right.
print(f'\nAll candidate windows of DIA-NN IDs: right windows with centroid offset <= 1 bin {np.mean(H[H.right].x_im_centroid_offset <= 1):.2f}, '
      f'wrong windows {np.mean(H[~H.right].x_im_centroid_offset <= 1):.2f}')
