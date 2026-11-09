"""Hypothesis-space check (diagnostic, raw 01, native byion arm). For OSW+Percolator human IDs (m_score<=0.01, rank 1,
reference only), and for DIA-NN human IDs (Q.Value<=0.01), report against the native candidates:
  (a) share whose reference RT lies inside our calibrated RT window (|RT - map(library RT)| <= tolerance_half_width_s);
  (b) share with a native candidate apex within 10 s of the reference RT (any candidate, in window);
  (c) share where our precursor winner is that candidate;
  (d) native candidates per precursor in window.
Reading: if (a) or (b) is far below 1 for the OSW IDs we miss, the window or candidate generation limits IDs."""
import json, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
LIB = B + '/astral-entrapment-20261107/lib/library.pqp'; R = 148983 / 151239
cal = json.load(open(D + '/window300k/cal.json')); rt = cal['rt']; tol = rt['tolerance_half_width_s']
kx, ky = np.array(rt['knots_library_rt']), np.array(rt['knots_run_rt_s'])
def rtmap(x):
    return np.interp(x, kx, ky, left=np.nan, right=np.nan) if False else np.interp(x, kx, ky) + np.where(x < kx[0], (x - kx[0]) * (ky[1] - ky[0]) / (kx[1] - kx[0]), 0) + np.where(x > kx[-1], (x - kx[-1]) * (ky[-1] - ky[-2]) / (kx[-1] - kx[-2]), 0)
c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
p = pd.read_sql('select p.ID precursor_id, p.CHARGE z, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0', c).drop_duplicates('precursor_id')
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
p = p.merge(ent, on='precursor_id'); p['key'] = p.ms + '/' + p.z.astype(str); p['pred'] = rtmap(p.lrt.values)
cr = pd.read_csv(D + '/stage2-byion-design/byion/crossfit.tsv', sep='\t')
win = pd.read_csv(D + '/pipeline/all_windows.tsv', sep='\t', usecols=['sample_id', 'rt_apex'])
cr = cr.merge(win, on='sample_id'); cr = cr[(cr.in_rt_window == 1) & (cr.is_decoy == 0)]
best = cr[cr.precursor_best == 1].set_index('precursor_id')
# our accepted set at FDP_c 1% (tool stage2_score ranking, entrapment-matched)
w = best.join(p.set_index('precursor_id')[['e']], how='inner').sort_values('stage2_score', ascending=False)
fdp = np.cumsum(w.e.values) * (1 + 1 / R) / np.arange(1, len(w) + 1); acc = set(w.index[:int(np.where(fdp <= 0.01)[0].max() + 1)])
by = cr.groupby('precursor_id').rt_apex.apply(np.array).to_dict()
print(f'RT window half width {tol:.1f} s (OSW: +-300 s). Native in-window candidates per target precursor: mean {cr.groupby("precursor_id").size().mean():.2f}')
def report(name, ref):
    ref = ref.merge(p[['key', 'precursor_id', 'pred']], on='key')
    a = np.abs(ref.RT - ref.pred) <= tol
    near = np.array([pid in by and np.any(np.abs(by[pid] - t) <= 10) for pid, t in zip(ref.precursor_id, ref.RT)])
    won = np.array([pid in best.index and abs(best.loc[pid, 'rt_apex'] - t) <= 10 for pid, t in zip(ref.precursor_id, ref.RT)])
    isacc = ref.precursor_id.isin(acc).values
    for sub, m in (('all', np.ones(len(ref), bool)), ('we accept', isacc), ('we miss', ~isacc)):
        print(f'{name:9s} {sub:9s} n={m.sum():5d} | in our RT window {a[m].mean():.3f} | candidate within 10 s {near[m].mean():.3f} | our winner is it {won[m].mean():.3f}')
    miss = ~isacc
    print(f'{name:9s} misses split: outside window {(~a & miss).sum()}, in window but no candidate {(a & ~near & miss).sum()}, candidate but other winner {(near & ~won & miss).sum()}, right winner below cut {(won & miss).sum()}')
o = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t',
                usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'decoy', 'm_score', 'ProteinName', 'peak_group_rank'])
o = o[o.run_name.str.endswith('_01') & (o.decoy == 0) & (o.m_score <= 0.01) & (o.peak_group_rank == 1)]
o = o[~o.ProteinName.str.split(';').map(lambda a: all('_p_target' in x for x in a))]
report('OSW', pd.DataFrame({'key': o.FullPeptideName + '/' + o.Charge.astype(str), 'RT': o.RT}))
dn = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                     columns=['Run', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Protein.Ids'])
dn = dn[dn.Run.str.endswith('_01') & (dn['Q.Value'] <= 0.01)]
dn = dn[~dn['Protein.Ids'].astype(str).str.contains('_p_target')]
report('DIA-NN', pd.DataFrame({'key': dn['Modified.Sequence'] + '/' + dn['Precursor.Charge'].astype(str), 'RT': dn.RT * 60}))
