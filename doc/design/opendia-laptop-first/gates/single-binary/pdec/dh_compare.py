"""D vs H (diagnostic): why does the matched filter gain 6.5 % on D and 2.0 % on H?
Per run: DIA-NN peak width and run length, OpenDIA cycle (MS2 spectra per isolation window), RT bins per peak,
calibration (RT tolerance, anchors), and the A / B1 / B2 / C funnel of DIA-NN IDs for CUR and MT (seed 0).
Uses only the tool outputs as written; DIA-NN labels for evaluation only."""
import sys, json, glob, sqlite3, numpy as np, pandas as pd
C = '/nfs/research/juan/singj/opendia'
LIB = f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid')
tg = pq[pq.d == 0]; sk = (tg.ms + '/' + tg.z.astype(str)).reset_index().drop_duplicates(0).set_index(0).pid
LOGS = {'D': (735, 8635), 'H': (None, 32966)}
for run in ('D', 'H'):
    d = pd.read_parquet(f'{C}/entlib-diann5-20261119/{run}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Decoy', 'RT', 'FWHM', 'RT.Start', 'RT.Stop', 'Precursor.Quantity'])
    d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy()
    d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
    d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
    d = d[pq.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid'); d['rt_s'] = d.RT * 60
    fw = d.FWHM * 60; base = (d['RT.Stop'] - d['RT.Start']) * 60
    rt_span = (d.rt_s.max() - d.rt_s.min()) / 60
    print(f'\n===== {run}: DIA-NN non-entrapment IDs {len(d)}, elution span {rt_span:.1f} min')
    print(f'  DIA-NN FWHM s median {fw.median():.2f} (p10 {fw.quantile(.1):.2f}, p90 {fw.quantile(.9):.2f}); base width s median {base.median():.2f}')
    for arm in ('cur', 'mt'):
        rd = f'{C}/entlib-mt-20261120/runs/{run}_{arm}_s0'
        cal = json.load(open(glob.glob(f'{rd}/*.calibration.json')[0]))
        P = pd.read_csv(glob.glob(f'{rd}/*.precursors.tsv')[0].replace('precursors.tsv', 'precursors.tsv'), sep='\t',
                        usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'rt_start', 'rt_end', 'score', 'q_value'])
        P = P[[not x for x in P.precursor_id.isna()]]
        allw = P.copy()
        P = P[P.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        ours = set(P[(P.is_decoy == 0) & (P.q_value <= 0.01)].index)
        if arm == 'cur':
            print(f'  calibration: RT tolerance {cal["rt"]["tolerance_half_width_s"]:.1f} s, anchors {cal["rt"]["anchors_used"]}, '
                  f'RT residual p50/p95 {cal["rt"]["abs_residual_p50_s"]:.1f}/{cal["rt"]["abs_residual_p95_s"]:.1f} s, '
                  f'MS2 tol {cal["ms2_mz"]["tolerance_half_width_ppm"]:.1f} ppm, IM sigma {cal["im"]["sigma_vssc"]:.4f}')
            # search peak widths and how many 0.94 s tensor bins they span
            pw = (P.rt_end - P.rt_start)
            acc = P.loc[list(ours)]
            print(f'  OpenDIA accepted peak width s median {(acc.rt_end - acc.rt_start).median():.2f}; tensor bin 0.94 s -> '
                  f'{(acc.rt_end - acc.rt_start).median() / 0.9375:.1f} bins per accepted peak')
        # funnel: right window = window of that precursor closest to DIA-NN RT within 6 s
        w = allw[allw.precursor_id.isin(d.index)].copy(); w['drt'] = (w.rt_apex - d.rt_s.reindex(w.precursor_id).values).abs()
        right = w.sort_values('drt').drop_duplicates('precursor_id').set_index('precursor_id'); right = right[right.drt <= 6]
        cls = pd.Series('C', index=d.index); hr = d.index.isin(right.index); cls[hr] = 'B2'
        br = P.candidate_row.reindex(d.index).values
        cls[hr & (br == right.candidate_row.reindex(d.index).values)] = 'B1'; cls[d.index.isin(ours)] = 'A'
        vc = cls.value_counts().reindex(['A', 'B1', 'B2', 'C']).fillna(0).astype(int)
        print(f'  {arm}: ' + '  '.join(f'{k} {v} ({v / len(d):.3f})' for k, v in vc.items()))
        if arm == 'cur': cur_cls = cls
        else:
            moved = pd.crosstab(cur_cls, cls)
            print('  CUR class -> MT class:\n' + moved.to_string())
            # where do B1 sit in DIA-NN quantity: faint share per run
            q = np.log10(d['Precursor.Quantity'].clip(lower=1))
            for k in ('A', 'B1'):
                print(f'  {k} (CUR) DIA-NN log10 quantity median {q[cur_cls == k].median():.2f}; FWHM s {fw[cur_cls == k].median():.2f}')
