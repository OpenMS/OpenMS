"""Test: are accepted entrapment winners sitting on real eluting peptides (same precursor m/z within the isolation
window and same RT), more often than accepted decoy winners? Real peptides = the ProteoBench DIA-NN reference IDs of
the same raw (8.65M library, q <= 0.01, any target, includes peptides not in our 100K library). Evaluation only.

For each pair winner above the q <= 0.10 cut (and q <= 0.01), count reference precursors with
|m/z - winner m/z| <= 10 ppm (same mass, so same isolation and same MS1 isotopes) and |RT - winner apex| <= 5 s.
Also report: share whose composition equals a reference peptide's composition (shuffle parent identified)."""
import sqlite3, sys, numpy as np, pandas as pd
sys.argv = [sys.argv[0], 'D11FWT', '../xicids-20261112', 'xic']
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
R = 148983 / 151239
exec(open('diagnose.py').read().split("for raw in ('raw01', 'raw02'):")[0])  # helpers, labels, tdc
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge', 'Precursor.Mz', 'RT', 'Q.Value', 'Decoy', 'Protein.Ids'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0)]
seq = pd.read_sql("""select p.ID precursor_id, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID""", c).set_index('precursor_id').seq
comp = lambda s: ''.join(sorted(s))
for raw, run in (('raw01', '_01'), ('raw02', '_02')):
    rr = ref[ref.Run.str.endswith(run)].copy()
    rr['RTs'] = rr.RT * 60.0
    rr = rr.sort_values('Precursor.Mz')
    ref_mz = rr['Precursor.Mz'].values; ref_rt = rr.RTs.values
    ref_comp = set(rr['Stripped.Sequence'].map(comp))
    ref_seq = set(rr['Stripped.Sequence'])
    pop = pd.read_csv(f'{ROOT}/{raw}/population.tsv', sep='\t')
    s = np.load(f'{ROOT}/{raw}/scores_{ARM}.npy').astype(np.float64)
    tgt_of_pair = pop[pop.is_decoy == 0].drop_duplicates('pair_key').set_index('pair_key').precursor_id
    pop['cls'] = np.where(pop.is_decoy == 0, pop.precursor_id.map(ent), 2 + pop.pair_key.map(tgt_of_pair).map(ent).fillna(-9))
    w, best = tdc(pop, s)
    best = best.join(pre[['charge', 'mz']], on='precursor_id')
    best['seq'] = best.precursor_id.map(seq)
    def near(mz, rt, ppm=10.0, drt=5.0):
        lo = np.searchsorted(ref_mz, mz * (1 - ppm * 1e-6)); hi = np.searchsorted(ref_mz, mz * (1 + ppm * 1e-6))
        return int(np.any(np.abs(ref_rt[lo:hi] - rt) <= drt))
    best['on_real'] = [near(m, r) for m, r in zip(best.mz.values, best.rt_apex.values)]
    best['on_real_mz'] = [near(m, r, 10.0, 1e9) for m, r in zip(best.mz.values, best.rt_apex.values)]
    best['comp_ref'] = best.seq.map(lambda x: comp(x) in ref_comp).astype(int)
    best['in_ref'] = best.seq.isin(ref_seq).astype(int)
    names = {0: 'human', 1: 'entrap', 2: 'dec(hum)', 3: 'dec(ent)'}
    print(f'==== {raw}: reference IDs {len(rr)}')
    for q in (0.01, 0.05, 0.10):
        cut = best[best.q <= q].s.min(); top = best[best.s >= cut]
        print(f'  q <= {q:.2f}: share of winners  [same m/z+RT as a reference ID | same m/z any RT | composition = a reference peptide | sequence in reference]')
        for k_, g in top.groupby('cls'):
            if k_ in names:
                print(f'    {names[k_]:9s} n {len(g):5d}  {g.on_real.mean():.3f} | {g.on_real_mz.mean():.3f} | {g.comp_ref.mean():.3f} | {g.in_ref.mean():.3f}')
        # Calibration split by whether the winner sits on a real peptide.
        for flag, g in top.groupby('on_real'):
            E = int((g.cls == 1).sum()); D = int((g.is_decoy == 1).sum()); N = int((g.is_decoy == 0).sum())
            print(f'    on_real={flag}:' + fmt(E, D, N))
    # Background: all population precursors (no scores).
    allp = pop.drop_duplicates('precursor_id').join(pre[['mz']], on='precursor_id')
    allp['on_real_mz'] = [near(m, 0, 10.0, 1e9) for m in allp.mz.values]
    print('  background, all population precursors, same m/z (10 ppm) as a reference ID:',
          allp.groupby('cls').on_real_mz.mean().round(3).to_dict())
