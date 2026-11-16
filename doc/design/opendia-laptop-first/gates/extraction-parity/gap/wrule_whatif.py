"""PRECHECK_PEAK part 2: window choice. For each precursor the W rule keeps the window with the largest linear part
z'w_z of its fold model. What-ifs (label-free; same fold models): keep the window with the largest full score h'w.
Then the same tail calibration and TDC (reimplemented in Python: the decoy-tail calibration is approximated by the
evaluation's per-fold decoy-z standardisation, as evaluate_v2.tdc already does). Reports, against DIA-NN as reference:
right-peak kept, and N@FDP_c 1 %."""
import sys, numpy as np, pandas as pd, sqlite3
RAW = sys.argv[1]
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'; B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
P = '/home/singjc/.jcode/scratch/e2e/peak/score'
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,pe.MODIFIED_SEQUENCE ms,p.DECOY decoy from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2)
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key')
REF = {'raw01': B + '/diann-pairlib-20261113/report_raw01.parquet', 'raw02': B + '/diann-decoys-20261114/own/raw02/report.parquet'}[RAW]
ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid'); ref['rt_s'] = ref.RT * 60
pop = pd.read_csv(f'{G}/armE_{RAW}/population.tsv', sep='\t'); pop['row'] = np.arange(len(pop))
dref = pop.precursor_id.map(ref.rt_s); right = (pop.rt_apex - dref).abs() <= 6
for s in range(3):
    kept = np.load(f'{P}/{RAW}_s{s}.npy').astype(float); raw = np.load(f'{P}/{RAW}_s{s}.npy.raw.npy').astype(float); lin = np.load(f'{P}/{RAW}_s{s}.npy.linear.npy').astype(float)
    out = []
    for name, sel in (('W rule (max linear)', lin), ('max full score', raw)):
        # keep one window per precursor by `sel`, score it with the full fold-model score
        o = pop.assign(sel=sel).sort_values(['sel', 'row'], ascending=[False, True]).drop_duplicates('precursor_id')
        sc = np.full(len(pop), -1e30); sc[o.row.values] = raw[o.row.values]
        w, best = tdc(pop, sc)
        tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
        tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values); n1q = int((tb.q <= 0.01).sum())
        Wk = w.set_index('precursor_id')
        hr = [p for p in ref.index if p in Wk.index]
        rk = int(np.sum(right.values[Wk.loc[hr, 'row'].astype(int).values]))
        out.append(f'{name}: right-peak kept {rk} | N@0.5/1/2% {[n_at(tw.s.values, tw.e.values, t) for t in (0.005, 0.01, 0.02)]} q<=0.01 {n1q}')
    # sanity: W-rule version must reproduce the C++ kept set
    same = np.array_equal(np.where(kept > -1e29)[0], np.sort(pop.assign(sel=lin).sort_values(['sel', 'row'], ascending=[False, True]).drop_duplicates('precursor_id').row.values))
    print(f'{RAW} s{s} (W reproduces C++ kept set: {same})\n   ' + '\n   '.join(out), flush=True)
