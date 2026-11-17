"""PRECHECK_FRAG (diagnostic; labels only in the grouping). Would choosing which 6 of a precursor's 12 predicted fragments
to search, from what this run actually observes, raise the evidence of the weak right-peak IDs more than that of decoys?
Inputs: arm-E raw01/raw02 xic12 [n,12,32] (12 predicted fragments, library order), lib12 [n,12], orig12 [n,12] (1 = one
of the 6 searched fragments), annotations from library12_pop.pqp; C++ scorer seed 0 scores.
1. Observed rate: share of windows whose central area (bins 12..19) of fragment slot k is > 0, by library rank, for
   A (accepted targets), R (DIA-NN-supported right-peak rejected), D (score-matched decoys), E/ED (entrapment band).
2. Run-specific ion-class rates, label-free and cross-fitted: for the calibration-anchor-like set = target winners at the
   current q <= 0.01 of the OTHER folds, the rate at which each ion class (type b/y x charge x ordinal bin x rank) is
   observed. Re-select per precursor the 6 fragments with the highest class rate x predicted intensity.
   Compare 'support' (number of the 6 with central signal) and the co-elution of the chosen 6 under original vs
   re-selected choice, for R vs D and E vs ED (specificity). Gain is only useful if R rises more than D and E does not.
Usage: frag_precheck.py"""
import sqlite3, sys, numpy as np, pandas as pd
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'; B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
X02 = '/home/singjc/.jcode/scratch/e2e/joint/raw02'
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
from sklearn.metrics import roc_auc_score
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,p.DECOY decoy,pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2)
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key')
c12 = sqlite3.connect(f'file:{G}/armE_raw01/library12_pop.pqp?mode=ro', uri=True)
t12 = pd.read_sql("select m.PRECURSOR_ID pid, t.ANNOTATION ann, t.TYPE ty, t.ORDINAL ord, t.CHARGE ch from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID order by m.PRECURSOR_ID, t.ID", c12)
t12['rank'] = t12.groupby('pid').cumcount(); t12 = t12[t12['rank'] < 12]
# Ion class: type, charge, ordinal bin (1-2, 3-5, 6-9, 10+)
t12['cls'] = t12.ty.fillna('?') + t12.ch.fillna(0).astype(int).astype(str) + '_' + pd.cut(t12.ord.fillna(0), [-1, 2, 5, 9, 999], labels=['o1', 'o3', 'o6', 'o10']).astype(str)
cls_codes, cls_names = pd.factorize(t12.cls)
t12['c'] = cls_codes
NC = len(cls_names)
for RAW, REF in (('raw01', B + '/diann-pairlib-20261113/report_raw01.parquet'), ('raw02', B + '/diann-decoys-20261114/own/raw02/report.parquet')):
    D = f'{G}/armE_{RAW}'; XD = D if RAW == 'raw01' else X02
    pop = pd.read_csv(f'{D}/population.tsv', sep='\t'); pop['row'] = np.arange(len(pop)); n = len(pop)
    X = np.load(f'{XD}/xic12.npy', mmap_mode='r'); L = np.load(f'{XD}/lib12.npy'); O = np.load(f'{D}/orig12.npy') if RAW == 'raw01' else None
    if O is None:
        # raw02 orig12 is not local: derive from the 6-fragment search library annotation set
        c6 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
        t6 = pd.read_sql("select m.PRECURSOR_ID pid, t.ANNOTATION ann from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID", c6)
        orig = set(zip(t6.pid, t6.ann))
        t12['orig'] = [(p, a) in orig for p, a in zip(t12.pid, t12.ann)]
        Om = t12.pivot_table(index='pid', columns='rank', values='orig', aggfunc='first').reindex(columns=range(12)).fillna(False).astype(np.uint8)
        O = Om.reindex(pop.precursor_id).fillna(0).values.astype(np.uint8)
    C = t12.pivot_table(index='pid', columns='rank', values='c', aggfunc='first').reindex(columns=range(12))
    Cw = C.reindex(pop.precursor_id).fillna(-1).values.astype(int)
    A = np.zeros((n, 12), np.float32)
    for b in range(0, n, 200000): A[b:b + 200000] = np.asarray(X[b:b + 200000], np.float32)[:, :, 12:20].sum(2)
    obs = A > 0
    sc = np.load(f'{G}/cppscore/{RAW}_s0.npy').astype(float)
    w, best = tdc(pop, sc); bq = best.set_index('precursor_id').q
    tw = w[w.is_decoy == 0].sort_values('s', ascending=False); e = tw.precursor_id.map(ent).values
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(tw) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1); cut = tw.s.iloc[n1 - 1]
    acc = set(tw.precursor_id.iloc[:n1]); W = w.set_index('precursor_id')
    ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
    ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
    ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid')
    right = [p for p in ref.index if p in W.index and abs(pop.rt_apex.values[int(W.row[p])] - ref.RT[p] * 60) <= 6]
    rR = W.loc[[p for p in right if p not in acc], 'row'].astype(int).values; sR = W.loc[[p for p in right if p not in acc], 's'].values
    rA = W.loc[[p for p in right if p in acc], 'row'].astype(int).values
    dw = w[w.is_decoy == 1].sort_values('s'); ds = dw.s.values
    def matched(v):
        ix = np.clip(np.searchsorted(ds, v), 1, len(ds) - 2); return np.unique(np.r_[dw.row.values[ix - 1], dw.row.values[ix], dw.row.values[ix + 1]].astype(int))
    rD = matched(sR)
    Ew = tw[(tw.precursor_id.map(ent) == 1) & (tw.s < cut) & (tw.s >= np.median(sR))]; rE = Ew.row.values.astype(int); rED = matched(Ew.s.values)
    print(f'== {RAW}: A {len(rA)} R {len(rR)} D {len(rD)} E {len(rE)} ED {len(rED)}')
    print('  observed share by library rank (1..12):')
    for nm, rows in (('A', rA), ('R', rR), ('D', rD), ('E', rE)):
        print(f'    {nm}: ' + ' '.join(f'{v:.2f}' for v in obs[rows].mean(0)) + f' | searched 6: {(obs[rows] & (O[rows] == 1)).sum(1).mean():.2f} of 6')
    # Cross-fitted class rates from the run's own confident targets (other folds), label-free.
    fold = pop.fold.values
    conf = w[(w.is_decoy == 0)].copy(); conf['q'] = conf.precursor_id.map(bq).fillna(1); conf = conf[conf.q <= 0.01]
    newsel = np.zeros((n, 12), bool)
    for k in range(3):
        rows = conf[conf.fold != k].row.values.astype(int)
        cc = Cw[rows].ravel(); oo = obs[rows].ravel(); m = cc >= 0
        num = np.bincount(cc[m], oo[m], NC); den = np.bincount(cc[m], minlength=NC)
        rate = (num + 1) / (den + 2)
        te = fold == k
        sc_ = np.where(Cw[te] >= 0, rate[np.maximum(Cw[te], 0)] * L[te], -1)
        top = np.argsort(-sc_, 1)[:, :6]
        sel = np.zeros_like(sc_, bool); np.put_along_axis(sel, top, True, 1); sel &= sc_ > 0
        newsel[te] = sel
    def support(rows, S): return (obs[rows] & S[rows]).sum(1)
    Sorig = O == 1
    for nm, a, b in (('R vs D', rR, rD), ('E vs ED', rE, rED), ('A vs D', rA, rD)):
        so_a, so_b = support(a, Sorig), support(b, Sorig); sn_a, sn_b = support(a, newsel), support(b, newsel)
        y = np.r_[np.ones(len(a)), np.zeros(len(b))]
        print(f'  {nm}: support orig {so_a.mean():.2f} vs {so_b.mean():.2f} (AUC {roc_auc_score(y, np.r_[so_a, so_b]):.3f}) | re-selected {sn_a.mean():.2f} vs {sn_b.mean():.2f} (AUC {roc_auc_score(y, np.r_[sn_a, sn_b]):.3f})')
    ch = (newsel != Sorig).any(1)
    print(f'  re-selection changes the fragment set of {ch.mean():.3f} of windows; mean fragments swapped {((newsel & ~Sorig).sum(1)[ch]).mean():.2f}')
