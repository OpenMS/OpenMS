"""PRECHECK_MS1 (diagnostic; labels only for grouping). Arm-E population, C++ scorer seed 0-4 scores.
MS1 evidence per window from ms1.npy [n, 4, 32] (M+0..M+3 summed within the MS1 ppm window, apex +- 15 s, 32 bins) and
xic.npy [n, 6, 32] (6 search fragments). Label-free features, computed identically for targets and decoys:
  m1_corr_xic   Pearson of the M+0 trace and the fragment-sum trace over the central +- 9.4 s (bins 6..26)
  m1_iso_cos    cosine of observed central areas M+0..M+3 against the theoretical envelope of the precursor's composition
  m1_iso_corr   mean Pearson of M+1 and M+2 traces with M+0 over the central bins
  m1_apex_dt    |apex bin of M+0 - apex bin of fragment sum| (bins)
  m1_log_m0     log1p central M+0 area
Why D11M (v4) raised FDP_c: unknown; candidates checked here:
  (a) decoy/target asymmetry: a decoy has the target's precursor m/z, so its MS1 trace is the TARGET's (same m/z, same RT
      window centre). MS1 features can only separate a decoy window from the target's real peak by RT agreement with the
      decoy's own fragment trace. Check: AUC of each feature for decoy winners vs entrapment winners vs null windows.
  (b) entrapment leakage: an entrapment has its own m/z (shuffled sequence, same composition -> SAME m/z and isotope
      envelope as its human parent!). Its MS1 trace is the parent's, so MS1 evidence of the parent's peak transfers to
      the co-eluting entrapment. Check: entrapment-in-band vs matched decoys AUC (specificity), and same-composition share.
Groups (score-matched, as PRECHECK_JOINT): R = DIA-NN-supported right-peak target winners rejected at today's cut;
D = 3 decoy winners nearest in score per R; E = entrapment winners in the same band; ED = their matched decoys.
Null windows: target and decoy windows with |library RT error| > 3 sigma of the calibration (no peptide expected there).
Usage: ms1_precheck.py"""
import json, re, sqlite3, sys, numpy as np, pandas as pd
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'; B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
M = '/home/singjc/.jcode/scratch/e2e/ms1'
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
from sklearn.metrics import roc_auc_score
AA = {'G': (2, 3, 1, 1, 0), 'A': (3, 5, 1, 1, 0), 'S': (3, 5, 1, 2, 0), 'P': (5, 7, 1, 1, 0), 'V': (5, 9, 1, 1, 0), 'T': (4, 7, 1, 2, 0),
      'C': (3, 5, 1, 1, 1), 'L': (6, 11, 1, 1, 0), 'I': (6, 11, 1, 1, 0), 'N': (4, 6, 2, 2, 0), 'D': (4, 5, 1, 3, 0), 'Q': (5, 8, 2, 2, 0),
      'K': (6, 12, 2, 1, 0), 'E': (5, 7, 1, 3, 0), 'M': (5, 9, 1, 1, 1), 'H': (6, 7, 3, 1, 0), 'F': (9, 9, 1, 1, 0), 'R': (6, 12, 4, 1, 0),
      'Y': (9, 9, 1, 2, 0), 'W': (11, 10, 2, 1, 0)}
MOD = {'UniMod:4': (2, 3, 1, 1, 0), 'UniMod:35': (0, 0, 0, 1, 0)}
ISO = [(0, [0.0107, 0.0]), (1, [0.000115, 0.0]), (2, [0.00364, 0.0]), (3, [0.00038, 0.00205]), (4, [0.0075, 0.0425])]
def envelope(comp, n=4):
    p = np.zeros(n); p[0] = 1.0
    for e, ab in ISO:
        k = int(comp[e])
        if k == 0: continue
        dist = np.zeros(n); dist[0] = 1 - sum(ab)
        for j, a in enumerate(ab, 1):
            if j < n: dist[j] = a
        for _ in range(k): p = np.convolve(p, dist)[:n]
    return p / p.sum()
def composition(ms):
    comp = np.array([0, 2, 0, 1, 0])
    for m in re.finditer(r'\((UniMod:\d+)\)', ms): comp += MOD.get(m.group(1), (0, 0, 0, 0, 0))
    for r in re.sub(r'\(.*?\)', '', ms): comp += AA[r]
    return comp
def corr(a, b):
    a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1) + 1e-12)
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,p.DECOY decoy,p.LIBRARY_RT lrt,pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2).drop_duplicates('precursor_id')
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key'); P = pre.set_index('precursor_id')
NAMES = ['m1_corr_xic', 'm1_iso_cos', 'm1_iso_corr', 'm1_apex_dt', 'm1_log_m0']
for RAW, REF in (('raw01', B + '/diann-pairlib-20261113/report_raw01.parquet'), ('raw02', B + '/diann-decoys-20261114/own/raw02/report.parquet')):
    pop = pd.read_csv(f'{G}/armE_{RAW}/population.tsv', sep='\t'); pop['row'] = np.arange(len(pop)); n = len(pop)
    env = {p: envelope(composition(P.ms[p])) for p in pop.precursor_id.unique()}
    E = np.stack([env[p] for p in pop.precursor_id.values])
    ms1 = np.load(f'{M}/{RAW}/ms1.npy', mmap_mode='r'); xic = np.load(f'{M}/{RAW}/xic.npy', mmap_mode='r')
    F = np.zeros((n, 5), np.float32)
    for b in range(0, n, 200000):
        m = np.asarray(ms1[b:b + 200000], np.float64); x = np.asarray(xic[b:b + 200000], np.float64).sum(1); e = E[b:b + 200000]
        mc, xc = m[:, :, 6:26], x[:, 6:26]
        has = mc[:, 0].max(1) > 0
        area = m[:, :, 12:20].sum(2); obs = area / np.maximum(area.sum(1, keepdims=True), 1e-12)
        F[b:b + 200000, 0] = np.where(has, corr(mc[:, 0], xc), 0)
        F[b:b + 200000, 1] = np.where(has, (obs * e).sum(1) / (np.linalg.norm(obs, axis=1) * np.linalg.norm(e, axis=1) + 1e-12), 0)
        F[b:b + 200000, 2] = np.where(has, 0.5 * (corr(mc[:, 0], mc[:, 1]) + corr(mc[:, 0], mc[:, 2])), 0)
        F[b:b + 200000, 3] = np.where(has, np.abs(mc[:, 0].argmax(1) - xc.argmax(1)), 20)
        F[b:b + 200000, 4] = np.log1p(area[:, 0])
    F = np.nan_to_num(F)
    np.save(f'{M}/{RAW}/ms1feat.npy', F)
    cal = json.load(open(f'{G}/armE_{RAW}/cal/cal.json'))['rt']
    kx, ky = np.array(cal['knots_library_rt']), np.array(cal['knots_run_rt_s'])
    lr = P.lrt.reindex(pop.precursor_id).values; mp = np.interp(lr, kx, ky)
    null = np.abs(pop.rt_apex.values - mp) > 3 * cal['sigma_s']
    ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
    ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
    ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid')
    print(f'== {RAW}: windows {n}; MS1 signal present {np.mean(F[:, 4] > 0):.3f}; null windows {null.sum()}')
    dec = pop.is_decoy.values == 1
    print('  null windows, AUC target > decoy (tell check; 0.5 = none): ' + ', '.join(f'{nm} {roc_auc_score(~dec[null], F[null, j]):.3f}' for j, nm in enumerate(NAMES)))
    for s in (0, 2):
        sc = np.load(f'{G}/cppscore/{RAW}_s{s}.npy').astype(float)
        w, best = tdc(pop, sc)
        tw = w[w.is_decoy == 0].sort_values('s', ascending=False); ee = tw.precursor_id.map(ent).values
        fdp = np.cumsum(ee) * (1 + 1 / R) / np.arange(1, len(tw) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1); cut = tw.s.iloc[n1 - 1]
        acc = set(tw.precursor_id.iloc[:n1]); W = w.set_index('precursor_id')
        right = [p for p in ref.index if p in W.index and abs(pop.rt_apex.values[int(W.row[p])] - ref.RT[p] * 60) <= 6]
        Rj = [p for p in right if p not in acc]; rR = W.loc[Rj, 'row'].astype(int).values; sR = W.loc[Rj, 's'].values
        dw = w[w.is_decoy == 1].sort_values('s'); ds = dw.s.values
        def matched(v):
            ix = np.clip(np.searchsorted(ds, v), 1, len(ds) - 2); return np.unique(np.r_[dw.row.values[ix - 1], dw.row.values[ix], dw.row.values[ix + 1]].astype(int))
        rD = matched(sR)
        Ew = tw[(tw.precursor_id.map(ent) == 1) & (tw.s < cut) & (tw.s >= np.median(sR))]; rE = Ew.row.values.astype(int); rED = matched(Ew.s.values)
        # accepted entrapments (the FDP_c driver) vs accepted decoy-free baseline
        Ea = tw.iloc[:n1]; rEa = Ea[Ea.precursor_id.map(ent) == 1].row.values.astype(int); rHa = Ea[Ea.precursor_id.map(ent) == 0].row.values.astype(int)
        print(f'  seed {s}: R {len(rR)} D {len(rD)} | E in band {len(rE)} ED {len(rED)} | accepted entrapments {len(rEa)}')
        for j, nm in enumerate(NAMES):
            a = roc_auc_score(np.r_[np.ones(len(rR)), np.zeros(len(rD))], np.r_[F[rR, j], F[rD, j]])
            b = roc_auc_score(np.r_[np.ones(len(rE)), np.zeros(len(rED))], np.r_[F[rE, j], F[rED, j]]) if len(rE) > 10 else np.nan
            ha = roc_auc_score(np.r_[np.ones(len(rHa)), np.zeros(len(rEa))], np.r_[F[rHa, j], F[rEa, j]]) if len(rEa) > 5 else np.nan
            print(f'    {nm:12s} AUC R>D {a:.3f} | E>ED {b:.3f} | accepted human > accepted entrapment {ha:.3f} | med R {np.median(F[rR, j]):.3f} D {np.median(F[rD, j]):.3f} E {np.median(F[rE, j]):.3f}')
