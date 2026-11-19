"""Fair simulation (diagnostic) of a label-free shared-signal rule. Among precursor winners with tool q<=0.05 (targets AND
decoys), drop a winner if >=k of its fragments lie within 20 ppm of a co-eluting (|dRT|<=10 s, |dmz|<=12.5) STRONGER winner
(any kind). Then recompute q on the remaining winners by plain TDC, (D+1)/T in score order, and for comparison on all winners
the same way (baseline). Report N, lower and combined entrapment bounds at the recomputed q<=0.01. The rule uses no labels;
the bounds use labels for evaluation only. Caveat: plain TDC rather than the tool's paired, fold-standardised TDC."""
import sqlite3, glob, bisect, numpy as np, pandas as pd
c = sqlite3.connect('file:hela_r5.predicted_decoys.pqp?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
tg = pq[pq.d == 0]; r = (tg.e == 1).sum() / (tg.e == 0).sum(); kk = 1 + 1 / r
def qplain(df):
    df = df.sort_values('score', ascending=False); d = df.is_decoy.cumsum().values; t = (1 - df.is_decoy).cumsum().values
    q = (d + 1) / np.maximum(t, 1); q = np.minimum.accumulate(q[::-1])[::-1]; return df.assign(q2=q)
def bounds(df):
    a = df[(df.q2 <= 0.01) & (df.is_decoy == 0)]; E = int(a.e.sum()); N = len(a); return N, E / max(N, 1), E * kk / max(N, 1)
for f in sorted(glob.glob('runs/[DH]_cur_s0/*.features.tsv.gz')):
    F = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    F = F[F.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    A = F[F.q_value <= 0.05].join(pq[['mz', 'e']], on='precursor_id').reset_index(drop=True)
    fr = pd.read_sql(f"select tm.PRECURSOR_ID p, tr.PRODUCT_MZ m from TRANSITION_PRECURSOR_MAPPING tm join TRANSITION tr on tr.ID=tm.TRANSITION_ID "
                     f"where tm.PRECURSOR_ID in ({','.join(map(str, A.precursor_id))})", c)
    frag = fr.groupby('p').m.apply(lambda s: np.sort(s.values)).to_dict()
    S = A.sort_values('rt_apex').reset_index(drop=True); rt = S.rt_apex.values; mz = S.mz.values; sc = S.score.values; pid = S.precursor_id.values
    best = np.zeros(len(S), int)
    for i in range(len(S)):
        a = frag.get(pid[i])
        if a is None: continue
        lo = bisect.bisect_left(rt, rt[i] - 10); hi = bisect.bisect_right(rt, rt[i] + 10); b0 = 0
        for j in range(lo, hi):
            if j == i or sc[j] <= sc[i] or abs(mz[j] - mz[i]) > 12.5: continue
            b = frag.get(pid[j])
            if b is None: continue
            idx = np.searchsorted(b, a); n = 0
            for v, q in zip(a, idx):
                if any(0 <= t < len(b) and abs(b[t] - v) <= v * 20e-6 for t in (q - 1, q)): n += 1
            b0 = max(b0, n)
        best[i] = b0
    S['shs'] = best; S['e'] = S.e.fillna(0)
    run = f.split('/')[1][0]
    base = qplain(S); N0, l0, c0 = bounds(base)
    tool = S[(S.q_value <= 0.01) & (S.is_decoy == 0)]; Nt = len(tool); Et = int(tool.e.sum())
    print(f'{run}: tool q<=0.01 N {Nt} lower {Et/Nt:.4f} combined {Et*kk/Nt:.4f} | plain TDC baseline N {N0} lower {l0:.4f} combined {c0:.4f}')
    for k in (3, 4):
        dropped = S[S.shs >= k]
        rule = qplain(S[S.shs < k]); N1, l1, c1 = bounds(rule)
        print(f'   rule k={k}: dropped {len(dropped)} (decoys {int(dropped.is_decoy.sum())}); N {N1} ({N1/N0-1:+.2%}) lower {l1:.4f} combined {c1:.4f}')
