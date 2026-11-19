"""Label-free claimant-strength rule (diagnostic): per accepted winner (q<=0.01), the largest number of its fragments within
20 ppm of a co-eluting accepted winner that scores HIGHER (any class; labels not used for the rule), within |dRT|<=10 s and
|dmz|<=12.5. Report the share with >=k such fragments by class (labels for evaluation only) and the effect of dropping them."""
import sqlite3, glob, bisect, numpy as np, pandas as pd
c = sqlite3.connect('file:hela_r5.predicted_decoys.pqp?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['pe'] = pq.key.map(pq[pq.d == 0].set_index('key').e)
k_ent = (pq[pq.d == 0].e == 1).sum() / (pq[pq.d == 0].e == 0).sum(); kk = 1 + 1 / k_ent
for f in sorted(glob.glob('runs/[DH]_cur_s0/*.features.tsv.gz')):
    F = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    F = F[F.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    A = F[F.q_value <= 0.01].join(pq[['d', 'mz', 'e', 'pe']], on='precursor_id').reset_index(drop=True)
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
    S['shs'] = best
    run = f.split('/')[1][0]
    cls = {'entrapment': S[(S.d == 0) & (S.e == 1)], 'decoy(entr pair)': S[(S.d == 1) & (S.pe == 1)], 'decoy(HeLa pair)': S[(S.d == 1) & (S.pe == 0)], 'HeLa': S[(S.d == 0) & (S.e == 0)]}
    print(run, 'accepted', len(S))
    for k in (2, 3, 4):
        print(f'  >= {k} fragments shared with a stronger co-eluting winner: ' + '  '.join(f'{n} {(v.shs >= k).mean():.3f} ({int((v.shs >= k).sum())})' for n, v in cls.items()))
    T = S[S.d == 0]; E = int(T.e.sum()); N = len(T)
    for k in (3, 4):
        keep = T[T.shs < k]; Ek = int(keep.e.sum())
        print(f'  drop >= {k}: targets {N} -> {len(keep)} ({len(keep)/N-1:+.2%}); lower bound {E/N:.4f} -> {Ek/len(keep):.4f}; combined {E*kk/N:.4f} -> {Ek*kk/len(keep):.4f} (q not recomputed)')
