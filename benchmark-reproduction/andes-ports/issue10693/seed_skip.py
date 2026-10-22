import sys, numpy as np, pandas as pd, pyopenms as oms
exec(open('/home/user/bench/scripts/seed_check.py').read().split("g['seed'] =")[0].replace("link, pairs, amap, seedf, ffidf = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]",
 "link, pairs, amap, seedf, ffidf = '/home/user/runs/diag/link.tsv','/home/user/runs/diag/cand_pairs.tsv.gz',0,'/home/user/runs/diag/seeds/debug_seeds_fraction_1_1.featureXML','/home/user/runs/diag/seeds/debug_fraction_1_1.featureXML'"))
prots, peps = [], oms.PeptideIdentificationList() if hasattr(oms, 'PeptideIdentificationList') else []
oms.IdXMLFile().load('/home/user/data/ids_hye/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.idXML', prots, peps)
rows=[]
for p in peps:
    h=p.getHits()[0]; td=h.getMetaValue('target_decoy') if h.metaValueExists('target_decoy') else '?'
    q=h.getMetaValue('q-value_score') if h.metaValueExists('q-value_score') else (h.getMetaValue('q-value') if h.metaValueExists('q-value') else float('nan'))
    rows.append((p.getMZ(), p.getRT(), str(td), float(q), h.getSequence().toString()+'_'+str(h.getCharge())))
rows.sort(); P=np.array([(a,b) for a,b,_,_,_ in rows]); META=[(c,d,e) for _,_,c,d,e in rows]; print('PSMs', len(P))
C = 1.0033548378
def skipped(mz, rt, z):
    for k in range(6):
        m = mz - k / z * C; tol = m * 10e-6
        lo, hi = np.searchsorted(P[:, 0], m - tol), np.searchsorted(P[:, 0], m + tol, 'right')
        for i in range(lo, hi):
            if abs(P[i, 1] - rt) <= 30: return (k,) + META[i]
    return None
m_, r_, c_, _ = S
res = []
for row in g[g.status.isin(['no_feature', 'outside_rt'])].itertuples():
    lo, hi = np.searchsorted(m_, row.mz_t * (1 - 10e-6)), np.searchsorted(m_, row.mz_t * (1 + 10e-6), 'right')
    cand = [i for i in range(lo, hi) if c_[i] == row.charge and abs(r_[i] - row.rt_raw) <= row.gate]
    if not cand: continue
    sk = [skipped(m_[i], r_[i], c_[i]) for i in cand]
    if all(sk):
        k, td, q, key = sk[0]
        same = key == row.key
        res.append(('skipped', k, td, q <= 0.01, same))
x = pd.DataFrame(res, columns=['what', 'isotope_offset', 'target_decoy', 'q<=0.01', 'same_precursor'])
print(len(x), 'skipped seeds; blocking PSM:')
for c in ['isotope_offset', 'target_decoy', 'q<=0.01', 'same_precursor']: print((x[c].value_counts(normalize=True) * 100).round(1).to_string(), '\n')
