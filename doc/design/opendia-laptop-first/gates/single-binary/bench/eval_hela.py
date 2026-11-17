"""HeLa diaPASEF 200 ng: OpenDIAWorkflow vs DIA-NN 2.2 vs OpenSWATH + PyProphet (SVM, XGBoost) on the same library.
No entrapment here, so every count is the tool's own q-value at <= 0.01 (decoy estimate), not a true FDR.
Precursors are matched by modified sequence (UniMod form) + charge."""
import sqlite3, re, numpy as np, pandas as pd
H = '/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
B = '/home/singjc/.jcode/scratch/e2e/hela'
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published'
lib = sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID pid, p.CHARGE z, p.DECOY decoy, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", lib)
pre['key'] = pre.ms + '/' + pre.z.astype(str)
key_of = dict(zip(pre.pid, pre.key))
od = pd.read_csv(f'{B}/opendia/precursors.tsv', sep='\t')
O = set(od[od.q_value <= 0.01].precursor_id.map(key_of))
dn = pd.read_parquet(f'{B}/diann_report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value'])
D = set((dn[dn['Q.Value'] <= 0.01]['Modified.Sequence'] + '/' + dn[dn['Q.Value'] <= 0.01]['Precursor.Charge'].astype(str)))
D_lib = D & set(pre[pre.decoy == 0].key)
def pyp(osw):
    c = sqlite3.connect(f'file:{osw}?mode=ro', uri=True)
    q = pd.read_sql("""select f.PRECURSOR_ID pid, min(s.QVALUE) q from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID
                       join PRECURSOR p on p.ID=f.PRECURSOR_ID where p.DECOY=0 and s.RANK=1 group by f.PRECURSOR_ID""", c)
    return set(q[q.q <= 0.01].pid.map(key_of))
S, X = pyp(f'{R}/svm.osw'), pyp(f'{R}/xgboost_t5.osw')
print(f'library targets {int((pre.decoy == 0).sum())}; DIA-NN keys found in library {len(D_lib)} of {len(D)}')
rows = [('OpenDIAWorkflow (this branch)', O), ('DIA-NN 2.2 (own decoys)', D), ('OpenSWATH + PyProphet SVM', S), ('OpenSWATH + PyProphet XGBoost', X)]
for name, s in rows:
    print(f'{name:32s} N {len(s):6d} | in DIA-NN {len(s & D) / max(len(s), 1):.3f} | in SVM {len(s & S) / max(len(s), 1):.3f} | of DIA-NN covered {len(s & D) / len(D):.3f}')
print(f'OpenDIA / DIA-NN {len(O) / len(D):.3f}; OpenDIA / SVM {len(O) / len(S):.3f}; OpenDIA / XGB {len(O) / len(X):.3f}')
print(f'OpenDIA only (not DIA-NN, not SVM): {len(O - D - S)}; DIA-NN and SVM but not OpenDIA: {len((D & S) - O)}')
