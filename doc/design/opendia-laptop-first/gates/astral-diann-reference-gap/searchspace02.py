"""Raw 02 (spent, descriptive): pooled paired q over stage2_score restricted to DIA-NN's search space (charge 1-4, m/z 300-1800)."""
import pandas as pd, numpy as np, sqlite3, sys
B='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
sys.path.insert(0,B+'/astral-entrapment-raw02-20261108'); import evaluate
LIB=B+'/astral-entrapment-20261107/lib/library.pqp'
lab=evaluate.library_labels(LIB); r=148983/151240
con=sqlite3.connect(LIB); pr=pd.read_sql('select ID precursor_id, TRAML_ID, CHARGE charge, PRECURSOR_MZ mz from PRECURSOR',con); con.close()
c=pd.read_csv(B+'/astral-raw02-stage2-check-20261108/crossfit/crossfit.tsv',sep='\t',usecols=['precursor_id','is_decoy','fold','precursor_best','q_value','stage2_score'])
b=c[c.precursor_best==1].merge(pr,on='precursor_id'); b['key']=b.TRAML_ID.str.replace('^DECOY_','',regex=True)
# Per-fold decoy-winner median/MAD standardisation, as the tool does before pooling.
def st(g):
    d=g[g.is_decoy==1].stage2_score; m=d.median(); s=1.4826*(d-m).abs().median(); return (g.stage2_score-m)/s
b['stage2_score']=b.groupby('fold',group_keys=False).apply(st)
ref=set(pd.read_parquet('ref_exact.parquet').query('Run.str.endswith("_02")',engine='python').precursor_id)
def run(df,name):
    df=df.sort_values(['stage2_score','is_decoy'],ascending=[False,False]).drop_duplicates('key').sort_values('stage2_score',ascending=False)
    T=np.cumsum(df.is_decoy.values==0); D=np.cumsum(df.is_decoy.values==1)
    q=np.minimum.accumulate(np.minimum(1,(D+1)/np.maximum(T,1))[::-1])[::-1]; df=df.assign(q2=q)
    for qq in (0.01,0.05):
        a=df[(df.is_decoy==0)&(df.q2<=qq)]; ne=int(a.precursor_id.map(lab).sum())
        print(f'{name:28s} pairs={len(df):7d} q<={qq}: accepted={len(a)} N_e={ne} FDP_c={ne*(1+1/r)/max(len(a),1):.4f} DIA-NN-confirmed={int(a.precursor_id.isin(ref).sum())}')
run(b,'all (reproduces tool q)')
run(b[b.charge.between(1,4)&b.mz.between(300,1800)],'charge 1-4, m/z 300-1800')
run(b[b.charge.between(2,3)&b.mz.between(300,1800)],'charge 2-3, m/z 300-1800')
print('library target precursors by charge', pr[~pr.TRAML_ID.str.startswith('DECOY_')].charge.value_counts().sort_index().to_dict())
