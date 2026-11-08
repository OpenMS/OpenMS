"""Offline replica of OpenDIACNN -stage2 (cross-fitted semi-supervised ridge LDA, pooled paired TDC) for raw 01 design ablations.
Target/decoy labels only. Reproduce the tool first, then add one feature group at a time."""
import numpy as np, pandas as pd, sqlite3, sys
D='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
LIB='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
BASE=['cnn_z','fragment_fraction','fragment_coelution','library_correlation','ms2_abs_ppm','ms1_abs_ppm','ms1_coelution','ms1_isotope']

def load():
    f=pd.read_csv(D+'/stage2-design/features.tsv',sep='\t')
    c=pd.read_csv(D+'/stage2-design/crossfit.tsv',sep='\t',usecols=['sample_id','group_id','precursor_best','q_value','in_rt_window'])
    s=pd.read_csv(D+'/pipeline/search.features.tsv',sep='\t')
    s=s.rename(columns={'candidate_row':'sample_id'}).drop(columns=['precursor_id','is_decoy'])
    f=f.merge(c,on='sample_id').merge(s,on='sample_id')
    con=sqlite3.connect(LIB); pr=pd.read_sql('select ID precursor_id, TRAML_ID, CHARGE charge, LIBRARY_RT lib_rt from PRECURSOR',con); con.close()
    pr['key']=pr.TRAML_ID.str.replace('^DECOY_','',regex=True)
    return f.merge(pr[['precursor_id','key','charge','lib_rt']],on='precursor_id')

def winners(df,score):
    o=np.lexsort((-df.sample_id.values, score))[::-1]  # highest score first
    sub=df.iloc[o]; keep=~sub.precursor_id.duplicated()
    return sub.index[keep.values]

def qvals(df,score_col):
    w=df.sort_values([score_col,'is_decoy'],ascending=[False,False]).drop_duplicates('key')
    T=np.cumsum(w.is_decoy.values==0); Dc=np.cumsum(w.is_decoy.values==1)
    q=np.minimum.accumulate(np.minimum(1,(Dc+1)/np.maximum(T,1))[::-1])[::-1]
    return w.assign(q=q)

def select(df,idx,score,fdr):
    w=df.loc[idx].assign(s=score[idx]).sort_values(['s','is_decoy'],ascending=[False,False]).drop_duplicates('key')
    T=np.cumsum(w.is_decoy.values==0); Dc=np.cumsum(w.is_decoy.values==1)
    ok=(Dc+1)/np.maximum(T,1)<=fdr
    k=np.where(ok)[0]; cut=k.max()+1 if len(k) else 0
    pos=w.index[:cut][w.is_decoy.values[:cut]==0]
    neg=df.loc[idx].index[df.loc[idx].is_decoy.values==1]
    return pos,neg

def stage2(df,feats,iters=3,f0=0.15,f1=0.05,ridge=1e-3):
    X=df[feats].to_numpy(float); out=np.zeros(len(df)); W=[]
    for fold in sorted(df.fold.unique()):
        tr=np.where(df.fold.values!=fold)[0]; te=np.where(df.fold.values==fold)[0]
        mu=X[tr].mean(0); sd=X[tr].std(0); sd[sd<=0]=1; Z=(X-mu)/sd
        score=np.full(len(df),np.nan); score[tr]=X[tr,0]
        trdf=df.iloc[tr]
        for it in range(iters):
            win=winners(trdf,score[tr]); pos,neg=select(df,win,score,f0 if it==0 else f1)
            pi=df.index.get_indexer(pos); ni=df.index.get_indexer(neg)
            mp=Z[pi].mean(0); mn=Z[ni].mean(0)
            S=((Z[pi]-mp).T@(Z[pi]-mp)+(Z[ni]-mn).T@(Z[ni]-mn))/(len(pi)+len(ni)-2)+ridge*np.eye(len(feats))
            w=np.linalg.solve(S,mp-mn); score[tr]=Z[tr]@w
        out[te]=Z[te]@w; W.append(w)
    return out,np.array(W)

def evaluate(df,score):
    d=df.assign(s2=score)
    win=d.loc[winners(d,score)]
    # per-fold decoy-winner median/MAD standardisation, as the tool does
    def st(g):
        x=g[g.is_decoy==1].s2; m=x.median(); s=1.4826*(x-m).abs().median(); return (g.s2-m)/s
    win['z']=win.groupby('fold',group_keys=False).apply(st)
    q=qvals(win,'z'); t=q[q.is_decoy==0]
    return {k:int((t.q<=k).sum()) for k in (0.01,0.02,0.05)}

if __name__=='__main__':
    df=load(); df=df[df.in_rt_window==1].reset_index(drop=True)
    tool=df[(df.precursor_best==1)&(df.is_decoy==0)]
    print('tool q<=0.01', int((pd.to_numeric(tool.q_value,errors='coerce')<=0.01).sum()))
    s,W=stage2(df,BASE); print('replica', evaluate(df,s)); print(np.round(W,2))
    df.to_parquet('raw01_rows.parquet')
