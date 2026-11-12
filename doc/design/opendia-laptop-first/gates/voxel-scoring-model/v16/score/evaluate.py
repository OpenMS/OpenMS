"""Evaluation-only joins, permitted only after both raw label-blind reconstruction gates pass."""
import json,sqlite3,time
from pathlib import Path
import numpy as np
import pandas as pd
from helpers import percentile,auc,membership_category,winners
OUT=Path(__file__).resolve().parent
V2=OUT.parents[2]
ROOT=V2.parents[1]
for raw in ['raw01','raw02']:
    assert (OUT/raw/'LABEL_BLIND_COMPLETE').exists()
    assert json.loads((OUT/raw/'verification.json').read_text())['pass']
(OUT/'EVALUATION_LABEL_JOIN.timestamp').write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
con=sqlite3.connect('file:'+str(ROOT/'astral-entrapment-20261107/lib/library.pqp')+'?mode=ro',uri=True)
ent=pd.read_sql("select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID",con).set_index('precursor_id').e
summary=[]
comparisons=[]
coverage=[]
support=[]
for raw,expected in [('raw01',(866,1413)),('raw02',(929,1362))]:
    P=pd.read_csv(V2/raw/'population.tsv',sep='\t')
    A=np.load(OUT/raw/'reconstructed.npz')
    rows=winners(P,A['linear'])
    W=P.iloc[rows].copy().reset_index(drop=True)
    W['row']=rows
    W['full']=A['scores'][rows].astype(float)
    W['linear']=A['linear'][rows].astype(float)
    W['nonlinear']=W.full-W.linear
    wt=np.load(V2/raw/'scores_D11FWT.npy')
    assert set(rows)==set(np.flatnonzero(wt>-1e20))
    W['WT']=wt[rows]
    stages=pd.read_csv(V2/'v3seq'/f'misses_D11FWT_{raw}.tsv',sep='\t').set_index('precursor_id').stage
    W['stage']=W.precursor_id.map(stages).fillna('not_reference')
    W['entrapment']=W.precursor_id.map(ent)
    assert W.entrapment.notna().all()
    W['group']=np.where(W.is_decoy==1,'DEC',np.where(W.entrapment==1,'ENT','H'))
    masks={'D':W.stage.eq('D right window, below cut'),'E':W.stage.eq('E accepted'),'allH':W.group.eq('H'),'ENT':W.group.eq('ENT'),'DEC':W.group.eq('DEC')}
    assert (int(masks['D'].sum()),int(masks['E'].sum()))==expected
    for f in range(3):
        ix=W.index[W.fold==f]
        dc=W[(W.fold==f)&masks['DEC']]
        for col in ['full','linear','nonlinear']:
            W.loc[ix,col+'_pct']=percentile(W.loc[ix,col],dc[col])
    W['rank_delta']=W.full_pct-W.linear_pct
    # Models for which each winner's family is eligible to train. Both exact window and any window of precursor.
    logs=[]
    for k in range(3):
        eligible=W[W.fold!=k]
        for scope in ['row','precursor']:
            for role in ['all_positive','training_positive','validation_positive']:
                bits=[]
                for r in range(3):
                    M=np.load(OUT/raw/f'membership_fold{k}_r{r}.npz')
                    ok=M['positive'].copy()
                    if role=='training_positive': ok &= ~M['validation']
                    if role=='validation_positive': ok &= M['validation']
                    selected=M['rows'][ok]
                    ids=selected if scope=='row' else P.precursor_id.values[selected]
                    query=eligible.row if scope=='row' else eligible.precursor_id
                    bits.append(np.isin(query,ids))
                B=np.column_stack(bits)
                L=eligible[['row','precursor_id','fold','stage','group','full_pct']].copy()
                L['model']=k
                L['scope']=scope
                L['role']=role
                L['r0'],L['r1'],L['r2']=B.T
                L['category']=membership_category(B)
                logs.append(L)
                for group,mask in masks.items():
                    ix=eligible.index[mask.loc[eligible.index]]
                    b=B[eligible.index.isin(ix)]
                    cat=membership_category(b)
                    coverage.append(dict(raw=raw,model=k,scope=scope,role=role,group=group,n=len(b),r0=int(b[:,0].sum()),r1=int(b[:,1].sum()),r2=int(b[:,2].sum()),never=int((cat=='never_admitted').sum()),lost=int((cat=='admitted_then_lost').sum()),retained=int((cat=='retained_r2').sum()),r1_to_r2_lost=int((b[:,1]&~b[:,2]).sum())))
    LL=pd.concat(logs,ignore_index=True)
    LL.to_csv(OUT/raw/'membership_joined.tsv.gz',sep='\t',index=False)
    # Collapse the two eligible models to distinct winners, preserving train versus validation distinction.
    for scope in ['row','precursor']:
        for role in ['all_positive','training_positive','validation_positive']:
            L=LL[(LL.scope==scope)&(LL.role==role)]
            pooled=L.groupby('row')[['r0','r1','r2']].max()
            assert len(pooled)==len(W)
            cats=pd.Series(membership_category(pooled.values),index=pooled.index)
            W[scope+'_'+role+'_coverage']=W.row.map(cats)
    # Held-out feature support is descriptive, thresholds learned without this fold.
    outside=np.zeros((len(P),17),np.uint8)
    for k in range(3):
        S=np.load(OUT/raw/f'support_fold{k}.npz')
        outside[S['rows']]=S['outside']
    W['outside_positive_support']=outside[rows].sum(1)
    for group,mask in masks.items():
        G=W[mask]
        summary.append(dict(raw=raw,group=group,n=len(G),full_pct_median=G.full_pct.median(),linear_pct_median=G.linear_pct.median(),nonlinear_pct_median=G.nonlinear_pct.median(),rank_delta_median=G.rank_delta.median(),down_any=int((G.rank_delta<0).sum()),down_01=int((G.rank_delta<=-.01).sum()),down_10=int((G.rank_delta<=-.1).sum()),up_01=int((G.rank_delta>=.01).sum()),cross99_down=int(((G.linear_pct>=.99)&(G.full_pct<.99)).sum()),cross99_up=int(((G.linear_pct<.99)&(G.full_pct>=.99)).sum()),support_outside_median=G.outside_positive_support.median()))
        for j in range(17): support.append(dict(raw=raw,group=group,feature_index=j,outside_fraction=float(outside[G.row,j].mean())))
        # Compare to all decoys and to full-score top-1% decoys within fold, not pooled fold scales.
        for region in ['all','top1pct','matched_full_strata']:
            for col in ['full','linear','nonlinear','outside_positive_support']:
                values=[]
                weights=[]
                for f in range(3):
                    gf=G[G.fold==f]
                    dc=W[(W.fold==f)&masks['DEC']]
                    if region=='top1pct': dc=dc[dc.full_pct>=.99]
                    edges=[0,.5,.9,.99,.999,1.00001] if region=='matched_full_strata' else [0,1.00001]
                    for lo,hi in zip(edges[:-1],edges[1:]):
                        a=gf[(gf.full_pct>=lo)&(gf.full_pct<hi)]
                        b=dc[(dc.full_pct>=lo)&(dc.full_pct<hi)] if region=='matched_full_strata' else dc
                        if len(a) and len(b): values.append(auc(a[col],b[col]));weights.append(len(a))
                comparisons.append(dict(raw=raw,group=group,control=region,measure=col,auc=float(np.average(values,weights=weights)) if weights else np.nan,covered=sum(weights),total=len(G)))
    W.to_csv(OUT/raw/'winner_diagnostics.tsv.gz',sep='\t',index=False)
    print(raw,'joined',len(W),'D/E',expected,flush=True)
pd.DataFrame(summary).to_csv(OUT/'score_summary.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(OUT/'training_coverage.tsv',sep='\t',index=False)
pd.DataFrame(comparisons).to_csv(OUT/'control_auc.tsv',sep='\t',index=False)
pd.DataFrame(support).to_csv(OUT/'feature_support.tsv',sep='\t',index=False)
