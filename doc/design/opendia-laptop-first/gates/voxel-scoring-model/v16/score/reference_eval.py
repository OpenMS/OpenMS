"""Use ONLY frozen reference scores. No failed reconstruction attribution."""
from pathlib import Path
import sqlite3
import numpy as np
import pandas as pd
from helpers import winners,percentile,auc
O=Path(__file__).resolve().parent
V2=O.parents[2]
C=sqlite3.connect('file:'+str(V2.parents[1]/'astral-entrapment-20261107/lib/library.pqp')+'?mode=ro',uri=True)
ent=pd.read_sql("select pp.PRECURSOR_ID precursor_id,min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID",C).set_index('precursor_id').e
summaries=[]
controls=[]
strata=[]
for raw,expected in [('raw01',(866,1413)),('raw02',(929,1362))]:
    P=pd.read_csv(V2/raw/'population.tsv',sep='\t')
    full=np.load(V2/raw/'scores_D11F.npy').astype(float)
    linear=np.load(V2/raw/'lin_D11F.npy').astype(float)
    wt=np.load(V2/raw/'scores_D11FWT.npy')
    rows=winners(P,linear)
    assert set(rows)==set(np.flatnonzero(wt>-1e20))
    W=P.iloc[rows].copy().reset_index(drop=True)
    W['row']=rows
    W['full'],W['linear'],W['nonlinear']=full[rows],linear[rows],full[rows]-linear[rows]
    W['WT']=wt[rows]
    st=pd.read_csv(V2/'v3seq'/f'misses_D11FWT_{raw}.tsv',sep='\t').set_index('precursor_id').stage
    W['stage']=W.precursor_id.map(st).fillna('not_reference')
    e=W.precursor_id.map(ent)
    assert e.notna().all()
    W['group']=np.where(W.is_decoy==1,'DEC',np.where(e==1,'ENT','H'))
    masks={'D':W.stage.eq('D right window, below cut'),'E':W.stage.eq('E accepted'),'allH':W.group.eq('H'),'ENT':W.group.eq('ENT'),'DEC':W.group.eq('DEC')}
    assert (int(masks['D'].sum()),int(masks['E'].sum()))==expected
    for f in range(3):
        ix=W.index[W.fold==f]
        dc=W[(W.fold==f)&masks['DEC']]
        for col in ['full','linear','nonlinear']: W.loc[ix,col+'_pct']=percentile(W.loc[ix,col],dc[col])
    W['rank_delta']=W.full_pct-W.linear_pct
    assert np.isfinite(W[['full','linear','nonlinear','full_pct','linear_pct','nonlinear_pct']]).all().all()
    for group,mask in masks.items():
        G=W[mask]
        summaries.append(dict(raw=raw,group=group,n=len(G),full_pct_median=G.full_pct.median(),linear_pct_median=G.linear_pct.median(),nonlinear_pct_median=G.nonlinear_pct.median(),rank_delta_median=G.rank_delta.median(),down_any=int((G.rank_delta<0).sum()),down_01=int((G.rank_delta<=-.01).sum()),down_10=int((G.rank_delta<=-.1).sum()),up_01=int((G.rank_delta>=.01).sum()),cross99_down=int(((G.linear_pct>=.99)&(G.full_pct<.99)).sum()),cross99_up=int(((G.linear_pct<.99)&(G.full_pct>=.99)).sum()),nonlinear_below_decoy_median=int((G.nonlinear_pct<.5).sum())))
        for lo,hi in zip([0,.5,.9,.99,.999],[.5,.9,.99,.999,1.00001]):
            g=G[(G.full_pct>=lo)&(G.full_pct<hi)]
            strata.append(dict(raw=raw,group=group,low=lo,high=hi,n=len(g),linear_pct_median=g.linear_pct.median(),nonlinear_pct_median=g.nonlinear_pct.median(),down01=int((g.rank_delta<=-.01).sum())))
        for region in ['all','top1pct','matched_full_strata']:
            for col in ['full','linear','nonlinear']:
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
                controls.append(dict(raw=raw,group=group,control=region,measure=col,auc=float(np.average(values,weights=weights)) if weights else np.nan,covered=sum(weights),total=len(G)))
    W.to_csv(O/f'reference_winners_{raw}.tsv.gz',sep='\t',index=False)
    print(raw,len(W),'D/E',expected,flush=True)
pd.DataFrame(summaries).to_csv(O/'reference_score_summary.tsv',sep='\t',index=False)
pd.DataFrame(controls).to_csv(O/'reference_control_auc.tsv',sep='\t',index=False)
pd.DataFrame(strata).to_csv(O/'reference_score_strata.tsv',sep='\t',index=False)
