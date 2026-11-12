"""Post-verification coverage and selector audit. No model fitting or label-selected inputs."""
import json,time
from pathlib import Path
import numpy as np
import pandas as pd
from helpers import winners,membership_category
from mismatch_helpers import compare_rows
O=Path(__file__).resolve().parent
V2=O.parents[3]
for raw in ['raw01','raw02']:
    assert (O/raw/'LABEL_BLIND_COMPLETE').exists()
    assert json.loads((O/raw/'verification.json').read_text())['pass']
(O/'EVALUATION_JOIN.timestamp').write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
coverage=[]
category_summary=[]
selector_summary=[]
member_mismatch=[]
final_summary=[]
feature_summary=[]
for raw in ['raw01','raw02']:
    P=pd.read_csv(V2/raw/'population.tsv',sep='\t')
    W=pd.read_csv(O.parent/f'reference_winners_{raw}.tsv.gz',sep='\t').sort_values('row').reset_index(drop=True)
    assert P.iloc[W.row].precursor_id.to_list()==W.precursor_id.to_list()
    wr=W.row.to_numpy()
    n=len(W)
    masks={'D':W.stage.eq('D right window, below cut').to_numpy(),'E':W.stage.eq('E accepted').to_numpy(),'allH':W.group.eq('H').to_numpy(),'ENT':W.group.eq('ENT').to_numpy(),'DEC':W.group.eq('DEC').to_numpy()}
    # Indexed [winner, model, round, role], eligibility retained separately to avoid confusing held-out with absent.
    roles=['all_positive','training_positive','validation_positive']
    exact=np.zeros((n,3,3,3),bool)
    precursor=np.zeros_like(exact)
    far=np.zeros_like(exact)
    eligible=np.column_stack([W.fold.ne(k).to_numpy() for k in range(3)])
    assert (eligible.sum(1)==2).all()
    pids=P.precursor_id.to_numpy()
    rt=P.rt_apex.to_numpy()
    fullref=np.load(V2/raw/'scores_D11F.npy')
    fr=winners(P,fullref)
    frows,diff,dist=compare_rows(P,wr,fr)
    W['final_full_row']=frows
    W['final_full_linear_row_diff']=diff
    W['final_full_linear_rt_distance']=dist
    for group,mask in masks.items():
        final_summary.append(dict(raw=raw,group=group,n=int(mask.sum()),different=int(diff[mask].sum()),far_gt10=int((dist[mask]>10).sum()),median_distance=float(np.median(dist[mask]))))
    membership_details=[]
    selector_details=[]
    for k in range(3):
        for r in range(3):
            M=np.load(O/raw/f'membership_fold{k}_r{r}.npz')
            mr=M['rows']
            assert np.all(P.fold.to_numpy()[mr]!=k)
            for j,role in enumerate(roles):
                ok=M['positive'].copy()
                if role=='training_positive': ok &= ~M['validation']
                if role=='validation_positive': ok &= M['validation']
                pos=mr[ok]
                exact[:,k,r,j]=np.isin(wr,pos)&eligible[:,k]
                precursor[:,k,r,j]=np.isin(W.precursor_id,pids[pos])&eligible[:,k]
                if len(pos):
                    # Actual admitted row vs the eventual held-out LINEAR-selected W row.
                    wi,md,dt=compare_rows(P,pos,wr)
                    far_pids=np.unique(pids[pos[dt>10]])
                    far[:,k,r,j]=np.isin(W.precursor_id,far_pids)&eligible[:,k]
                    L=pd.DataFrame({'row':pos,'winner_row':wi,'model':k,'round':r,'role':role,'different':md,'rt_distance':dt})
                    L=L.merge(W[['row','stage','group','full_pct']],left_on='winner_row',right_on='row',suffixes=('','_winner'),validate='many_to_one')
                    membership_details.append(L)
                    for group,mask in masks.items():
                        use=np.isin(wi,wr[mask])
                        member_mismatch.append(dict(raw=raw,model=k,round=r,role=role,group=group,admitted_rows=int(use.sum()),different=int(md[use].sum()),far_gt10=int((dt[use]>10).sum())))
                for group,mask in masks.items():
                    use=mask&eligible[:,k]
                    coverage.append(dict(raw=raw,model=k,round=r,role=role,group=group,eligible=int(use.sum()),exact=int(exact[use,k,r,j].sum()),precursor=int(precursor[use,k,r,j].sum()),admitted_other_apex_gt10=int(far[use,k,r,j].sum())))
            if r:
                S=np.load(O/raw/f'selectors_fold{k}_r{r}.npz')
                fl=S['full_rows']
                lr,sd,ds=compare_rows(P,fl,S['linear_rows'])
                assert np.all(P.fold.to_numpy()[fl]!=k)
                infer,hd,hdist=compare_rows(P,fl,wr)
                selected=np.isin(fl,mr[M['positive']])
                selected_train=np.isin(fl,mr[M['positive']&~M['validation']])
                selected_val=np.isin(fl,mr[M['positive']&M['validation']])
                T=pd.DataFrame({'full_row':fl,'linear_row':lr,'winner_row':infer,'model':k,'round':r,'same_model_different':sd,'same_model_rt_distance':ds,'heldout_different':hd,'heldout_rt_distance':hdist,'selected_positive':selected,'selected_training_positive':selected_train,'selected_validation_positive':selected_val})
                # Save label-free row identities plus evaluation group only after baseline verification.
                T=T.merge(W[['row','stage','group','full_pct']],left_on='winner_row',right_on='row',validate='one_to_one')
                selector_details.append(T)
                for group,mask in masks.items():
                    g=np.isin(infer,wr[mask])
                    for subset,sub in [('all_eligible',np.ones(len(fl),bool)),('selected_positive',selected),('selected_training_positive',selected_train),('selected_validation_positive',selected_val)]:
                        use=g&sub
                        selector_summary.append(dict(raw=raw,model=k,round=r,group=group,subset=subset,n=int(use.sum()),same_model_different=int(sd[use].sum()),same_model_far_gt10=int((ds[use]>10).sum()),heldout_different=int(hd[use].sum()),heldout_far_gt10=int((hdist[use]>10).sum())))
    np.savez_compressed(O/raw/'coverage_arrays.npz',winner_rows=wr,eligible=eligible,exact=exact,precursor=precursor,far=far,roles=np.array(roles))
    for scope,B in [('exact',exact),('precursor',precursor)]:
        for j,role in enumerate(roles):
            bits=B[:,:,:,j].any(axis=1)
            cat=membership_category(bits)
            W[scope+'_'+role+'_category']=cat
            for r in range(3): W[f'{scope}_{role}_r{r}']=bits[:,r]
            W[scope+'_'+role+'_models_retained']=B[:,:,2,j].sum(1)
            for group,mask in masks.items():
                for c in ['never_admitted','admitted_then_lost','retained_r2']:
                    use=mask&(cat==c)
                    category_summary.append(dict(raw=raw,scope=scope,role=role,group=group,category=c,n=int(use.sum()),group_total=int(mask.sum()),heldout_below_decoy50=int((W.full_pct.to_numpy()[use]<.5).sum()),heldout_below_decoy90=int((W.full_pct.to_numpy()[use]<.9).sum()),heldout_above_decoy99=int((W.full_pct.to_numpy()[use]>=.99).sum()),full_pct_median=float(np.median(W.full_pct.to_numpy()[use])) if use.any() else np.nan,linear_pct_median=float(np.median(W.linear_pct.to_numpy()[use])) if use.any() else np.nan))
    W['ever_precursor_never_exact_positive']=precursor[:,:,:,0].any(axis=(1,2))&~exact[:,:,:,0].any(axis=(1,2))
    W['ever_far_positive']=far[:,:,:,0].any(axis=(1,2))
    W['r2_far_positive']=far[:,:,2,0].any(axis=1)
    # Label-blind final-training support bounds evaluated ONLY after exact model reproduction.
    outside=np.zeros((len(P),17),np.uint8)
    for k in range(3):
        S=np.load(O/raw/f'support_fold{k}.npz')
        outside[S['rows']]=S['outside']
    W['features_outside_training_positive_1_99']=outside[wr].sum(1)
    for group,mask in masks.items():
        for f in range(17): feature_summary.append(dict(raw=raw,group=group,feature_index=f,n=int(mask.sum()),outside_fraction=float(outside[wr[mask],f].mean())))
    W.to_csv(O/raw/'coverage_winners.tsv.gz',sep='\t',index=False)
    pd.concat(membership_details,ignore_index=True).to_csv(O/raw/'admitted_row_mismatches.tsv.gz',sep='\t',index=False)
    pd.concat(selector_details,ignore_index=True).to_csv(O/raw/'selector_rows.tsv.gz',sep='\t',index=False)
    print(raw,'coverage and selector diagnostics complete',flush=True)
for name,rows in [('coverage_by_model_round',coverage),('coverage_categories',category_summary),('selector_mismatch',selector_summary),('admitted_row_mismatch',member_mismatch),('final_selector_mismatch',final_summary),('feature_support',feature_summary)]:
    pd.DataFrame(rows).to_csv(O/(name+'.tsv'),sep='\t',index=False)
