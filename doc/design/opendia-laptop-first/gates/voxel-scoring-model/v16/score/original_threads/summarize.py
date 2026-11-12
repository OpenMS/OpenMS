"""Summarize already validated coverage outputs; no model or feature choice."""
from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parent
C=pd.read_csv(O/'coverage_categories.tsv',sep='\t')
print('COVERAGE ALL POSITIVE, DISTINCT WINNERS')
print(C[(C.role=='all_positive')&C.group.isin(['D','E','allH','ENT','DEC'])][['raw','scope','group','category','n','heldout_below_decoy90','full_pct_median']].to_string(index=False))
print('\nTRAINING ONLY DISTINCT WINNERS D/E')
print(C[(C.role=='training_positive')&C.group.isin(['D','E'])][['raw','scope','group','category','n','heldout_below_decoy90']].to_string(index=False))
S=pd.read_csv(O/'selector_mismatch.tsv',sep='\t')
A=S.groupby(['raw','round','group','subset'],as_index=False)[['n','same_model_different','same_model_far_gt10','heldout_different','heldout_far_gt10']].sum()
A.to_csv(O/'selector_mismatch_pooled.tsv',sep='\t',index=False)
print('\nR2 SAME MODEL FULL->LINEAR AND TRAINING FULL->HELDOUT LINEAR')
print(A[(A['round']==2)&A.subset.isin(['all_eligible','selected_positive'])].to_string(index=False))
print('\nFINAL HELDOUT FULL->LINEAR')
print(pd.read_csv(O/'final_selector_mismatch.tsv',sep='\t').to_string(index=False))
rows=[]
for raw in ['raw01','raw02']:
    W=pd.read_csv(O/raw/'coverage_winners.tsv.gz',sep='\t')
    masks={'D':W.stage.eq('D right window, below cut'),'E':W.stage.eq('E accepted'),'allH':W.group.eq('H'),'ENT':W.group.eq('ENT'),'DEC':W.group.eq('DEC')}
    for g,m in masks.items():
        G=W[m]
        rows.append(dict(raw=raw,group=g,n=len(G),ever_precursor_never_exact=int(G.ever_precursor_never_exact_positive.sum()),ever_far_positive=int(G.ever_far_positive.sum()),r2_far_positive=int(G.r2_far_positive.sum()),median_features_outside=G.features_outside_training_positive_1_99.median()))
    print('\n',raw,'D HELDOUT SCORE <decoy90 BY EXACT COVERAGE')
    print(pd.crosstab(W.loc[masks['D'],'exact_all_positive_category'],W.loc[masks['D'],'full_pct']<.9).to_string())
    print(raw,'D R2 retained model counts',W.loc[masks['D'],'exact_all_positive_models_retained'].value_counts().to_dict())
pd.DataFrame(rows).to_csv(O/'distinct_winner_mismatch.tsv',sep='\t',index=False)
print('\nDISTINCT WINNER OTHER-APEX ADMISSION')
print(pd.DataFrame(rows).to_string(index=False))
