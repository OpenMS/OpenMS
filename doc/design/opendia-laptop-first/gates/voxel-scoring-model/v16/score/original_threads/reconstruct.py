"""Exact existing D11F recipe instrumentation. Never calls run or writes reference arrays."""
import os,sys,inspect,json,hashlib,time,importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from helpers import winners
OUT=Path(__file__).resolve().parent
V4=OUT.parents[2]
RAW=sys.argv[1]
os.environ['TRAIN_THREADS']='10'
sys.argv=['train_v4.py',RAW,'D11F']
spec=importlib.util.spec_from_file_location('original_train',V4/'train_v4.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
D=OUT/RAW
D.mkdir(exist_ok=True)
assert m.pop.sample_id.is_unique
assert (m.pop.groupby('pair_key').fold.nunique()==1).all()
Hx=np.c_[m.H,m.F12]
assert np.isfinite(Hx).all()
scores=np.zeros(len(m.pop),np.float32)
linear=np.zeros_like(scores)
roundlog=[]
for k in range(3):
    tr_all=np.where(m.fold!=k)[0]
    mu,sd=Hx[tr_all].mean(0),Hx[tr_all].std(0)
    sd[sd==0]=1
    Z=(Hx-mu)/sd
    t=m.train_rows[(m.train_rows.fold!=k)&m.train_rows.group.isin(['anchor','anchor_decoy','decoy'])]
    rows=m.row_of[t.sample_id].values
    y=t.label.values.astype(np.float32)
    valid=np.array([m.fnv('v'+str(s))%100<15 for s in t.family.values])
    for r in range(3):
        if r:
            win,q,dec,wtab=m.pooled_q(tr_all,sc(tr_all))
            # Label-blind snapshots of the actual full-score selector feeding this round.
            # Also compare the linear selector of that exact same prior-round model.
            previous=inspect.getclosurevars(sc).nonlocals
            prior_linear=(m.embed(previous['net'],tr_all,Z,False)[:,:Z.shape[1]]@previous['w'][:Z.shape[1]]).numpy()
            linear_winners=winners(m.pop.iloc[tr_all].reset_index(drop=True),prior_linear)
            np.savez_compressed(D/f'selectors_fold{k}_r{r}.npz',full_rows=wtab.r.values,linear_rows=tr_all[linear_winners],competition_rows=win,q=q,decoy=dec)
            pos=win[(dec==0)&(q<=(.15 if r==1 else .05))]
            negs=wtab[wtab.dec==1].r.values
            rows=np.concatenate([pos,negs])
            y=np.r_[np.ones(len(pos)),np.zeros(len(negs))].astype(np.float32)
            valid=np.array([m.fnv('v'+str(s))%100<15 for s in m.pair_key[rows]])
        assert np.all(m.fold[rows]!=k)
        np.savez_compressed(D/f'membership_fold{k}_r{r}.npz',rows=rows,sample_id=m.pop.sample_id.values[rows],positive=y==1,validation=valid)
        seed=100+k if r==0 else 200+10*k+r
        sc,ep,vo=m.fit(rows,y,Z,False,True,seed,valid,100,0.)
        cv=inspect.getclosurevars(sc).nonlocals
        torch.save({'state_dict':cv['net'].state_dict(),'readout':cv['w'],'mu':mu,'sd':sd,'seed':seed,'fold':k,'round':r,'epochs':ep,'validation_objective':vo},D/f'checkpoint_fold{k}_r{r}.pt')
        roundlog.append(dict(fold=k,round=r,seed=seed,epochs=ep,objective=vo,train_positive=int(((y==1)&~valid).sum()),validation_positive=int(((y==1)&valid).sum()),negative=int((y==0).sum())))
        print(RAW,roundlog[-1],flush=True)
    te=np.where(m.fold==k)[0]
    scores[te]=sc(te)
    linear[te]=(m.embed(cv['net'],te,Z,False)[:,:Z.shape[1]]@cv['w'][:Z.shape[1]]).numpy()
    # Label-blind feature support against final training positives only, using held-out rows.
    tp=rows[(y==1)&~valid]
    lo,hi=np.quantile(Hx[tp],[.01,.99],axis=0)
    np.savez_compressed(D/f'support_fold{k}.npz',rows=te,outside=((Hx[te]<lo)|(Hx[te]>hi)).astype(np.uint8),lower=lo,upper=hi)
np.savez_compressed(D/'reconstructed.npz',scores=scores,linear=linear)
pd.DataFrame(roundlog).to_csv(D/'rounds.tsv',sep='\t',index=False)
ref=np.load(Path(m.O)/RAW/'scores_D11F.npy')
ref_lin=np.load(Path(m.O)/RAW/'lin_D11F.npy')
w=winners(m.pop,linear)
wr=winners(m.pop,ref_lin)
checks={'spearman':float(spearmanr(scores,ref).statistic),'max_abs_score_error':float(np.max(np.abs(scores-ref))),'max_abs_linear_error':float(np.max(np.abs(linear-ref_lin))),'W_row_agreement':float(np.isin(w,wr).mean()),'threads':torch.get_num_threads(),'rows':len(scores)}
checks['pass']=checks['spearman']>=.99999 and checks['max_abs_score_error']<=.001 and checks['W_row_agreement']==1
(D/'verification.json').write_text(json.dumps(checks,indent=2))
print(RAW,checks,flush=True)
if not checks['pass']: raise RuntimeError('Frozen reconstruction tolerance failed. No attribution authorized.')
(D/'LABEL_BLIND_COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
