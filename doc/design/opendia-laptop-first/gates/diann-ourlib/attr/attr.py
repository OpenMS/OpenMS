"""Which inputs make the D11F nonlinear part push DIA-NN-confirmed targets down? Raw01, exact V16 checkpoints (round 2,
per fold). For each held-out winner: nonlinear contribution n = tanh(g(z)).w[17:]; for each feature j, n_j = n with z_j
replaced by the fold's decoy-winner median. Effect_j = n - n_j (negative = feature j lowers the score). Groups:
DIA-NN(ourlib)-only at our W window (DN), both-accepted (BOTH), decoy winners (DEC), entrapment winners (ENT)."""
import sys, numpy as np, pandas as pd, torch, importlib.util, os
V='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
CK=V+'/v4/stageD_v16/score/original_threads/raw01'
os.environ['TRAIN_THREADS']='10'; sys.argv=['train_v4.py','raw01','D11F']
spec=importlib.util.spec_from_file_location('t',V+'/v4/train_v4.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
Hx=np.c_[m.H,m.F12]; names=pd.read_csv(V+'/v4/stageD_v16/score/original_threads/feature_names.tsv',sep='\t').feature.tolist()
W=pd.read_csv(V+'/v4/stageD_v16/score/reference_winners_raw01.tsv.gz',sep='\t')
dn=pd.read_csv('../diann_ourlib_raw01.tsv.gz',sep='\t'); dn=dn[dn['Q.Value']<=0.01]
dn=dn.set_index('precursor_id')
W['dn']=W.precursor_id.isin(dn.index)
W['dn_rt']=W.precursor_id.map(dn.RT*60)
W['at_rt']=(W.rt_apex-W.dn_rt).abs()<=10
grp=np.full(len(W),'',object)
grp[(W.dn&W.at_rt&(W.full_pct<0.999)).values]='DN'
grp[(W.dn&W.at_rt&(W.full_pct>=0.999)).values]='BOTH'
grp[(W.group=='DEC').values]='DEC'; grp[(W.group=='ENT').values]='ENT'
W['g']=grp; W=W[W.g!='']
rows=[]
for k in range(3):
    ck=torch.load(f'{CK}/checkpoint_fold{k}_r2.pt',weights_only=False)
    net=m.Net(False,True,17); net.load_state_dict(ck['state_dict']); net.eval(); w=ck['readout'].double()
    Z=(Hx-ck['mu'])/ck['sd']
    Wk=W[W.fold==k]; r=Wk.row.values; z=torch.from_numpy(Z[r].astype(np.float32))
    dmed=np.median(Z[W[(W.fold==k)&(W.g=='DEC')].row.values],0)
    with torch.no_grad():
        nl=lambda zz:(torch.tanh(net.g(zz)).double()@w[17:]).numpy()
        lin=(z.double()@w[:17]).numpy(); n0=nl(z)
        eff=np.zeros((len(r),17))
        for j in range(17):
            z2=z.clone(); z2[:,j]=float(dmed[j]); eff[:,j]=n0-nl(z2)
    d=pd.DataFrame(eff,columns=names); d['g']=Wk.g.values; d['nl']=n0; d['lin']=lin; d['fold']=k
    d=d.join(pd.DataFrame(Z[r],columns=['z_'+n for n in names]))
    rows.append(d)
D=pd.concat(rows,ignore_index=True)
# scale: decoy-winner SD of full score per fold
for k in range(3):
    x=D[(D.fold==k)]; sd=(x[x.g=='DEC'].nl+x[x.g=='DEC'].lin).std()
    D.loc[D.fold==k,names+['nl','lin']]/=sd
print('groups',D.g.value_counts().to_dict())
print('\nmedian linear / nonlinear contribution (decoy-SD units):'); print(D.groupby('g')[['lin','nl']].median().round(2).to_string())
print('\nmean effect of each feature on the nonlinear part vs decoy-median value (decoy-SD units):')
print(D.groupby('g')[names].mean().T.round(2).to_string())
print('\nmedian z of each feature:')
print(D.groupby('g')[['z_'+n for n in names]].median().T.round(2).to_string())
D.to_csv('attr_raw01.tsv.gz',sep='\t',index=False)
