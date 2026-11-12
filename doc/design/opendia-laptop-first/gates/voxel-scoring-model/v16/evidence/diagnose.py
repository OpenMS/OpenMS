"""Frozen descriptive measurement audit. No fitting, threshold optimization or label-defined features.
Execute from anywhere with the supplied Python. All writes are next to this file.
Fixed diagnostic comparisons: same fold, score bins width 0.5 anchored at zero,
minimum 5 finite rows per arm. AUC uses midranks. No area/support matching.
High-confidence DIA-NN: Q.Value<=0.001 and PEP<=0.01, marginal is its complement.
These descriptive bins are not acceptance rules. Raw02 is descriptive replication.
"""
from pathlib import Path
import hashlib,json,datetime,warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata
OUT=Path(__file__).resolve().parent
V4=OUT.parents[1]
V2=V4.parent
ROOT=V2.parents[1]
H_NAMES=['transition_support','log_matched_scans','log_peak_width','abs_library_rt_error','boundary_log_fwhm','boundary_log_width','boundary_coelution','boundary_library_correlation','boundary_apex_agreement','byion_b','byion_y']
F_NAMES=['f12_support','f12_coelution','f12_libcos','f12_extra_support','f12_extra_coelution','f12_extra_share']
P_NAMES=['p_n_tight','p_ppm_abs_med','p_tight_share','c_med_top3','c_med_k']
N_NAMES=['apex_logp','sum_logp','scans_logp_ge2','repeat_logp','libcos']

def auc(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    a=a[np.isfinite(a)];b=b[np.isfinite(b)]
    if not len(a) or not len(b): return np.nan
    return float((rankdata(np.r_[a,b])[:len(a)].sum()-len(a)*(len(a)+1)/2)/(len(a)*len(b)))

def paircorr(x,mask=None):
    # Median pairwise Pearson, excluding constant fragment profiles.
    if mask is None: mask=np.ones((len(x),x.shape[2]),bool)
    n=mask.sum(1)[:,None]
    mu=(x*mask[:,None,:]).sum(2)/np.maximum(n,1)
    y=(x-mu[:,:,None])*mask[:,None,:]
    norm=np.sqrt((y*y).sum(2))
    c=np.einsum('nft,ngt->nfg',y,y)
    den=norm[:,:,None]*norm[:,None,:]
    c=np.divide(c,den,out=np.full_like(c,np.nan),where=den>1e-12)
    i,j=np.triu_indices(x.shape[1],1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        return np.nanmedian(c[:,i,j],axis=1)

def profiles(x):
    x=np.asarray(x,float); core=x[:,:,12:20]; summed=x.sum(1)
    occ=summed>0; oc=core.sum(1)>0
    area=core.sum(2); support=(area>0).sum(1)
    z={}
    z['core_supported']=support
    z['core_empty_time_fraction']=1-oc.mean(1)
    # Conditional absence among fragments that exist in the core, at times with any signal.
    denom=support*oc.sum(1)
    z['conditional_fragment_dropout']=np.divide(denom-(core>0).sum((1,2)),denom,out=np.full(len(x),np.nan),where=denom>0)
    z['single_bin_fragments']=((core>0).sum(2)==1).sum(1)
    z['corr_zerofill']=paircorr(x)
    z['corr_occupied_time']=paircorr(x,occ)
    z['shared_zero_corr_inflation']=z['corr_zerofill']-z['corr_occupied_time']
    z['corr_rebin2_origin0']=paircorr(x.reshape(len(x),x.shape[1],16,2).sum(3))
    z['corr_rebin2_origin1']=paircorr(x[:,:,1:31].reshape(len(x),x.shape[1],15,2).sum(3))
    z['rebin_phase_sensitivity']=abs(z['corr_rebin2_origin0']-z['corr_rebin2_origin1'])
    z['rebin_corr_change']=z['corr_rebin2_origin0']-z['corr_zerofill']
    # Occupied bins above half total-profile maximum, not a native interpolated FWHM.
    z['halfheight_occupied_width_s']=((summed>=summed.max(1)[:,None]*.5)&occ).sum(1)*30/32
    z['core_peak_offset_s']=(np.argmax(summed[:,12:20],axis=1)+12+.5)*30/32-15
    z['global_peak_outside_core']=((np.argmax(summed,axis=1)<12)|(np.argmax(summed,axis=1)>=20)).astype(float)
    z['core_fraction']=area.sum(1)/np.maximum(summed.sum(1),1e-12)
    return pd.DataFrame(z)

def stratified(tab,features,raw):
    records=[]; counts=[]
    bins=np.floor(tab.s.to_numpy()/0.5).astype(int)
    tab=tab.assign(score_bin=bins)
    for left,right in [('D','DEC'),('HUMAN','DEC'),('ENT','DEC'),('D','ENT'),('E','DEC')]:
        a=tab[tab[left]];b=tab[tab[right]]
        ag=list(a.groupby(['fold','score_bin']))
        bg=dict(tuple(b.groupby(['fold','score_bin'])))
        for feature in features:
            vals=[]; weights=[]; covered_a=0; covered_b=0; gaps=[]
            for key,aa in ag:
                bb=bg.get(key,b.iloc[:0])
                if len(aa)<5 or len(bb)<5:continue
                aa=aa[np.isfinite(aa[feature])];bb=bb[np.isfinite(bb[feature])]
                if len(aa)<5 or len(bb)<5:continue
                vals.append(auc(aa[feature],bb[feature]));weights.append(len(aa))
                covered_a+=len(aa);covered_b+=len(bb);gaps.append(abs(aa.s.mean()-bb.s.mean()))
            records.append(dict(raw=raw,left=left,right=right,feature=feature,auc_high=np.average(vals,weights=weights) if weights else np.nan,n_left=len(a),n_right=len(b),finite_left=int(np.isfinite(a[feature]).sum()),finite_right=int(np.isfinite(b[feature]).sum()),matched_left=covered_a,matched_right=covered_b,coverage_left=covered_a/max(len(a),1),coverage_right=covered_b/max(len(b),1),strata=len(vals),weighted_mean_score_gap=np.average(gaps,weights=weights) if weights else np.nan))
        for key,aa in a.groupby(['fold','score_bin']):
            bb=b[(b.fold==key[0])&(b.score_bin==key[1])]
            counts.append(dict(raw=raw,left=left,right=right,fold=key[0],score_bin=key[1],n_left=len(aa),n_right=len(bb)))
    return pd.DataFrame(records),pd.DataFrame(counts)

def tests():
    assert auc([1,1],[1,1])==.5
    assert auc([2,3],[0,1])==1
    assert auc([0,1],[2,3])==0
    assert auc([0,1],[1,2])==.125
    assert np.isnan(auc([np.nan],[1]))
    x=np.zeros((2,6,32));x[:,:,14:17]=np.array([1,2,1])
    f=profiles(x)
    assert np.allclose(f.corr_zerofill,1)
    assert np.allclose(f.conditional_fragment_dropout,0)
    assert np.allclose(f.core_empty_time_fraction,5/8)
    assert np.allclose(f.halfheight_occupied_width_s,3*30/32)
    t=pd.DataFrame({'s':[0.1]*10+[0.1]*10,'fold':[0]*10+[1]*10,'D':[True]*10+[False]*10,'DEC':[False]*10+[True]*10,'HUMAN':[True]*10+[False]*10,'ENT':[False]*20,'E':[False]*20,'v':range(20)})
    r,_=stratified(t,['v'],'synthetic')
    assert (r.matched_left==0).all(), 'Cross-fold matches prohibited'
    print('TESTS PASS: ties, direction, missingness, support, width, fold isolation',flush=True)

def main():
    tests()
    import sqlite3
    con=sqlite3.connect('file:'+str(ROOT/'astral-entrapment-20261107/lib/library.pqp')+'?mode=ro',uri=True)
    ent=pd.read_sql("select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID",con).set_index('precursor_id').e
    summaries=[]; rtstats=[]; aucs=[]; match=[]; inventory=[]; consistency=[]
    expected={'raw01':(866,1413),'raw02':(929,1362)}
    for raw in ['raw01','raw02']:
        print('BEGIN',raw,flush=True)
        folder=V2/raw
        pop=pd.read_csv(folder/'population.tsv',sep='\t')
        s=np.load(folder/'scores_D11FWT.npy')
        assert len(s)==len(pop) and np.isfinite(s).all() and pop.sample_id.is_unique
        w=pop.assign(s=s,row=np.arange(len(pop))).sort_values(['s','sample_id'],ascending=[False,True]).drop_duplicates('precursor_id')
        w=w[w.s>-1e20].copy().sort_values('row').reset_index(drop=True)
        stage=pd.read_csv(V2/f'v3seq/misses_D11FWT_{raw}.tsv',sep='\t').set_index('precursor_id').stage
        w['D']=w.precursor_id.map(stage).eq('D right window, below cut')
        w['E']=w.precursor_id.map(stage).eq('E accepted')
        w['DEC']=w.is_decoy.eq(1)
        labels=w.precursor_id.map(ent)
        assert labels[w.is_decoy==0].notna().all()
        w['ENT']=w.is_decoy.eq(0)&labels.eq(1)
        w['HUMAN']=w.is_decoy.eq(0)&labels.eq(0)
        assert (int(w.D.sum()),int(w.E.sum()))==expected[raw]
        rows=w.row.to_numpy();features=[]
        specs=[('heuristics',H_NAMES),('feat12',F_NAMES),('featPC',P_NAMES),('nat',[f'nat_{z}' for z in N_NAMES]),('natE',[f'natE_{z}' for z in N_NAMES]),('natU',[f'natU_{z}' for z in N_NAMES]+['natU_n_unique'])]
        for filename,names in specs:
            arr=np.load(folder/f'{filename}.npy',mmap_mode='r'); assert len(arr)==len(pop)
            val=np.asarray(arr[rows],float);assert val.shape[1]==len(names)
            w[names]=val;features.extend(names)
            inventory.append(dict(raw=raw,file=str(folder/f'{filename}.npy'),shape=str(arr.shape),bytes=(folder/f'{filename}.npy').stat().st_size,mtime_ns=(folder/f'{filename}.npy').stat().st_mtime_ns))
        # Existing tabular measurements, verify row identity instead of trusting concatenation.
        for name in ['dd_feats','dd_band']:
            old=pd.read_csv(V4/f'{name}_{raw}.tsv',sep='\t')
            old=old.drop_duplicates('row').set_index('row')
            if 'precursor_id' in old:
                assert np.array_equal(pop.iloc[old.index].precursor_id,old.precursor_id)
            consistency.append(dict(raw=raw,check=name,n_rows=len(old),n_winners=int(np.isin(rows,old.index).sum())))
        for name in ['xic','xic12']:
            arr=np.load(folder/f'{name}.npy',mmap_mode='r');assert len(arr)==len(pop)
            parts=[]
            for start in range(0,len(rows),10000):
                parts.append(profiles(arr[rows[start:start+10000]]))
            f=pd.concat(parts,ignore_index=True).add_prefix(name+'_')
            w=pd.concat([w,f],axis=1);features.extend(f.columns)
        w['native_binned_coherence_gap']=w.boundary_coelution-w.xic_corr_zerofill
        w['native_binned_libcos_gap']=w.natE_libcos-w.f12_libcos
        features+=['native_binned_coherence_gap','native_binned_libcos_gap']
        dn=pd.read_csv(V4/f'dd_diann_{raw}.tsv',sep='\t')
        assert dn.row.is_unique
        assert np.array_equal(pop.iloc[dn.row].precursor_id,dn.precursor_id)
        cols=['row','rt_diff','fwhm_s','PEP','Q.Value','best_fr_rank','RT.Start','RT.Stop','RT','lq']
        w=w.merge(dn[cols],on='row',how='left',validate='one_to_one')
        w['abs_rt_fwhm']=abs(w.rt_diff)/w.fwhm_s.where(w.fwhm_s>0)
        w['abs_rt_s']=abs(w.rt_diff)
        w['diann_high_conf']=w['Q.Value'].le(.001)&w.PEP.le(.01)
        w['outside_diann_boundaries']=(w.rt_apex<w['RT.Start']*60)|(w.rt_apex>w['RT.Stop']*60)
        # Read geometry only, no assumptions that this is equivalent to native FWHM.
        w['diann_fwhm_bins']=w.fwhm_s/(30/32)
        w.to_csv(OUT/f'window_measurements_{raw}.tsv',sep='\t',index=False)
        for group in ['D','E','HUMAN','ENT','DEC']:
            g=w[w[group]]
            for f in features:
                v=g[f].replace([np.inf,-np.inf],np.nan).dropna()
                summaries.append(dict(raw=raw,group=group,feature=f,n=len(g),finite=len(v),q10=v.quantile(.1),median=v.median(),q90=v.quantile(.9),mean=v.mean()))
        for group in ['D','E']:
            for subset in ['all','high_conf','marginal']:
                g=w[w[group]& (True if subset=='all' else w.diann_high_conf if subset=='high_conf' else ~w.diann_high_conf)]
                good=g[np.isfinite(g.abs_rt_fwhm)]
                rtstats.append(dict(raw=raw,group=group,subset=subset,n=len(g),valid=len(good),median_fwhm_s=good.fwhm_s.median(),rt_abs_q50=good.abs_rt_s.median(),rt_abs_q90=good.abs_rt_s.quantile(.9),rt_fwhm_q50=good.abs_rt_fwhm.median(),rt_fwhm_q90=good.abs_rt_fwhm.quantile(.9),over_half_fwhm=(good.abs_rt_fwhm>.5).mean(),over_one_fwhm=(good.abs_rt_fwhm>1).mean(),over_two_fwhm=(good.abs_rt_fwhm>2).mean(),over_10s=(good.abs_rt_s>10).mean(),outside_diann_boundaries=good.outside_diann_boundaries.mean(),median_qty_log10=good.lq.median(),median_native_coherence=good.boundary_coelution.median(),median_binned_coherence=good.xic_corr_zerofill.median(),median_core_dropout=good.xic_conditional_fragment_dropout.median(),best_fragment_not_in_six=good.best_fr_rank.eq(0).mean()))
        a,b=stratified(w,features,raw);aucs.append(a);match.append(b)
        print('DONE',raw,'winners',len(w),'D/E',w.D.sum(),w.E.sum(),flush=True)
    pd.DataFrame(summaries).to_csv(OUT/'group_summary.tsv',sep='\t',index=False)
    pd.DataFrame(rtstats).to_csv(OUT/'rt_fwhm_summary.tsv',sep='\t',index=False)
    pd.concat(aucs).to_csv(OUT/'score_stratified_auc.tsv',sep='\t',index=False)
    pd.concat(match).to_csv(OUT/'score_strata_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(inventory).to_csv(OUT/'array_inventory.tsv',sep='\t',index=False)
    pd.DataFrame(consistency).to_csv(OUT/'alignment_checks.tsv',sep='\t',index=False)
    # Raw01 acquisition cadence, metadata only, already extracted and unit-checked in V13b.
    idx=pd.read_parquet(ROOT/'dd-plots-20261111/scan_index.parquet');idx=idx[idx.ms==2].sort_values('rt')
    dt=idx.groupby(['lo','hi']).rt.diff().dropna()
    pd.Series({'raw':'raw01','n_intervals':len(dt),'median_revisit_s':dt.median(),'q10_s':dt.quantile(.1),'q90_s':dt.quantile(.9),'bin_width_s':30/32}).to_csv(OUT/'sampling_cadence.tsv',sep='\t')
    plot(pd.DataFrame(rtstats),pd.DataFrame(summaries),pd.concat(aucs))

def plot(rt,gs,aa):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for raw,style in [('raw01','-'),('raw02','--')]:
        r=rt[(rt.raw==raw)&(rt.subset=='all')]
        for grp in ['D','E']:
            z=r[r.group==grp].iloc[0]
            axs[0,0].plot([.5,1,2],[z.over_half_fwhm,z.over_one_fwhm,z.over_two_fwhm],style,marker='o',label=raw+' '+grp)
    axs[0,0].set(xlabel='Absolute RT error / DIA-NN FWHM threshold',ylabel='Fraction exceeding threshold',title='RT disagreement, all D/E rows');axs[0,0].legend()
    groups=['D','E','HUMAN','ENT','DEC']
    for raw,style in [('raw01','-'),('raw02','--')]:
        for f,c in [('xic_corr_zerofill','C0'),('xic_corr_occupied_time','C1'),('boundary_coelution','C2')]:
            z=gs[(gs.raw==raw)&(gs.feature==f)].set_index('group').reindex(groups)
            axs[0,1].plot(groups,z['median'],style,color=c,marker='.',label=raw+' '+f)
    axs[0,1].set(title='Median coherence: differing time representations',ylabel='Correlation');axs[0,1].legend(fontsize=6)
    for raw,style in [('raw01','-'),('raw02','--')]:
        for f,c in [('xic_core_empty_time_fraction','C0'),('xic_conditional_fragment_dropout','C1')]:
            z=gs[(gs.raw==raw)&(gs.feature==f)].set_index('group').reindex(groups)
            axs[1,0].plot(groups,z['median'],style,color=c,marker='.',label=raw+' '+f)
    axs[1,0].set(title='Empty time bins versus fragment dropout',ylabel='Median fraction');axs[1,0].legend(fontsize=6)
    fs=['natE_apex_logp','xic_corr_zerofill','xic_corr_occupied_time','xic_shared_zero_corr_inflation','xic_rebin_phase_sensitivity','native_binned_coherence_gap']
    for raw,mark in [('raw01','o'),('raw02','x')]:
        for left,c in [('D','C0'),('HUMAN','C1'),('ENT','C2')]:
            z=aa[(aa.raw==raw)&(aa.left==left)&(aa.right=='DEC')].set_index('feature').reindex(fs)
            axs[1,1].plot(np.arange(len(fs)),z.auc_high,marker=mark,color=c,linestyle='-' if raw=='raw01' else '--',label=raw+' '+left)
    axs[1,1].axhline(.5,color='gray',lw=1);axs[1,1].set(xticks=np.arange(len(fs)),xticklabels=['native L','zero corr','occupied corr','zero inflation','phase sensitivity','native gap'],ylabel='Tie-correct AUC, high-value orientation',title='Same-fold score strata, coverage in table');axs[1,1].tick_params(axis='x',rotation=35);axs[1,1].legend(fontsize=6)
    fig.savefig(OUT/'mechanisms.png',dpi=180);fig.savefig(OUT/'mechanisms.pdf');plt.close(fig)

if __name__=='__main__':
    import sys
    if '--test' in sys.argv:tests()
    else:main()
