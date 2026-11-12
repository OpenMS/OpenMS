"""Supplement fixed before reading primary results: alignment checks, boundary-width units,
confidence subsets, and equal score-stratum weights to expose reference-selection bias.
No tuning: reuse diagnose.py's width-0.5, >=5 finite rows, same-fold strata.
"""
from pathlib import Path
import numpy as np,pandas as pd
from diagnose import OUT,V4,V2,auc

def main():
    checks=[]; common=[]; widths=[]
    features=['natE_apex_logp','natU_apex_logp','f12_libcos','p_ppm_abs_med','xic_corr_zerofill','xic_corr_occupied_time','xic_shared_zero_corr_inflation','xic_rebin_phase_sensitivity','xic_rebin_corr_change','xic_conditional_fragment_dropout','xic_core_empty_time_fraction','native_binned_coherence_gap','native_binned_libcos_gap']
    for raw in ['raw01','raw02']:
        t=pd.read_csv(OUT/f'window_measurements_{raw}.tsv',sep='\t')
        old=pd.read_csv(V4/f'dd_feats_{raw}.tsv',sep='\t').drop_duplicates('row').set_index('row')
        a=t.set_index('row').loc[old.index]
        for new,prev in [('xic_core_supported','n_sig'),('xic_corr_zerofill','corr_med'),('xic_core_fraction','frac_core')]:
            v=abs(a[new]-old[prev]);err=v.max()
            assert err<1e-6,(raw,new,err)
            checks.append(dict(raw=raw,feature=new,old_feature=prev,n=len(v),finite=int(v.notna().sum()),max_abs_error=err))
        t['score_bin']=np.floor(t.s/.5).astype(int)
        t['D_high']=t.D&t.diann_high_conf
        t['D_marginal']=t.D&~t.diann_high_conf
        # D controls must share exactly the same fold/score stratum weights as HUMAN/ENT.
        for feature in features:
            for refgroup in ['D','D_high','D_marginal']:
                vals={g:[] for g in [refgroup,'HUMAN','ENT']};wt=[];coverage={g:0 for g in [refgroup,'HUMAN','ENT','DEC']}
                for key,g in t.groupby(['fold','score_bin']):
                    arms={n:g[g[n]&np.isfinite(g[feature])] for n in coverage}
                    if any(len(v)<5 for v in arms.values()):continue
                    wt.append(len(arms[refgroup]))
                    for name in vals:vals[name].append(auc(arms[name][feature],arms['DEC'][feature]))
                    for name in coverage:coverage[name]+=len(arms[name])
                for name,v in vals.items():
                    common.append(dict(raw=raw,reference=refgroup,group=name,feature=feature,auc_high=np.average(v,weights=wt) if wt else np.nan,n_strata=len(wt),reference_covered=coverage[refgroup],reference_total=int(t[refgroup].sum()),group_covered=coverage[name],dec_covered=coverage['DEC']))
        for group in ['D','E','HUMAN','ENT','DEC']:
            g=t[t[group]]
            # Stored logs are retained as logs. expm1 must not be called physical width without source verification.
            for feature in ['boundary_log_fwhm','boundary_log_width','xic_halfheight_occupied_width_s','xic12_halfheight_occupied_width_s','xic_rebin_phase_sensitivity','xic_shared_zero_corr_inflation']:
                widths.append(dict(raw=raw,group=group,feature=feature,n=len(g),median=g[feature].median(),q90=g[feature].quantile(.9),share_gt_0_1=(g[feature]>.1).mean()))
    pd.DataFrame(checks).to_csv(OUT/'numerical_reproduction.tsv',sep='\t',index=False)
    pd.DataFrame(common).to_csv(OUT/'common_score_weight_auc.tsv',sep='\t',index=False)
    pd.DataFrame(widths).to_csv(OUT/'width_phase_checks.tsv',sep='\t',index=False)
    print('PASS: existing dd_feats reproduced within 1e-6 on both raws')
    print('Common-stratum AUC and width diagnostics written')

if __name__=='__main__':main()
