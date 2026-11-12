from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parent
pd.set_option('display.width',200)
pd.set_option('display.max_columns',25)
print('RT/FWHM CONFIDENCE SUBSETS')
print(pd.read_csv(O/'rt_fwhm_summary.tsv',sep='\t').round(4).to_string(index=False))
g=pd.read_csv(O/'group_summary.tsv',sep='\t')
fs=['transition_support','boundary_log_fwhm','boundary_log_width','boundary_coelution','f12_support','p_n_tight','p_ppm_abs_med','natE_apex_logp','natE_libcos','f12_libcos','xic_core_supported','xic_core_empty_time_fraction','xic_conditional_fragment_dropout','xic_single_bin_fragments','xic_corr_zerofill','xic_corr_occupied_time','xic_shared_zero_corr_inflation','xic_rebin_phase_sensitivity','xic_rebin_corr_change','xic_halfheight_occupied_width_s','xic_global_peak_outside_core','native_binned_coherence_gap','native_binned_libcos_gap']
print('\nGROUP MEDIANS')
print(g[g.feature.isin(fs)].pivot(index=['raw','feature'],columns='group',values='median').round(4).to_string())
print('\nSHAPE MEANS AND UPPER TAILS')
print(g[g.feature.isin(['xic_rebin_phase_sensitivity','xic_shared_zero_corr_inflation','xic_global_peak_outside_core','xic_core_empty_time_fraction'])][['raw','group','feature','finite','mean','q90']].round(4).to_string(index=False))
a=pd.read_csv(O/'score_stratified_auc.tsv',sep='\t')
print('\nCONTROLLED AUC, HIGH-VALUE ORIENTATION')
print(a[a.feature.isin(fs)&a.right.eq('DEC')].pivot(index=['raw','feature'],columns='left',values='auc_high').round(4).to_string())
print('\nCOVERAGE AND SCORE GAPS')
print(a[a.feature.eq('xic_corr_zerofill')].round(5).to_string(index=False))
print('\nSAMPLING CADENCE')
print((O/'sampling_cadence.tsv').read_text())
if (O/'common_score_weight_auc.tsv').exists():
 a=pd.read_csv(O/'common_score_weight_auc.tsv',sep='\t')
 print('\nCOMMON D-REFERENCE STRATUM WEIGHTS')
 print(a[a.reference.eq('D')].pivot(index=['raw','feature'],columns='group',values='auc_high').round(4).to_string())
 print('\nCOMMON STRATUM COVERAGE')
 print(a[a.feature.eq('xic_corr_zerofill')].round(4).to_string(index=False))
 print('\nHIGH CONFIDENCE AND MARGINAL, SAME STRATA CONTROLS')
 print(a[a.feature.isin(['natE_apex_logp','xic_corr_occupied_time','xic_rebin_phase_sensitivity'])].round(4).to_string(index=False))
