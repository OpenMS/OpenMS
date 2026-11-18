"""Disjoint OpenSWATH sub-score groups added to OpenDIA's 24 features; labelled GBM ceiling (diagnostic)."""
import pandas as pd, numpy as np, sqlite3
from sklearn.ensemble import HistGradientBoostingClassifier
R='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published'
c = sqlite3.connect(f'file:{R}/svm.osw?mode=ro', uri=True)
cols=[r[1] for r in c.execute('pragma table_info(FEATURE_MS2)') if r[1].startswith('VAR_')]
X = pd.read_sql(f"select f.PRECURSOR_ID pid, {','.join('m.'+x for x in cols)} from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID join FEATURE_MS2 m on m.FEATURE_ID=f.ID where s.RANK=1", c)
fs = pd.read_csv('img/funnel_imgrid_svm.tsv', sep='\t'); fd = pd.read_csv('img/funnel_imgrid_diann.tsv', sep='\t'); pos=set(fd.pid)&set(fs.pid)
F = pd.read_csv('img/pre_imgrid/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.features.tsv', sep='\t')
oc=[x for x in F.columns if x not in ('candidate_row','precursor_id','is_decoy','fold','rt_apex','score','q_value')]
w=F.sort_values('score',ascending=False).drop_duplicates('precursor_id')
M=w.merge(X, left_on='precursor_id', right_on='pid', how='inner')
yy=np.where(M.is_decoy==1,0,np.where(M.precursor_id.isin(pos),1,-1)); dec=M.is_decoy.values
def n_at(s):
    o=np.argsort(-s); d=np.cumsum(dec[o]); t=np.cumsum(1-dec[o]); q=d/np.maximum(t,1); qm=np.minimum.accumulate(q[::-1])[::-1]; return int(((qm<=0.01)&(dec[o]==0)).sum())
def ceil(cc):
    s=np.zeros(len(M))
    for k in range(3):
        tr=(M.fold!=k).values&(yy>=0); te=(M.fold==k).values
        s[te]=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,random_state=0).fit(M.loc[tr,cc].fillna(0).values,yy[tr]).decision_function(M.loc[te,cc].fillna(0).values)
    return n_at(s)
groups = {
 'elution model fit': ['VAR_ELUTION_MODEL_FIT_SCORE'],
 'isotope (fragment) corr/overlap': ['VAR_ISOTOPE_CORRELATION_SCORE','VAR_ISOTOPE_OVERLAP_SCORE'],
 'xcorr coelution (+w)': ['VAR_XCORR_COELUTION','VAR_XCORR_COELUTION_WEIGHTED'],
 'xcorr shape (+w)': ['VAR_XCORR_SHAPE','VAR_XCORR_SHAPE_WEIGHTED'],
 'intensity + log S/N': ['VAR_INTENSITY_SCORE','VAR_LOG_SN_SCORE'],
 'mutual information (3)': [c for c in cols if 'MI_' in c and 'IM_' not in c],
 'mass deviation (+w)': ['VAR_MASSDEV_SCORE','VAR_MASSDEV_SCORE_WEIGHTED'],
}
base=ceil(oc); print('OpenDIA 24 alone:', base)
for g, cc in groups.items():
    cc=[x for x in cc if x in M.columns]; n=ceil(oc+cc); print(f'+ {g:34s} {cc}: {n} ({n-base:+d})')
allg=sum(groups.values(), []); n=ceil(oc+allg); print(f'+ all of the above ({len(allg)}): {n} ({n-base:+d})')
