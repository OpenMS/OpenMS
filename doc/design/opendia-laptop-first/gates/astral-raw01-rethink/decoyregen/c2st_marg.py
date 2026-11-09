"""Direction of the C2ST separators: medians and AUC (entrapment vs decoy winners; also human-target winners as a
reference) for the top features, all winners and winners above the q<=0.20 cut. Also: are entrapment/decoy log_width
differences present among 'no-signal' winners (log_sn below its decoy median)? Family O, MLP seed 0."""
import os, sys, numpy as np
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
sc = S.crossfit(S.ALL, kind='mlp', seed=0)
w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc); t2 = sc[kept[(~S.dec[kept]) & (q <= 0.20)]].min()
dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
print(f'run {os.environ.get("OPENDIA_RUN")}')
for n in ('log_width', 'lib_corr', 'log_sn', 'intensity', 'xcorr_coel_w', 'lib_manhattan', 'norm_rt'):
    j = S.NAMES.index(n); a, b = S.X0[ew, j], S.X0[dw, j]
    ok = np.isfinite(a).all() and np.isfinite(b).all()
    auc = roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.nan_to_num(np.r_[a, b]))
    at, bt = a[sc[ew] >= t2], b[sc[dw] >= t2]
    print(f'  {n:14s} median e {np.nanmedian(a):.4f} d {np.nanmedian(b):.4f} AUC {auc:.3f} | above q<=0.20: e {np.nanmedian(at):.4f} d {np.nanmedian(bt):.4f}')
j = S.NAMES.index('log_sn'); thr = np.nanmedian(S.X0[dw, j]); jw = S.NAMES.index('log_width')
a = S.X0[ew, jw][S.X0[ew, j] < thr]; b = S.X0[dw, jw][S.X0[dw, j] < thr]
print(f'  low-S/N winners: log_width median e {np.median(a):.4f} d {np.median(b):.4f} AUC {roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.r_[a, b]):.3f}')
