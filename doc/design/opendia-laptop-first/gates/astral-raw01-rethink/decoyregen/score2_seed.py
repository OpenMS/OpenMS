"""Post-hoc variant of score2.py (written after S failed to seed on raw 02): if fewer than 200 target winners pass the
semi-supervised q cut, the positives are the top 200 target winners by the current score instead. Applied to every
family alike. Everything else as score2.py; the predeclared rule is read on this variant only as a secondary result."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]
        pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
sys.argv = ['score2.py']
exec(open(os.path.dirname(os.path.abspath(__file__)) + '/score2.py').read().split('import oswlike_score as S', 1)[1])
