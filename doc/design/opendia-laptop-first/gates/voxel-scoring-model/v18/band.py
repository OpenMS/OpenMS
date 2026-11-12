import sys, pandas as pd, numpy as np
raw = sys.argv[1]
W = pd.read_csv(f'sibling_{raw}.tsv.gz', sep='\t'); W['sz'] = W.sib_z_near.replace(-np.inf, -50).fillna(-99)
hum = W.grp.isin(['D', 'E', 'H'])
for lo, hi in [(0.99, 0.999), (0.999, 0.9999), (0.99, 0.9999)]:
    s = (W.full_pct >= lo) & (W.full_pct < hi) & (W.sz > 2)
    print(f'[{lo},{hi}) sib z>2: human {int((s & hum).sum())} ent {int((s & (W.grp == "ENT")).sum())} dec {int((s & (W.grp == "DEC")).sum())} D {int((s & (W.grp == "D")).sum())}')
