"""Overfitting checks for a single CNN -> MLP discriminator (raw 01, design run, diagnostic only).

1. Fold hygiene: are target and paired decoy always in the same fold (no pair leak)?
2. Memorisation: within each fitted model, do decoys it trained on score lower than held-out decoys? Measured as
   AUC(held-out decoy winners vs training decoy winners) on the model's own score; 0.5 = no memorisation.
3. External error check: entrapment FDP (Wen et al. combined, r = 0.9851) at q<=0.01 / 0.05 for LDA vs MLP ensemble,
   and an MLP on the minimal input (cnn_z only + 4 search features) to see whether the CNN alone feeds the gain.
Entrapment labels come from library protein accessions only and are never used for fitting.
"""
import sqlite3
import numpy as np, pandas as pd, torch
src = open('lag.py').read().split('res = {}')[0]
g = {}
exec(src, g)
df = g['df']
c = sqlite3.connect('file:' + g['LIB'] + '?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) all_ent
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID
  join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
df = df.merge(ent, on='precursor_id', how='left')
g['df'] = df
R_ENT = 148983 / 151239

# 1. fold hygiene
fk = df.groupby('pair_key').fold.nunique()
print(f'pairs spanning >1 fold: {(fk > 1).sum()} of {len(fk)}')

def fdp(best):
    out = []
    for q in (0.01, 0.05):
        t = best[(best.is_decoy == 0) & (best.q <= q)]
        ne = int(t.all_ent.sum()); nh = len(t) - ne
        out.append(f'q<={q}: N={len(t)} Nh={nh} Ne={ne} FDPc={ne * (1 + 1 / R_ENT) / max(1, len(t)):.4f}')
    return ' | '.join(out)

# 2. memorisation probe: refit one fold model explicitly to know which decoys it trained on
def memo(cols, kind, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    X = df[cols].to_numpy(np.float64)
    res = []
    for k in g['folds']:
        tr = np.where(g['fold_of'] != k)[0]
        score = X[:, 0].copy()
        for fdr in [0.15, 0.05, 0.05]:
            pos, neg = g['select'](tr, score, fdr)
            score = g['lda'](X, pos, neg)
        if kind == 'mlp':
            for _ in range(2):
                pos, neg = g['select'](tr, score, 0.05)
                score = g['mlp_fit'](X, pos, neg)
        te = np.where(g['fold_of'] == k)[0]
        held = g['winners_idx'](te, score); held = held[g['is_dec'][held]]
        res.append(g['auc'](pd.Series(score[held]), pd.Series(score[neg])))
    return np.round(res, 3)

SB = g['SEARCH'] + g['BOUND']
print('memorisation AUC(held-out decoys vs trained decoys), per fold:')
print('  lda', memo(SB, 'lda'), '  mlp', memo(SB, 'mlp'))

def ens(cols, seeds):
    acc = 0
    for s in seeds:
        torch.manual_seed(s); np.random.seed(s)
        sc = g['crossfit'](cols, 'mlp')
        for k in g['folds']:
            m = g['fold_of'] == k
            sc[m] = (sc[m] - sc[m].mean()) / sc[m].std()
        acc = acc + sc
    return acc

def report(name, sc):
    tmp = df[['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key', 'all_ent']].assign(score=sc)
    best, w = g['paired_pooled'](tmp, 'score')
    print(f'{name:34s} {fdp(best)}', flush=True)

report('tool stage2 (search+boundary LDA)', df.stage2_score.values if 'stage2_score' in df else g['crossfit'](SB, 'lda'))
report('LDA search+boundary (python)', g['crossfit'](SB, 'lda'))
report('MLP ens5 search+boundary', ens(SB, range(10, 15)))
report('MLP ens5 search only (CNN+4)', ens(g['SEARCH'], range(10, 15)))
