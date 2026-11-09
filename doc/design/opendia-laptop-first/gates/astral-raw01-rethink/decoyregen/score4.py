"""Two confound checks after score3 (declared before running).
(a) Where is O's excess? Score O on the FULL set (MLP, seeds 0/1/2). Split pairs into 'shufshare' (the S decoy shares >=2
    top-6 fragments with its target, i.e. low-complexity or repetitive targets; dropped in score3) and 'clean'. At the
    full-set q<=0.05 cut, report S_e/S_d separately within each stratum (entrapment and decoy winners of that stratum).
    Reading: the excess is concentrated in shufshare if its ratio band lies above the clean one for >=2 of 3 seeds.
(b) MS1 borrowing. P and S decoys share their target's precursor m/z, so their MS1 subscores read the target's own
    precursor trace; O, entrapment and false targets do not. Rescore O/P/S on the clean set WITHOUT the MS1 group (MLP,
    seeds 0/1/2). Reading: P/S calibration (score3) is attributed to MS1 borrowing if without MS1 their median q<=0.05
    FDP_c rises to O's level (within 0.01 of O's median) on this run.
Usage: OPENDIA_RUN=R02."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
meta = S.m; fam = meta.family.values
u = meta[fam == 'S'][['precursor_id', 'pair']].drop_duplicates('precursor_id')
bad = set(q for p, q in zip(u.precursor_id, u.pair)
          if p in F and q in F and (np.abs(F[p][:, None] - F[q][None, :]) <= 34e-6 * F[p][:, None]).any(1).sum() >= 2)
isbad = np.isin(meta.pair.values, list(bad))
full = {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
def use(keep):
    for k, v in full.items(): setattr(S, k, v[keep])
    S.m = meta[keep].reset_index(drop=True)
def ratio(ke, ne, kd, nd):
    d = max(1e-12, kd / nd); return f'{ke / ne / d:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / d:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / d:.2f}] ({ke}/{kd})'
print('(a) O full set, stratified 5% tail', flush=True)
keep = (fam == 'T') | (fam == 'O'); use(keep); b = isbad[keep]
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min(); out = []
    for nm, msk in (('shufshare', b), ('clean', ~b)):
        ww = w[msk[w]]; dw = ww[S.dec[ww]]; ew = ww[(~S.dec[ww]) & S.entr[ww]]
        out.append(f'{nm} (pairs share {msk[w].mean():.3f}): {ratio(int((sc[ew] >= t).sum()), len(ew), int((sc[dw] >= t).sum()), len(dw))}')
    print(f'  seed{seed}: ' + ' | '.join(out), flush=True)
print('(b) clean set without MS1 group', flush=True)
NOMS1 = [c_ for c_ in S.ALL if c_ not in S.GROUPS['ms1']]
for Fm in ('O', 'P', 'S'):
    use((~isbad) & ((fam == 'T') | (fam == Fm)))
    for seed in (0, 1, 2):
        S.report(f'{Fm} clean noMS1 MLP seed{seed}', S.crossfit(NOMS1, kind='mlp', seed=seed * 10))
