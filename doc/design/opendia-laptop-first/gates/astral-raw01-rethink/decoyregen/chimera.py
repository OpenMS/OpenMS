"""Chimera test (k4). Hypothesis: false targets (entrapment) pick up fragment signal from co-eluting real peptides more
often than decoys do, so they survive into the score tail more often (S_e/S_d > 1).

Real-signal map: every DIA-NN precursor at Q.Value <= 0.01 in this run (reference/diagnostic only, never a training
signal). Its theoretical b/y ions (charge 1 and 2, carbamidomethyl C, 150-2000 Th) are taken as the fragments of a real
co-eluting peptide. Co-eluting means: the same isolation window (Astral 2 Th windows from 380 Th; approximate, the
real window bounds were not parsed), and the winner's apex within [RT.Start - 5 s, RT.Stop + 5 s] of that precursor.
The precursor itself (same modified sequence and charge) is excluded.
k = the number of the winner's top-6 library fragments that lie within +-17 ppm (the extraction half-width) of any
fragment of a co-eluting real precursor.

Scorer: family O (the library's decoys) on the full set, MLP seeds 0/1/2, with seed fallback as in score2_seed.py.
Winners = best candidate per precursor. Reported:
  (1) mean k and the share k>=1 / k>=2 for entrapment winners against decoy winners (all winners);
  (2) S_e/S_d at the q<=0.05 cut, within strata of k (0, 1, >=2);
  (3) S_e/S_d at the q<=0.05 cut with decoy winners reweighted to the entrapment winners' k distribution.
Reading, declared before running: the chimera hypothesis is supported on a run if, for >=2 of 3 seeds, the reweighted
5% ratio (3) removes at least half of the excess over 1 (r_w - 1 <= (r - 1) / 2) AND its band includes 1. It is
supported overall only if it holds on both runs. It is rejected on a run if entrapment and decoy winners differ by less
than 0.02 in the k>=1 share (no mechanism), or (3) leaves more than half of the excess for >=2 seeds.
Usage: OPENDIA_RUN=R02 (DIA-NN run suffix taken from it)."""
import os, sys, re, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S

RUN = os.environ['OPENDIA_RUN']; RUNNO = RUN[-2:]
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select

# ---- theoretical fragments of DIA-NN 1% precursors in this run
AA = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=103.00919 + 57.02146, L=113.08406,
          I=113.08406, N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891,
          F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
PROTON, H2O = 1.007276, 18.010565
def frags(seq):
    r = np.array([AA[a] for a in seq]); b = np.cumsum(r)[:-1]; y = np.cumsum(r[::-1])[:-1] + H2O
    f1 = np.r_[b + PROTON, y + PROTON]; f = np.r_[f1, (f1 + PROTON) / 2]
    return f[(f >= 150) & (f <= 2000)]
dn = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                     columns=['Run', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Charge', 'Precursor.Mz', 'RT.Start', 'RT.Stop', 'Q.Value'])
dn = dn[dn.Run.str.endswith('_' + RUNNO) & (dn['Q.Value'] <= 0.01)].reset_index(drop=True)
dn['key'] = dn['Modified.Sequence'].str.replace('(UniMod:4)', '', regex=False) + '/' + dn['Precursor.Charge'].astype(str)
dn['win'] = np.floor((dn['Precursor.Mz'] - 380) / 2).astype(int)
DF = [frags(s) for s in dn['Stripped.Sequence']]
lo = dn['RT.Start'].values * 60 - 5; hi = dn['RT.Stop'].values * 60 + 5
bywin = {w: g.index.values for w, g in dn.groupby('win')}
print(f'run {RUN}: DIA-NN 1% precursors {len(dn)}', flush=True)

# ---- winner library fragments, m/z, identity
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
pp = pd.read_sql('select p.ID pid, p.PRECURSOR_MZ mz, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID', c).drop_duplicates('pid')
pp['key'] = pp.ms.str.replace('(UniMod:4)', '', regex=False) + '/' + pp.z.astype(str)
MZ = dict(zip(pp.pid, pp.mz)); KEY = dict(zip(pp.pid, pp.key))

def kshare(pid, apex):
    f = F.get(pid)
    if f is None: return 0
    w = int(np.floor((MZ[pid] - 380) / 2)); idx = bywin.get(w)
    if idx is None: return 0
    idx = idx[(lo[idx] <= apex) & (hi[idx] >= apex)]
    me = KEY.get(pid)
    allf = [DF[i] for i in idx if dn.key.values[i] != me]
    if not allf: return 0
    a = np.sort(np.concatenate(allf)); j = np.searchsorted(a, f)
    d = np.minimum(np.abs(a[np.clip(j, 0, len(a) - 1)] - f), np.abs(a[np.clip(j - 1, 0, len(a) - 1)] - f))
    return int((d <= 17e-6 * f).sum())

meta = S.m; fam = meta.family.values
keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
S.m = meta[keep].reset_index(drop=True); apexv = S.m.apex.values

def band(ke, ne, sd):
    return f'{ke / ne / sd:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / sd:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / sd:.2f}]'
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    kd = np.array([kshare(S.prec[i], apexv[i]) for i in dw]); ke = np.array([kshare(S.prec[i], apexv[i]) for i in ew])
    st = lambda k: np.minimum(k, 2)
    print(f'seed{seed}: winners e {len(ew)} d {len(dw)} | mean k e {ke.mean():.3f} d {kd.mean():.3f} | k>=1 e {(ke >= 1).mean():.3f} d {(kd >= 1).mean():.3f} | k>=2 e {(ke >= 2).mean():.3f} d {(kd >= 2).mean():.3f}', flush=True)
    se, sd_ = sc[ew] >= t, sc[dw] >= t
    r = (se.sum() / len(ew)) / (sd_.sum() / len(dw))
    line = [f'pooled {band(int(se.sum()), len(ew), sd_.sum() / len(dw))}']
    for s_ in (0, 1, 2):
        me_, md_ = st(ke) == s_, st(kd) == s_
        if md_.sum() and (sd_ & md_).sum():
            line.append(f'k={s_ if s_ < 2 else ">=2"} ({int((se & me_).sum())}/{int((sd_ & md_).sum())}) {band(int((se & me_).sum()), me_.sum(), (sd_ & md_).sum() / md_.sum())}')
    # reweight decoys to the entrapment k-stratum distribution
    pe = np.array([(st(ke) == s_).mean() for s_ in (0, 1, 2)]); pd_ = np.array([(st(kd) == s_).mean() for s_ in (0, 1, 2)])
    wts = (pe / np.maximum(pd_, 1e-12))[st(kd)]
    sdw = (wts * sd_).sum() / wts.sum()
    rw = (se.sum() / len(ew)) / sdw
    line.append(f'reweighted {band(int(se.sum()), len(ew), sdw)} | excess removed {(r - rw) / max(1e-9, r - 1):.2f}')
    print('   5% cut: ' + ' | '.join(line), flush=True)
