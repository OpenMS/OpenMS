"""q-value calibration diagnosis (evaluation only; reads entrapment labels; no model, score or rule is changed).
Question: at q <= 0.01 the xic D11FWT arm has entrapment FDP_c 1.6-1.7 % (rule 1c limit 1.5 %). Where do decoys
under-represent false targets?

Calibration ratio of a set of accepted target winners: rho = E (1 + 1/r) / D, where E = entrapment target winners and
D = paired decoy winners above the same standardised score cut (same TDC as evaluate_v2.py). rho = 1 means the decoys
count false targets correctly; rho > 1 means q is optimistic. Poisson 90 % interval on rho from E and D.
Usage: diagnose.py <arm> <root of population/scores> <tag> [search features root pattern]"""
import json, sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import beta as betad
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pre = pd.read_sql('select ID precursor_id, CHARGE charge, PRECURSOR_MZ mz, LIBRARY_RT library_rt from PRECURSOR', c).set_index('precursor_id')

ARM, ROOT, TAG = sys.argv[1:4]
FEAT = {'xic': B + '/xicids-20261112/{raw}/pipeline/search.features.tsv',
        'old': {'raw01': B + '/astral-design-raw01-20261107/pipeline/search.features.tsv',
                'raw02': B + '/astral-entrapment-raw02-20261108/search.features.tsv'}}[TAG]
CAL = {'xic': B + '/xicids-20261112/{raw}/cal/cal.json',
       'old': {'raw01': B + '/astral-design-raw01-20261107/window300k/cal.json', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/cal.json'}}[TAG]


def path(p, raw):
    return p.format(raw=raw) if isinstance(p, str) else p[raw]


def tdc(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return w, best


def rho(E, D):
    """rho and a 90 % interval from the binomial split of E vs D (conditional on E + D)."""
    k = 1 + 1 / R
    if E + D == 0:
        return np.nan, np.nan, np.nan
    lo = betad.ppf(0.05, E, D + 1) if E > 0 else 0.0
    hi = betad.ppf(0.95, E + 1, D) if D > 0 else 1.0
    f = lambda p: k * p / (1 - p) if p < 1 else np.inf
    return (k * E / D if D else np.inf), f(lo), f(hi)


def fmt(E, D, N=None):
    r_, lo, hi = rho(E, D)
    n = f' N {N:5d}' if N is not None else ''
    return f'{n} E {E:4d} D {D:4d} rho {r_:5.2f} [{lo:4.2f}, {hi:5.2f}]'


for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{ROOT}/{raw}/population.tsv', sep='\t')
    s = np.load(f'{ROOT}/{raw}/scores_{ARM}.npy').astype(np.float64)
    feat = pd.read_csv(path(FEAT, raw), sep='\t').rename(columns={'candidate_row': 'sample_id'}).set_index('sample_id')
    cal = json.load(open(path(CAL, raw)))
    half = float(cal['rt']['tolerance_half_width_s'])
    run_lo, run_hi = feat.rt_start.min(), feat.rt_end.max()
    pop = pop.join(feat.drop(columns=['precursor_id', 'is_decoy', 'rt_apex']), on='sample_id')
    pop['mapped_rt'] = pop.rt_apex - pop.apex_rt_error_seconds
    pop['edge'] = ((pop.mapped_rt - half < run_lo) | (pop.mapped_rt + half > run_hi)).astype(int)
    pop['nwin'] = pop.groupby('precursor_id').sample_id.transform('size')
    # Class of every precursor: 0 human target, 1 entrapment target, 2 decoy of human, 3 decoy of entrapment.
    tgt_of_pair = pop[pop.is_decoy == 0].drop_duplicates('pair_key').set_index('pair_key').precursor_id
    pop['partner_e'] = pop.pair_key.map(tgt_of_pair).map(ent)
    pop['e'] = pop.precursor_id.map(ent)
    pop['cls'] = np.where(pop.is_decoy == 0, pop.e, 2 + pop.partner_e.fillna(-9))
    print(f'==== {TAG} {raw} {ARM}: population {len(pop)} windows, {pop.precursor_id.nunique()} precursors; '
          f'RT window half width {half:.1f} s, run {run_lo:.0f}-{run_hi:.0f} s')

    # 1. Candidate windows per precursor by class (no scores involved).
    per = pop.drop_duplicates('precursor_id')
    names = {0: 'human', 1: 'entrap', 2: 'dec(hum)', 3: 'dec(ent)'}
    print('1. windows per precursor (mean, p50, p90, share with 5+), edge-window share, max support >= 3 share')
    for k_, g in per.groupby('cls'):
        if k_ not in names:
            continue
        ms = pop[pop.cls == k_].groupby('precursor_id').transition_support.max()
        print(f'   {names[k_]:9s} n {len(g):6d}  windows {g.nwin.mean():.2f} {g.nwin.median():.0f} {g.nwin.quantile(0.9):.0f} '
              f'{(g.nwin >= 5).mean():.3f}  edge {pop[pop.cls == k_].edge.mean():.3f}  support>=3 {(ms >= 3).mean():.3f}')

    w, best = tdc(pop, s)
    best = best.join(pre[['charge', 'mz']], on='precursor_id').sort_values('s', ascending=False).reset_index(drop=True)
    # Entrapment targets and decoys among the accepted pair winners.
    print('2. overall: E = entrapment target winners, D = decoy pair winners, both at the target cut')
    for q in (0.01, 0.02, 0.05, 0.10, 0.20):
        acc = best[best.q <= q]
        cut = acc.s.min() if len(acc) else np.inf
        top = best[best.s >= cut]
        E = int((top.cls == 1).sum()); D = int((top.is_decoy == 1).sum()); N = int((top.is_decoy == 0).sum())
        print(f'   q <= {q:.2f} cut {cut:7.2f}:' + fmt(E, D, N))

    # 3. Strata. Use the q <= 0.10 cut for counts (more events), report q <= 0.01 too where possible.
    acc = best[best.q <= 0.10]; cut10 = acc.s.min()
    acc1 = best[best.q <= 0.01]; cut1 = acc1.s.min()
    top10 = best[best.s >= cut10]; top1 = best[best.s >= cut1]

    def strata(name, f):
        print(f'3. by {name} (q <= 0.10 cut | q <= 0.01 cut)')
        for lab, g in top10.groupby(f(top10)):
            g1 = top1[f(top1) == lab]
            print(f'   {str(lab):>14s}' + fmt(int((g.cls == 1).sum()), int((g.is_decoy == 1).sum()), int((g.is_decoy == 0).sum())) +
                  ' |' + fmt(int((g1.cls == 1).sum()), int((g1.is_decoy == 1).sum())))

    strata('edge window', lambda d: d.edge)
    strata('windows of the precursor', lambda d: pd.cut(d.nwin, [0, 1, 2, 4, 8, 1000]).astype(str))
    if 'peak_rank' in best:
        strata('picker peak rank of the winner', lambda d: d.peak_rank.clip(upper=3))
    strata('transition support', lambda d: d.transition_support.clip(upper=6))
    strata('|apex RT error| / half window', lambda d: pd.cut((d.apex_rt_error_seconds.abs() / half), [-1, 0.1, 0.25, 0.5, 0.75, 2]).astype(str))
    strata('charge', lambda d: d.charge.clip(upper=4))
    strata('precursor m/z', lambda d: pd.cut(d.mz, [0, 500, 650, 800, 3000]).astype(str))
    strata('mapped RT tertile', lambda d: pd.cut(d.mapped_rt, [-1e9, run_lo + (run_hi - run_lo) / 3, run_lo + 2 * (run_hi - run_lo) / 3, 1e9],
                                                 labels=['early', 'mid', 'late']).astype(str))
    strata('fold', lambda d: d.fold)
    strata('decoy partner class (targets: own class)', lambda d: d.cls)
    print()
