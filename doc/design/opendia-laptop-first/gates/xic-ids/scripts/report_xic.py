"""PLAN_XIC_IDS reporting (evaluation only, reads entrapment labels). Compares accepted precursor sets of D11FWT across
pipeline builds whose populations differ, so the comparison is on precursors, not windows.
Per raw and build: precursor winners, pooled paired TDC (as evaluate_v2.py / gains.py), q <= 0.01 accepted targets
(human / entrapment / paired decoys above the cut), N at FDP_c 1 %, and gained / lost precursors against the first build.
Raw01 only: overlap of each accepted set with DIA-NN 2.3.1 on our library at q <= 0.01 (diann-ourlib-20261112).
Usage: report_xic.py <label>=<root>:<arm> [...]; the first entry is the reference."""
import sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
dn = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'e', 'Q.Value'])
dn = dn[dn['Q.Value'] <= 0.01].drop_duplicates('precursor_id')
DN = set(dn.precursor_id)

def tdc(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return w, best

def fdp_set(w, t=0.01):
    tw = w[w.is_decoy == 0].sort_values('s', ascending=False, kind='stable')
    e = tw.precursor_id.map(ent).fillna(0).values
    f = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(f <= t)[0]
    return set(tw.precursor_id.values[: int(ok.max() + 1) if len(ok) else 0])

builds = []
for spec in sys.argv[1:]:
    label, rest = spec.split('=', 1); root, arm = rest.rsplit(':', 1); builds.append((label, root, arm))
E = lambda S: sum(int(ent.get(p, 0)) for p in S)
for raw in ('raw01', 'raw02'):
    sets = {}
    for label, root, arm in builds:
        pop = pd.read_csv(f'{root}/{raw}/population.tsv', sep='\t')
        w, best = tdc(pop, np.load(f'{root}/{raw}/scores_{arm}.npy'))
        k = best[best.q <= 0.01]; acc = set(k[k.is_decoy == 0].precursor_id); fs = fdp_set(w)
        sets[label] = (acc, fs)
        line = (f'{raw} {label:8s} population windows {len(pop)} precursors {pop.precursor_id.nunique()} | q<=0.01 targets {len(acc)} '
                f'human {len(acc) - E(acc)} entrap {E(acc)} paired decoys {int((k.is_decoy == 1).sum())} | N@FDP_c1% {len(fs)} entrap {E(fs)}')
        if raw == 'raw01':
            inp = DN & set(pop.precursor_id)
            line += (f' | DIA-NN q<=1% ({len(DN)}): in population {len(inp)}, accepted q<=0.01 {len(acc & DN)} '
                     f'(entrap {E(acc & DN)}), accepted not DIA-NN {len(acc - DN)} (entrap {E(acc - DN)})')
        print(line, flush=True)
    ref_label = builds[0][0]; (ra, rf) = sets[ref_label]
    for label, _, _ in builds[1:]:
        a, f = sets[label]
        for name, x, y in (('q<=0.01', ra, a), ('FDP_c1%', rf, f)):
            g, l = y - x, x - y
            print(f'{raw} {label} vs {ref_label} {name}: gained {len(g)} (human {len(g) - E(g)}, entrap {E(g)}'
                  + (f', DIA-NN {len(g & DN)}' if raw == 'raw01' else '') + f') | lost {len(l)} (human {len(l) - E(l)}, entrap {E(l)}'
                  + (f', DIA-NN {len(l & DN)}' if raw == 'raw01' else '') + ')', flush=True)
