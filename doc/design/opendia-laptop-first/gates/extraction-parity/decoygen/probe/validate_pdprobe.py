"""Validation: the in-tree OpenMS PeptDeep inference (pdprobe, built against OpenMS-r1-full-ninja-build) against DIALibGen 0.13.1.
Samples precursors stratified by charge (1-5), length bin and modification (none / CAM / Oxidation), from
  - decoys_p6.parquet (decoys, NCE 23, fragment m/z 100-1800, top 6) and
  - pred12.parquet (PAIRLIB targets, NCE 25, fragment m/z 100-2000, top 12),
predicts them with pdprobe (Lumos instrument slot, as DIALibGen's Astral) and compares:
  RT |difference|; per precursor, DIALibGen's fragment intensities against pdprobe's at the same annotations (max |diff|
  after normalising both to DIALibGen's top fragment); whether pdprobe's top-k inside the same m/z window equals DIALibGen's
  set, and if not whether the difference is a near-tie (intensity gap at the k-th rank < 1e-3 relative).
Output: probe/validate_pdprobe.txt"""
import subprocess, numpy as np, pandas as pd, pyarrow.parquet as pq, re, sys
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
PROBE = '/home/singjc/.jcode/scratch/pdprobe/build/pdprobe'
MODELS = '/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/share/OpenMS/models'
AA = dict(G=57.021464, A=71.037114, S=87.032028, P=97.052764, V=99.068414, T=101.047679, C=103.009185, L=113.084064, I=113.084064,
          N=114.042927, D=115.026943, Q=128.058578, K=128.094963, E=129.042593, M=131.040485, H=137.058912, F=147.068414, R=156.101111,
          Y=163.06332, W=186.079313)
MOD = {'Carbamidomethyl': 57.021464, 'Oxidation': 15.994915}
P, H2O = 1.007276467, 18.010564684
out = open(D + '/decoygen-20261113/probe/validate_pdprobe.txt', 'w')
pr = lambda *a: (print(*a, flush=True), print(*a, file=out, flush=True))


def residues(s):
    return np.array([AA[a] + (MOD[m] if m else 0.0) for a, m in re.findall(r'([A-Z])(?:\(([A-Za-z]+)\))?', s)])


def sample(path, n_per_stratum, seed):
    t = pq.read_table(path, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Product.Mz', 'Relative.Intensity', 'Fragment.Type',
                                     'Fragment.Charge', 'Fragment.Series.Number']).to_pandas()
    t['L'] = t['Modified.Sequence'].str.replace(r'\([A-Za-z]+\)', '', regex=True).str.len()
    t['lbin'] = pd.cut(t.L, [6, 9, 14, 20, 28, 35]).astype(str)
    t['mod'] = np.where(t['Modified.Sequence'].str.contains('Oxidation'), 'ox', np.where(t['Modified.Sequence'].str.contains('Carbamidomethyl'), 'cam', 'none'))
    t['stratum'] = t['Precursor.Charge'].astype(str) + '|' + t.lbin + '|' + t['mod']
    rng = np.random.default_rng(seed)
    keep = np.concatenate([rng.choice(idx, min(len(idx), n_per_stratum), replace=False) for idx in t.groupby('stratum').indices.values()])
    return t.iloc[np.sort(keep)]


def run(name, path, nce, mz_lo, mz_hi, k, per):
    s = sample(path, per, 0).reset_index(drop=True)
    inp = '\n'.join(f'{a} {z}' for a, z in zip(s['Modified.Sequence'], s['Precursor.Charge']))
    res = subprocess.run([PROBE, MODELS, str(nce)], input=inp, capture_output=True, text=True, check=True)
    rows = [l.split('\t') for l in res.stdout.strip().split('\n')]
    assert len(rows) == len(s)
    rt_d, imax, same, tie, diff_examples = [], [], 0, 0, []
    for (seq, z, rt, *v), (_, r) in zip(rows, s.iterrows()):
        assert seq == r['Modified.Sequence']
        z = int(float(z)); v = np.array(v, float); res_m = residues(seq); L = len(res_m); v = v.reshape(L - 1, 8)
        pred = {}
        for i in range(L - 1):
            b, y = res_m[:i + 1].sum(), res_m[i + 1:].sum() + H2O
            for ch, (ion, o, m, fz) in enumerate([('b', i + 1, b, 1), ('b', i + 1, b, 2), ('y', L - 1 - i, y, 1), ('y', L - 1 - i, y, 2)]):
                if fz > min(z, 2): continue  # DIALibGen max_fragment_charge 2 and not above the precursor charge
                mz = (m + fz * P) / fz
                if mz_lo <= mz <= mz_hi: pred[(ion, o, fz)] = v[i, ch]
        keys = [({2: 'b', 5: 'y'}[int(a)], int(o), int(c)) for a, o, c in zip(r['Fragment.Type'], r['Fragment.Series.Number'], r['Fragment.Charge'])]
        ref_i = np.asarray(r['Relative.Intensity'], float)
        p = np.array([pred.get(q, np.nan) for q in keys])
        top = max(pred.values()) if pred else 1.0
        imax.append(np.nanmax(np.abs(p / top - ref_i / ref_i.max())) if np.isfinite(p).any() else np.nan)
        rt_d.append(abs(float(rt) - r['RT']))
        mine = sorted(pred, key=pred.get, reverse=True)
        kk = min(k, len(keys))
        if set(mine[:kk]) == set(keys[:kk]):
            same += 1
        else:
            vals = sorted(pred.values(), reverse=True)
            gap = abs(vals[kk - 1] - vals[kk]) / top if len(vals) > kk else 1.0
            tie += gap < 1e-3
            if len(diff_examples) < 5: diff_examples.append((seq, z, round(gap, 5)))
    rt_d, imax = np.array(rt_d), np.array(imax)
    n = len(s)
    pr(f'== {name}: n {n} (strata {s.groupby(["Precursor.Charge", "lbin", "mod"]).ngroups}), charges {sorted(s["Precursor.Charge"].unique())}, '
       f'length {s.L.min()}-{s.L.max()}, oxidised {np.mean(s["mod"] == "ox"):.2f}')
    pr(f'   RT |diff| max {rt_d.max():.2e} | intensity |diff| at DIALibGen fragments: median {np.nanmedian(imax):.2e}, max {np.nanmax(imax):.2e}')
    pr(f'   top-{k} set identical {same / n:.4f}; differing {n - same}, of which near-ties (gap < 1e-3) {tie}; examples {diff_examples}')


run('decoys, NCE 23, 100-1800 m/z, top 6 (decoys_p6)', D + '/decoygen-20261113/lib/decoys_p6.parquet', 23, 100, 1800, 6, 400)
run('targets, NCE 25, 100-2000 m/z, top 12 (pred12)', D + '/pairlib-20261113/lib/pred12.parquet', 25, 100, 2000, 12, 400)
