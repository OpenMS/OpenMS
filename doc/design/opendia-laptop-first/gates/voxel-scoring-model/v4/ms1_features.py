"""MS1 isotope features per population window (label-free; from ms1.npy [n, 4, 32] = M+0..M+3 over apex +- 15 s, and
xic.npy [n, 6, 32]). Theoretical isotope envelope from the peptide's elemental composition (averagine-free: exact C/H/N/O/S
counts of the modified sequence, binomial-convolution approximation). Usage: ms1_features.py <raw01|raw02>.
Output: v2/<raw>/ms1feat.npy [n, 6] float32, columns MS1_FEATS."""
import sys, sqlite3, re, numpy as np, pandas as pd
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
MS1_FEATS = ['ms1_log_m0', 'ms1_m0_xic_corr', 'ms1_iso_cos', 'ms1_iso_corr_m1', 'ms1_apex_shift', 'ms1_m0_ratio']
# Residue compositions (C, H, N, O, S) of the residue in a chain (minus H2O).
AA = {'G': (2, 3, 1, 1, 0), 'A': (3, 5, 1, 1, 0), 'S': (3, 5, 1, 2, 0), 'P': (5, 7, 1, 1, 0), 'V': (5, 9, 1, 1, 0), 'T': (4, 7, 1, 2, 0),
      'C': (3, 5, 1, 1, 1), 'L': (6, 11, 1, 1, 0), 'I': (6, 11, 1, 1, 0), 'N': (4, 6, 2, 2, 0), 'D': (4, 5, 1, 3, 0), 'Q': (5, 8, 2, 2, 0),
      'K': (6, 12, 2, 1, 0), 'E': (5, 7, 1, 3, 0), 'M': (5, 9, 1, 1, 1), 'H': (6, 7, 3, 1, 0), 'F': (9, 9, 1, 1, 0), 'R': (6, 12, 4, 1, 0),
      'Y': (9, 9, 1, 2, 0), 'W': (11, 10, 2, 1, 0)}
MOD = {'UniMod:4': (2, 3, 1, 1, 0), 'UniMod:35': (0, 0, 0, 1, 0)}
# Heavy-isotope abundances: (element, [p(+1), p(+2)]).
ISO = [(0, [0.0107, 0.0]), (1, [0.000115, 0.0]), (2, [0.00364, 0.0]), (3, [0.00038, 0.00205]), (4, [0.0075, 0.0425])]

def envelope(comp, n=4):
    p = np.zeros(n); p[0] = 1.0
    for e, ab in ISO:
        k = comp[e]
        if k == 0: continue
        dist = np.zeros(n); dist[0] = 1 - sum(ab)
        for j, a in enumerate(ab, 1):
            if j < n: dist[j] = a
        for _ in range(int(k)): p = np.convolve(p, dist)[:n]
    return p / p.sum()

def composition(ms):
    comp = np.array([0, 2, 0, 1, 0])  # H2O
    for m in re.finditer(r'\((UniMod:\d+)\)', ms): comp += MOD.get(m.group(1), (0, 0, 0, 0, 0))
    for r in re.sub(r'\(.*?\)', '', ms): comp += AA[r]
    return comp

def corr(a, b):
    a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1) + 1e-12)

c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
seq = pd.read_sql("select m.PRECURSOR_ID precursor_id, pe.MODIFIED_SEQUENCE ms from PRECURSOR_PEPTIDE_MAPPING m join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c).set_index('precursor_id').ms
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
uniq = pop.precursor_id.unique(); env = {p: envelope(composition(seq[p])) for p in uniq}
E = np.stack([env[p] for p in pop.precursor_id.values]).astype(np.float32)
ms1 = np.load(f'{O}/{RAW}/ms1.npy', mmap_mode='r'); xic = np.load(f'{O}/{RAW}/xic.npy', mmap_mode='r')
n = len(pop); out = np.zeros((n, len(MS1_FEATS)), np.float32); C0, C1 = 12, 20  # central bins: apex +- 4.7 s
for b in range(0, n, 200000):
    m = np.asarray(ms1[b:b + 200000], np.float64); x = np.asarray(xic[b:b + 200000], np.float64).sum(1); e = E[b:b + 200000]
    m0 = m[:, 0]; area = m[:, :, C0:C1].sum(2)
    has = m0.max(1) > 0
    f = np.zeros((len(m), len(MS1_FEATS)))
    f[:, 0] = np.log10(area[:, 0] + 1)
    f[:, 1] = np.where(has, corr(m0, x), 0)
    obs = area / np.maximum(area.sum(1, keepdims=True), 1e-12)
    f[:, 2] = np.where(has, (obs * e).sum(1) / (np.linalg.norm(obs, axis=1) * np.linalg.norm(e, axis=1) + 1e-12), 0)
    f[:, 3] = np.where(has, corr(m0, m[:, 1]), 0)
    f[:, 4] = np.where(has, np.abs(m0.argmax(1) - x.argmax(1)), 16)
    f[:, 5] = np.where(has, np.log10((area[:, 0] + 1) / (m0.sum(1) + 1)), 0)  # share of M+0 signal near the apex
    out[b:b + 200000] = np.nan_to_num(f).astype(np.float32)
    print(RAW, b, '/', n, flush=True)
np.save(f'{O}/{RAW}/ms1feat.npy', out)
print(RAW, 'ms1 present', float((out[:, 0] > 0).mean()), 'feature medians', dict(zip(MS1_FEATS, np.median(out, 0).round(3))))
