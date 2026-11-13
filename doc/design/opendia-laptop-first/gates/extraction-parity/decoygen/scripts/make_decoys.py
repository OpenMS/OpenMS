"""PLAN_DECOYGEN step 1 (library only, label-free): decoy peptide sequences for each arm, from the PAIRLIB targets.

Arms (one decoy per target modified peptide; every charge of that peptide gets the same decoy sequence):
  A  shuffle: residues 1..n-1 permuted, C-terminal residue kept (same composition and precursor mass). Up to 100 seeded
     permutations (seed = FNV-1a 64 of the target sequence); the one with the fewest singly charged b/y ions within
     0.01 Da of the target's own b/y ions is kept (y1 and b(n-1) are equal by construction and not counted); first
     permutation with 0 stops the search. (DIA-NN's aim that shuffling changes every fragment mass, --dg-min-shuffle.)
  B  mutate: DIA-NN's rule. Residues 2 and n-1 (adjacent to the termini) are mutated GAVLIFMPWSCTYHKRQEND ->
     LLLVVLLLLTSSSSLLNDQE. A modified residue (C carbamidomethyl, M oxidation) is skipped for the adjacent residue farther
     from the terminus, then the terminal residue. The mutated residue carries no modification. Precursor mass changes.
  C  pseudo-reverse without the K/R swap: residues 1..n-1 reversed, C-terminal residue kept (same mass).
  D  pseudo-reverse with the K/R swap (the current library rule, OpenSwathDecoyGenerator default). Its sequences must equal
     the current library decoys; this is checked.
Modifications move with their residue. A decoy whose unmodified sequence (I = L) equals any library target sequence, or
another target's decoy in the same arm, is dropped (no decoy for that target).
Outputs (decoygen-20261113/lib/): decoy_map.tsv (arm, target ModifiedPeptideSequence, decoy modified sequence),
decoy_peptides.fasta (union of decoy unmodified sequences over arms, for prediction), make_decoys.txt (counts, checks)."""
import re, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
SRC = B + '/pairlib-20261113/lib/pairlib.decoys.tsv'
OUT = B + '/decoygen-20261113/lib'
H2O, TOL, TRIES = 18.010564684, 0.01, 100
AA = dict(G=57.021464, A=71.037114, S=87.032028, P=97.052764, V=99.068414, T=101.047679, C=103.009185, L=113.084064,
          I=113.084064, N=114.042927, D=115.026943, Q=128.058578, K=128.094963, E=129.042593, M=131.040485, H=137.058912,
          F=147.068414, R=156.101111, Y=163.06332, W=186.079313)
MOD = {'UniMod:4': 57.021464, 'UniMod:35': 15.994915}
MUT = dict(zip('GAVLIFMPWSCTYHKRQEND', 'LLLVVLLLLTSSSSLLNDQE'))


def parse(ms):
    """Modified sequence -> list of (residue, modification or '')."""
    return [(a, m or '') for a, m in re.findall(r'([A-Z])(?:\((UniMod:\d+)\))?', ms)]


def fmt(tok):
    return ''.join(a + (f'({m})' if m else '') for a, m in tok)


def masses(tok):
    return np.array([AA[a] + (MOD[m] if m else 0.0) for a, m in tok])


def by_ions(res):
    """Singly charged neutral b/y ladder masses for ordinals 1..n-1 (b: prefix sums, y: suffix sums + water)."""
    c = np.cumsum(res)
    b = c[:-1]; y = (c[-1] - c[:-1]) + H2O  # y_(n-k) for k = 1..n-1
    return b, y


def shared(tt, dt):
    """Decoy b/y ions within TOL of any target b/y ion, excluding y1 and b(n-1) (equal by construction)."""
    tb, ty = by_ions(masses(tt)); db, dy = by_ions(masses(dt))
    t = np.sort(np.r_[tb, ty])
    # b ladder ends with b(n-1) and the y ladder (y(n-1)..y1) ends with y1: both are dropped.
    d = np.r_[db[:-1], dy[:-1]]
    i = np.searchsorted(t, d)
    lo = np.abs(d - t[np.clip(i - 1, 0, len(t) - 1)]); hi = np.abs(t[np.clip(i, 0, len(t) - 1)] - d)
    return int((np.minimum(lo, hi) <= TOL).sum())


def fnv(s):
    h = 0xcbf29ce484222325
    for ch in s.encode():
        h = ((h ^ ch) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def arm_A(tok):
    rng = np.random.default_rng(fnv(''.join(a for a, _ in tok)))
    best = None
    for _ in range(TRIES):
        p = rng.permutation(len(tok) - 1)
        d = [tok[i] for i in p] + [tok[-1]]
        if d == tok:
            continue
        s = shared(tok, d)
        if best is None or s < best[0]:
            best = (s, d)
        if s == 0:
            break
    return best[1] if best else None


def arm_B(tok):
    d = list(tok); n = len(d)
    for order in ((1, 2, 0), (n - 2, n - 3, n - 1)):
        for i in order:
            if not d[i][1]:
                break
        else:
            i = order[2]
        d[i] = (MUT[d[i][0]], '')
    return d


def arm_C(tok):
    return tok[:-1][::-1] + [tok[-1]]


def arm_D(tok):
    a, m = tok[-1]
    return tok[:-1][::-1] + [({'K': 'R', 'R': 'K'}.get(a, a), m)]


cols = ['PeptideSequence', 'ModifiedPeptideSequence', 'Decoy']
lib = pd.read_csv(SRC, sep='\t', usecols=cols + ['TransitionGroupId'], dtype=str)
tg = lib[lib.Decoy == '0'].drop_duplicates('ModifiedPeptideSequence')
old = lib[lib.Decoy == '1'].drop_duplicates('TransitionGroupId')
target_seqs = set(tg.PeptideSequence.str.replace('I', 'L'))
log = open(OUT + '/make_decoys.txt', 'w')
pr = lambda *a: (print(*a), print(*a, file=log))
pr('target modified peptides', len(tg), 'unmodified', tg.PeptideSequence.nunique())
rows = []
for arm, fn in (('A', arm_A), ('B', arm_B), ('C', arm_C), ('D', arm_D)):
    seen = set(); drop_t = drop_d = 0; ov = []
    # One decoy sequence per unmodified target sequence keeps oxidised and plain forms on the same permutation (arm A).
    for ms in tg.ModifiedPeptideSequence:
        tok = parse(ms); d = fn(tok)
        if d is None:
            drop_t += 1; continue
        u = ''.join(a for a, _ in d).replace('I', 'L')
        if u in target_seqs:
            drop_t += 1; continue
        key = fmt(d)
        if key in seen:
            drop_d += 1; continue
        seen.add(key)
        rows.append((arm, ms, key, ''.join(a for a, _ in d)))
        if len(ov) < 20000: ov.append(shared(tok, d))
    ov = np.array(ov)
    pr(f'arm {arm}: decoys {len(seen)}, dropped (equals a target) {drop_t}, dropped (duplicate decoy) {drop_d}, '
       f'shared b/y ions with own target (first 20000): mean {ov.mean():.2f}, share 0 {np.mean(ov == 0):.3f}, >=3 {np.mean(ov >= 3):.3f}')
m = pd.DataFrame(rows, columns=['arm', 'target_ms', 'decoy_ms', 'decoy_seq'])
# Check arm D against the current library decoys: decoy TransitionGroupId = DECOY_<target ModifiedPeptideSequence>/<z>.
old['target_ms'] = old.TransitionGroupId.str.replace('^DECOY_', '', regex=True).str.rsplit('/', n=1).str[0]
od = old.drop_duplicates('target_ms').set_index('target_ms').ModifiedPeptideSequence
d = m[m.arm == 'D'].set_index('target_ms').decoy_ms
common = d.index.intersection(od.index)
pr('arm D vs current library decoys: compared', len(common), 'identical', int((d.loc[common] == od.loc[common]).sum()))
m.to_csv(OUT + '/decoy_map.tsv', sep='\t', index=False)
u = sorted(set(m.decoy_seq))
with open(OUT + '/decoy_peptides.fasta', 'w') as f:
    for s in u:
        f.write(f'>{s}\n{s}\n')
pr('prediction FASTA entries', len(u))
# Reference (library only): shared b/y ions between each FDRBench entrapment and its human parent (same rule as above).
prot = lib[lib.Decoy == '0'].drop_duplicates('ModifiedPeptideSequence')
names, seqs = [], []
for line in open(B + '/pairlib-20261113/lib/pairlib.fasta'):
    line = line.strip()
    if line.startswith('>'):
        names.append(line[1:].split()[0]); seqs.append('')
    else:
        seqs[-1] += line
ov = [shared(parse(seqs[i]), parse(seqs[i + 1])) for i in range(0, min(len(seqs), 40000), 2)
      if '_p_target' in names[i + 1] and set(seqs[i]) <= set(AA)]
ov = np.array(ov)
pr(f'reference, entrapment vs parent (first 20000 pairs, unmodified): shared b/y mean {ov.mean():.2f}, share 0 {np.mean(ov == 0):.3f}, '
   f'>=3 {np.mean(ov >= 3):.3f}')
