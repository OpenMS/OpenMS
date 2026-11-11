"""Peptide sequence featurizer, ported from ReDeeM foundation chemistry (redeem-properties/src/foundation/chemistry.rs,
featurize.rs, commit 9db3927): residue-local heavy-atom graphs (N-CA-C(=O) backbone + side chain), 12 atom features,
bond-order adjacency with self-loops, exact graph PTMs for UniMod:4 on C and UniMod:35 on M (the only mods in this library).
No ReDeeM weights or predictions are used. Output per precursor (padded to L = 64 residues): gid [L] uint8 (residue-graph id,
0 = padding), res [L] uint8 (ReDeeM token ids), rmask [L] bool; and seq_graph_table.npz with atoms [G, A = 20, 12], adj [G, A, A],
amask [G, A] for the distinct residue graphs (residue, mod, N-terminal, C-terminal).
Usage: seq_featurize.py <raw>  ->  v2/<raw>/seq_*.npy, one row per precursor in population order + seq_index.tsv."""
import re, sqlite3, sys, numpy as np, pandas as pd
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
L, A, F = 64, 20, 12
EL = {'C': (6, 12.011, 2.55, 4), 'N': (7, 14.007, 3.04, 3), 'O': (8, 15.999, 3.44, 2), 'S': (16, 32.06, 2.58, 2)}
TOK = {a: i + 1 for i, a in enumerate('ACDEFGHIKLMNPQRSTVWY')}  # ReDeeM residue_token_id order
assert TOK['K'] == 9 and TOK['R'] == 15 and TOK['S'] == 16

def residue_graph(r):
    """Atoms as (element, backbone, aromatic, donor, acceptor); bonds (i, j, order). Same templates as ReDeeM."""
    at = [('N', 1, 0, 1, 1), ('C', 1, 0, 0, 0), ('C', 1, 0, 0, 0), ('O', 1, 0, 0, 1)]; bd = [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 2.0)]
    def add(e, ar=0, do=0, ac=0): at.append((e, 0, ar, do, ac)); return len(at) - 1
    def con(i, j, o=1.0): bd.append((i, j, o))
    if r == 'A': con(1, add('C'))
    elif r == 'V': b = add('C'); c1 = add('C'); c2 = add('C'); con(1, b); con(b, c1); con(b, c2)
    elif r == 'L': b = add('C'); c = add('C'); d1 = add('C'); d2 = add('C'); con(1, b); con(b, c); con(c, d1); con(c, d2)
    elif r == 'I': b = add('C'); c1 = add('C'); c2 = add('C'); d = add('C'); con(1, b); con(b, c1); con(b, c2); con(c1, d)
    elif r == 'S': b = add('C'); o = add('O', 0, 1, 1); con(1, b); con(b, o)
    elif r == 'T': b = add('C'); o = add('O', 0, 1, 1); c = add('C'); con(1, b); con(b, o); con(b, c)
    elif r == 'C': b = add('C'); s = add('S', 0, 1, 1); con(1, b); con(b, s)
    elif r == 'M': b = add('C'); c = add('C'); s = add('S', 0, 0, 1); e = add('C'); con(1, b); con(b, c); con(c, s); con(s, e)
    elif r == 'D': b = add('C'); c = add('C'); o1 = add('O', 0, 0, 1); o2 = add('O', 0, 0, 1); con(1, b); con(b, c); con(c, o1, 2.0); con(c, o2)
    elif r == 'E':
        b = add('C'); c = add('C'); d = add('C'); o1 = add('O', 0, 0, 1); o2 = add('O', 0, 0, 1)
        con(1, b); con(b, c); con(c, d); con(d, o1, 2.0); con(d, o2)
    elif r == 'N': b = add('C'); c = add('C'); o = add('O', 0, 0, 1); n = add('N', 0, 1, 1); con(1, b); con(b, c); con(c, o, 2.0); con(c, n)
    elif r == 'Q':
        b = add('C'); c = add('C'); d = add('C'); o = add('O', 0, 0, 1); n = add('N', 0, 1, 1)
        con(1, b); con(b, c); con(c, d); con(d, o, 2.0); con(d, n)
    elif r == 'K':
        p = 1
        for _ in range(4): c = add('C'); con(p, c); p = c
        con(p, add('N', 0, 1, 1))
    elif r == 'R':
        b = add('C'); c = add('C'); d = add('C'); n1 = add('N', 0, 1, 1); cz = add('C'); n2 = add('N', 0, 1, 1); n3 = add('N', 0, 1, 1)
        for x, y in ((1, b), (b, c), (c, d), (d, n1), (n1, cz)): con(x, y)
        con(cz, n2, 1.5); con(cz, n3, 1.5)
    elif r == 'H':
        b = add('C'); c = add('C', 1); n1 = add('N', 1, 1, 1); c2 = add('C', 1); n2 = add('N', 1, 1, 1); c3 = add('C', 1)
        con(1, b); con(b, c)
        for x, y in ((c, n1), (n1, c2), (c2, n2), (n2, c3), (c3, c)): con(x, y, 1.5)
    elif r in 'FY':
        b = add('C'); ring = [add('C', 1) for _ in range(6)]; con(1, b); con(b, ring[0])
        for k in range(6): con(ring[k], ring[(k + 1) % 6], 1.5)
        if r == 'Y': con(ring[3], add('O', 0, 1, 1))
    elif r == 'W':
        b = add('C'); c1 = add('C', 1); c2 = add('C', 1); n = add('N', 1, 1, 1); c3, c4, c5, c6, c7 = (add('C', 1) for _ in range(5))
        con(1, b); con(b, c1)
        for x, y in ((c1, c2), (c2, n), (n, c3), (c3, c4), (c4, c1), (c4, c5), (c5, c6), (c6, c7), (c7, c3)): con(x, y, 1.5)
    elif r == 'P': b = add('C'); c = add('C'); d = add('C'); con(1, b); con(b, c); con(c, d); con(d, 0)
    elif r != 'G': raise ValueError(r)
    return at, bd

def apply_mod(r, uid, at, bd):
    """Exact ReDeeM graph PTMs used by this library."""
    side_s = next(i for i, a in enumerate(at) if a[0] == 'S' and not a[1])
    if uid == 4 and r == 'C':  # carbamidomethyl: S-CH2-C(=O)-NH2
        base = len(at); at += [('C', 0, 0, 0, 0), ('C', 0, 0, 0, 0), ('O', 0, 0, 0, 1), ('N', 0, 0, 1, 0)]
        bd += [(side_s, base, 1.0), (base, base + 1, 1.0), (base + 1, base + 2, 2.0), (base + 1, base + 3, 1.0)]
    elif uid == 35 and r == 'M':  # sulfoxide
        at.append(('O', 0, 0, 0, 1)); bd.append((side_s, len(at) - 1, 1.0))
    else: raise ValueError((r, uid))

def atom_feats(a, nterm, cterm):
    z, m, en, val = EL[a[0]]
    return [z / 34, m / 80, en / 4, val / 5, 0.0, a[2], a[3], a[4], a[1], 1 - a[1], 0.0, 1.0 if nterm else (-1.0 if cterm else 0.0)]

def parse(modseq):
    """'AC(UniMod:4)DM(UniMod:35)K' -> residues and {index: unimod}."""
    res, mods = [], {}
    for aa, uid in re.findall(r'([A-Z])(?:\(UniMod:(\d+)\))?', modseq):
        res.append(aa)
        if uid: mods[len(res) - 1] = int(uid)
    return res, mods

c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c).set_index('precursor_id').ms
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
pids = np.array(sorted(pop.precursor_id.unique()))
n = len(pids)
# Storage: the graph layers are residue-local, so each residue graph is stored once (table) and indexed per position.
# This is the same computation as dense per-position tensors (implementation note in PLAN_V3_SEQ.notes).
mm = lambda name, shape, dt: np.lib.format.open_memmap(f'{O}/{RAW}/seq_{name}.npy', mode='w+', dtype=dt, shape=shape)
GID, RID, RM = mm('gid', (n, L), np.uint8), mm('res', (n, L), np.uint8), mm('rmask', (n, L), np.bool_)
cache = {}
for k, pid in enumerate(pids):
    res, mods = parse(pre[pid]); assert len(res) <= L
    for i, r in enumerate(res):
        key = (r, mods.get(i), i == 0, i == len(res) - 1)
        if key not in cache:
            at, bd = residue_graph(r); at, bd = list(at), list(bd)
            if mods.get(i): apply_mod(r, mods[i], at, bd)
            assert len(at) <= A, (r, len(at))
            fa = np.zeros((A, F), np.float32); adj = np.zeros((A, A), np.float32); am = np.zeros(A, bool)
            for j, a in enumerate(at): fa[j] = atom_feats(a, key[2], key[3]); adj[j, j] = 1.0; am[j] = True
            for x, y, o in bd: adj[x, y] = o; adj[y, x] = o
            cache[key] = (len(cache) + 1, fa, adj, am)  # graph id 0 = padding
        GID[k, i] = cache[key][0]; RID[k, i] = TOK[r]; RM[k, i] = True
    if k % 50000 == 0: print(RAW, k, '/', n, flush=True)
for a in (GID, RID, RM): a.flush()
g = len(cache) + 1
T_atoms, T_adj, T_am = np.zeros((g, A, F), np.float32), np.zeros((g, A, A), np.float32), np.zeros((g, A), bool)
for gid, fa, adj, am in cache.values(): T_atoms[gid], T_adj[gid], T_am[gid] = fa, adj, am
np.savez(f'{O}/{RAW}/seq_graph_table.npz', atoms=T_atoms, adj=T_adj, amask=T_am, keys=np.array([str(k) for k in cache], dtype=object))
pd.DataFrame({'precursor_id': pids, 'seq_row': np.arange(n)}).to_csv(f'{O}/{RAW}/seq_index.tsv', sep='\t', index=False)
print(RAW, 'precursors', n, 'distinct residue graphs', len(cache), 'max atoms', int(T_am.sum(1).max()))
