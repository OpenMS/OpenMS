#!/usr/bin/env python3
"""STRICT list (library only, no run data): every entrapment-labelled library target peptide (after the step A relabel)
with minimum Hamming distance (I=L) <= 2 to any same-length window of the human FASTA. These are 'near-homolog'
entrapments: a real human or human-variant peptide can produce almost the same fragments. A strict homolog-free
evaluation counts them as neither human nor entrapment (removed from N, E and from r).
Output: <out>.tsv (seq, hd for hd <= 2) and a summary with the library counts.
Usage: strict_list.py <lib.pqp> <human.fasta> <out prefix>   env RELABEL"""
# env HUMAN=1: classify the human (non-entrapment) library peptides instead, by min NONZERO Hamming distance to a
# human window (nearest distinct human sequence), so entrapment and human neighbour density can be compared.
import sqlite3, sys, os, numpy as np, pandas as pd
lib, fa, out = sys.argv[1:4]
HUMAN = os.environ.get('HUMAN') == '1'
# env DECOYS=<n>: classify a fixed random sample (seed 0) of n decoy peptides instead (min Hamming to a human window).
DECOYS = int(os.environ.get('DECOYS', '0'))
KMAX = 12; AA = 'ACDEFGHIKLMNPQRSTVWY'; code = {a: i + 1 for i, a in enumerate(AA)}
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE useq, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PEPTIDE pe join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  where pe.DECOY=0 group by pe.ID""", c)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[P.useq.isin(rl), 'e'] = 0
P = P.groupby('useq').e.min().reset_index()
ent = sorted(P[P.e == (0 if HUMAN else 1)].useq); nh = int((P.e == 0).sum())
if DECOYS:
    dd_ = pd.read_sql("select distinct UNMODIFIED_SEQUENCE useq from PEPTIDE where DECOY=1", c).useq
    ent = sorted(dd_.sample(min(DECOYS, len(dd_)), random_state=0))
# env SEQS=<file>: classify the peptide sequences listed one per line instead.
if os.environ.get('SEQS'): ent = sorted(set(open(os.environ['SEQS']).read().split()))
seqs, s = [], []
for l in open(fa):
    if l.startswith('>'):
        if s: seqs.append(''.join(s)); s = []
    else: s.append(l.strip())
if s: seqs.append(''.join(s))
H = np.frombuffer(('#'.join(seqs) + '#').replace('I', 'L').encode(), dtype=np.uint8)
lut = np.zeros(256, np.int64)
for a, i in code.items(): lut[ord(a)] = i
lut[ord('I')] = code['L']
Hc = lut[H]; idx = {}
def kidx(k):
    if k not in idx:
        n = len(Hc) - k + 1; v = np.zeros(n, np.int64); bad = np.zeros(n, bool)
        for j in range(k): v = v * 32 + Hc[j:j + n]; bad |= Hc[j:j + n] == 0
        pos = np.where(~bad)[0]; v = v[pos]; o = np.argsort(v, kind='stable'); idx[k] = (v[o], pos[o])
    return idx[k]
def mind(p):
    q = np.array([code['L'] if a == 'I' else code.get(a, 31) for a in p], np.int64); L = len(q); best = 3
    cuts = [0, L // 3, (2 * L) // 3, L]
    for g in range(3):
        a, b = cuts[g], cuts[g + 1]; k = min(b - a, KMAX)
        if k < 1: continue
        key = 0
        for x in q[a:a + k]: key = key * 32 + int(x)
        v, pos = kidx(k); lo, hi = np.searchsorted(v, key), np.searchsorted(v, key, side='right')
        st = pos[lo:hi] - a; st = st[(st >= 0) & (st + L <= len(Hc))]
        if len(st) == 0: continue
        W = Hc[st[:, None] + np.arange(L)[None, :]]; W = W[~(W == 0).any(axis=1)]
        if len(W) == 0: continue
        dd = (W != q[None, :]).sum(axis=1)
        if HUMAN: dd = dd[dd > 0]
        if len(dd) == 0: continue
        best = min(best, int(dd.min()))
        if best == 0: return 0
    return best
hd = []
for i, p in enumerate(ent):
    hd.append(mind(p))
    if i % 50000 == 0: print(f'  {i}/{len(ent)}', flush=True)
R = pd.DataFrame({'seq': ent, 'hd': hd})
if os.environ.get('SEQS'): R.to_csv(out + '.all.tsv', sep='\t', index=False)
R[R.hd <= 2].to_csv(out + '.tsv', sep='\t', index=False)
vc = R.hd.value_counts().sort_index()
print(f'library: human peptides {nh}, entrapment peptides {len(ent)}; by min Hamming: ' + ', '.join(f'{k}: {v}' for k, v in vc.items()))
print(f'strict homolog-free entrapment peptides (hd >= 3): {int((R.hd >= 3).sum())}; near-homolog list (hd <= 2): {int((R.hd <= 2).sum())}')
