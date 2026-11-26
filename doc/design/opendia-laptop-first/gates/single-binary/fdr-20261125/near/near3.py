#!/usr/bin/env python3
"""NEAR3 (diagnosis): mechanism of the near-human (Hamming-1) entrapment excess.
For accepted (q <= 0.01) Hamming-1 entrapments in each run: is the human neighbour window (I=L) itself a library target
peptide (any charge), is it accepted in the same run, and how far apart are the apexes? Also, for the GBM question, the
H class (0/1/2/3+) of entrapments promoted by the GBM (G0 accepted, same-seed NET not) vs kept, and of decoys promoted.
Usage: near3.py <lib.pqp> <human.fasta> <hamming.tsv> <label> <net glob> <g0 glob>   env RELABEL"""
import sqlite3, sys, os, glob, numpy as np, pandas as pd
lib, fa, hmt, label, gn, gg = sys.argv[1:7]
KMAX = 12; AA = 'ACDEFGHIKLMNPQRSTVWY'; code = {a: i + 1 for i, a in enumerate(AA)}; inv = {i + 1: a for i, a in enumerate(AA)}
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[(P.d == 0) & P.useq.isin(rl), 'e'] = 0
Pi = P.set_index('id'); T = P[P.d == 0]
hum_pep = {}
for pid, s in zip(T[T.e == 0].id, T[T.e == 0].useq): hum_pep.setdefault(s.replace('I', 'L'), []).append(pid)
HD = dict(pd.read_csv(hmt, sep='\t').values)
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
def neighbours(p):
    q = np.array([code['L'] if a == 'I' else code.get(a, 31) for a in p], np.int64); L = len(q); out = set()
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
        for w in W[(W != q[None, :]).sum(axis=1) == 1]: out.add(''.join(inv[int(x)] for x in w))
    return out
def load(f):
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    return t.join(Pi[['e', 'useq', 'key']])
fn, fg = sorted(glob.glob(gn)), sorted(glob.glob(gg))
stats = dict(n=0, in_lib=0, acc=0, scored=0); dts = []
prom = {}; keep = {}; dprom = 0; dkeep = 0
for a, b in zip(fn, fg):
    tn, tg = load(a), load(b)
    accg = tg[(tg.is_decoy == 0) & (tg.q_value <= 0.01)]; accn = tn[(tn.is_decoy == 0) & (tn.q_value <= 0.01)]
    acc_h = set(accg[accg.e == 0].index)
    for pid, row in accg[accg.e == 1].iterrows():
        h = HD.get(row.useq, 3)
        (keep if pid in accn.index else prom)[h] = (keep if pid in accn.index else prom).get(h, 0) + 1
        if h != 1: continue
        stats['n'] += 1
        nb = [x for s_ in neighbours(row.useq) for x in hum_pep.get(s_, [])]
        if nb: stats['in_lib'] += 1
        sc = [x for x in nb if x in tg.index]
        if sc: stats['scored'] += 1
        ac = [x for x in nb if x in acc_h]
        if ac:
            stats['acc'] += 1; dts.append(min(abs(tg.rt_apex[x] - row.rt_apex) for x in ac))
    decg = set(tg.index[(tg.is_decoy == 1) & (tg.q_value <= 0.01)]); decn = set(tn.index[(tn.is_decoy == 1) & (tn.q_value <= 0.01)])
    dprom += len(decg - decn); dkeep += len(decg & decn)
dts = np.array(dts)
print(f'== {label}: runs {len(fn)}')
print(f'  accepted Hamming-1 entrapments {stats["n"]}: human neighbour is a library target {stats["in_lib"]} ({stats["in_lib"] / max(stats["n"], 1):.3f}), '
      f'scored in run {stats["scored"]}, accepted in run {stats["acc"]} ({stats["acc"] / max(stats["n"], 1):.3f})')
if len(dts): print(f'  apex gap to the accepted neighbour: median {np.median(dts):.1f} s, <= 3 s {np.mean(dts <= 3):.3f}, <= 10 s {np.mean(dts <= 10):.3f}')
tp, tk = sum(prom.values()), sum(keep.values())
print('  accepted entrapments by Hamming class, GBM-promoted vs kept: ' + ' | '.join(
    f'H{h}: promoted {prom.get(h, 0)} ({prom.get(h, 0) / max(tp, 1):.3f}) kept {keep.get(h, 0)} ({keep.get(h, 0) / max(tk, 1):.3f})' for h in (1, 2, 3)))
print(f'  counted decoys promoted {dprom}, kept {dkeep}')
