#!/usr/bin/env python3
"""NEAR audit (diagnosis): are the non-borrowing excess entrapments near-copies of human sequence that is not in the library?
For every entrapment peptide accepted (tool q <= QA) in any given run, plus a fixed random sample of library entrapment
peptides (control, seed 0), the minimum Hamming distance (I=L) to any same-length window of the human FASTA, capped at 3
(classes 0, 1, 2, 3+). Exact pigeonhole search: a window within Hamming <= 2 matches one of 3 segments exactly.
Then per set: accepted entrapments by class at q <= 0.01, and the bound when classes <= c are treated as ambiguous
(dropped from N, E and from the library entrapment count used for r; library class fractions from the control sample).
Usage: near_audit.py <lib.pqp> <human.fasta> <out prefix> <set>=<glob> ...   env RELABEL as elsewhere (applied first)"""
import sqlite3, sys, os, glob, random, numpy as np, pandas as pd
lib, fa, out = sys.argv[1:4]; specs = sys.argv[4:]
QA, NCTRL, KMAX = 0.05, 20000, 12
AA = 'ACDEFGHIKLMNPQRSTVWY'; code = {a: i + 1 for i, a in enumerate(AA)}
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[(P.d == 0) & P.useq.isin(rl), 'e'] = 0
T = P[P.d == 0]; Pi = P.set_index('id')
# human proteome as one int8 array, 0 = separator, I mapped to L
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
Hc = lut[H]
idx = {}
def kidx(k):
    if k not in idx:
        n = len(Hc) - k + 1; v = np.zeros(n, np.int64); bad = np.zeros(n, bool)
        for j in range(k): v = v * 32 + Hc[j:j + n]; bad |= Hc[j:j + n] == 0
        pos = np.where(~bad)[0]; v = v[pos]; o = np.argsort(v, kind='stable'); idx[k] = (v[o], pos[o])
    return idx[k]
def enc(p): return np.array([code['L'] if a == 'I' else code.get(a, 31) for a in p], np.int64)
def mind(p):
    q = enc(p); L = len(q); best = 3
    cuts = [0, L // 3, (2 * L) // 3, L]
    for g in range(3):
        a, b = cuts[g], cuts[g + 1]; k = min(b - a, KMAX)
        if k < 1: continue
        key = 0
        for x in q[a:a + k]: key = key * 32 + int(x)
        v, pos = kidx(k); lo, hi = np.searchsorted(v, key), np.searchsorted(v, key, side='right')
        st = pos[lo:hi] - a; st = st[(st >= 0) & (st + L <= len(Hc))]
        if len(st) == 0: continue
        W = Hc[st[:, None] + np.arange(L)[None, :]]
        ok = ~(W == 0).any(axis=1)
        if not ok.any(): continue
        d = (W[ok] != q[None, :]).sum(axis=1).min(); best = min(best, int(d))
        if best == 0: return 0
    return best
# accepted entrapment peptides across runs
runs = []
for spec in specs:
    name, g = spec.split('=', 1)
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        t = t[(t.is_decoy == 0)].join(Pi[['e', 'useq']], on='precursor_id')
        runs.append((name, f, t))
need = sorted(set().union(*[set(t[(t.e == 1) & (t.q_value <= QA)].useq) for _, _, t in runs]))
ent_lib = sorted(set(T[T.e == 1].useq)); rng = random.Random(0); ctrl = sorted(rng.sample(ent_lib, min(NCTRL, len(ent_lib))))
print(f'human residues {int((H != 35).sum())}; accepted entrapment peptides (q <= {QA}) {len(need)}; control sample {len(ctrl)} of {len(ent_lib)}', flush=True)
D = {}
for i, p in enumerate(sorted(set(need) | set(ctrl))):
    D[p] = mind(p)
    if i % 5000 == 0: print(f'  {i} done', flush=True)
pd.DataFrame({'seq': list(D), 'hd': list(D.values())}).to_csv(out + '_hamming.tsv', sep='\t', index=False)
cf = pd.Series([D[p] for p in ctrl]).value_counts(normalize=True).sort_index()
print('library entrapment peptides by min Hamming to human (control sample): ' + ', '.join(f'{k}: {v:.4f}' for k, v in cf.items()), flush=True)
nE = (T.e == 1).sum(); nH = (T.e == 0).sum()
for name in dict.fromkeys(n for n, _, _ in runs):
    rr = [t for n, _, t in runs if n == name]
    acc = pd.concat([t[t.q_value <= 0.01] for t in rr]); acc = acc.assign(hd=acc.useq.map(D))
    N, E = len(acc), int(acc.e.sum()); r = nE / nH
    byc = acc[acc.e == 1].hd.value_counts().sort_index()
    line = f'{name}: runs {len(rr)} N {N} E {E} lower {E / N:.4f} combined {E * (1 + 1 / r) / N:.4f} | accepted E by class ' + \
           ', '.join(f'{k}: {int(byc.get(k, 0))} ({byc.get(k, 0) / max(E, 1):.3f} vs lib {cf.get(k, 0):.3f})' for k in range(4))
    print(line, flush=True)
    for cmax in (1, 2):
        amb = acc[(acc.e == 1) & (acc.hd <= cmax)]
        N2, E2 = N - len(amb), E - len(amb)
        r2 = nE * (1 - cf[cf.index <= cmax].sum()) / nH
        print(f'    classes <= {cmax} ambiguous: N {N2} E {E2} lower {E2 / N2:.4f} combined {E2 * (1 + 1 / r2) / N2:.4f} (r {r2:.4f})', flush=True)
