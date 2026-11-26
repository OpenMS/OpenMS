#!/usr/bin/env python3
"""NEAR2 (diagnosis): do Hamming-1 near-human entrapments explain the residual excess, and are they modified human peptides?
For entrapment peptides that are either accepted (tool q <= 0.05) or members of a competitive own-decoy pair (best member
q <= 0.05), plus the near_audit control sample: min Hamming (I=L) to a same-length human FASTA window, and for Hamming 1
the substitution (human residue -> entrapment residue) and its mass delta.
Reports per arm (pooled over runs):
  (a) acceptance rate per library member by Hamming class within length bins (length-controlled enrichment);
  (b) P(entrapment beats own decoy) by Hamming class (pairs with best q <= 0.05);
  (c) Hamming-1 substitutions among accepted entrapments vs library control, by mass delta, with known modification
      equivalents marked (deamidation N->D/Q->E +0.984, methyl +14.016 G->A S->T D->E V->L/I A->? N->Q K->?..., oxidation-like
      +15.995 A->S F->Y P->? ...): enrichment ratio;
  (d) bounds after setting aside entrapments whose Hamming-1 substitution equals a common modification (pre-listed below),
      with r corrected for the library fraction set aside.
Usage: near2.py <lib.pqp> <human.fasta> <control hamming tsv> <set>=<glob> ...   env RELABEL"""
import sqlite3, sys, os, glob, numpy as np, pandas as pd
lib, fa, ctrl_tsv = sys.argv[1:4]; specs = sys.argv[4:]
KMAX = 12
MONO = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=103.00919, L=113.08406, I=113.08406,
            N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
            R=156.10111, Y=163.06333, W=186.07931)
# pre-listed common modification equivalents (mass delta of human -> entrapment residue, Da, tolerance 0.02)
MODS = {'deamidation': 0.98402, 'methyl': 14.01565, 'oxidation': 15.99491, 'formyl': 27.99491, 'dimethyl/ethyl': 28.0313,
        'acetyl': 42.01057, 'carbamyl': 43.00581, 'trimethyl': 42.04695, 'K/Q isobaric': 0.03638}
AA = 'ACDEFGHIKLMNPQRSTVWY'; code = {a: i + 1 for i, a in enumerate(AA)}; inv = {i + 1: a for i, a in enumerate(AA)}
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[(P.d == 0) & P.useq.isin(rl), 'e'] = 0
T = P[P.d == 0]; Pi = P.set_index('id')
pairs = T[T.e == 1].set_index('key')[['id', 'useq']].join(P[P.d == 1].set_index('key')[['id']].rename(columns={'id': 'did'}), how='inner')
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
def near(p):
    """(min Hamming capped at 3, set of (human aa, entrapment aa) substitutions at Hamming 1)"""
    q = np.array([code['L'] if a == 'I' else code.get(a, 31) for a in p], np.int64); L = len(q); best = 3; subs = set()
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
        d = (W != q[None, :]).sum(axis=1)
        if d.min() < best: best = int(d.min()); subs = set()
        if best == 1:
            for w in W[d == 1]:
                j = int(np.where(w != q)[0][0]); subs.add((inv[int(w[j])], p[j] if p[j] != 'I' else 'L'))
        if best == 0: return 0, set()
    return best, subs
def delta_class(subs):
    out = set()
    for h, e in subs:
        dm = MONO[e] - MONO[h]
        hit = [m for m, v in MODS.items() if abs(dm - v) <= 0.02]
        out.add(hit[0] if hit else 'other')
    return out
runs = []
for spec in specs:
    name, g = spec.split('=', 1)
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        runs.append((name, t))
ctrl = pd.read_csv(ctrl_tsv, sep='\t')
need = set()
pair_rows = []
for name, t in runs:
    acc = t[(t.is_decoy == 0) & (t.q_value <= 0.05)].index
    need |= set(Pi.useq.reindex(acc)[Pi.e.reindex(acc) == 1])
    x = pairs[pairs.id.isin(t.index) & pairs.did.isin(t.index)]
    qb = np.fmin(t.q_value.reindex(x.id).values, t.q_value.reindex(x.did).values)
    x = x[qb <= 0.05].assign(win=(t.score.reindex(x[qb <= 0.05].id).values > t.score.reindex(x[qb <= 0.05].did).values), arm=name)
    pair_rows.append(x); need |= set(x.useq)
rng = np.random.default_rng(0)
ctrl_seqs = list(ctrl.seq)
need |= set(ctrl_seqs)
print(f'peptides to classify {len(need)}', flush=True)
NB = {}
for i, p in enumerate(sorted(need)):
    NB[p] = near(p)
    if i % 10000 == 0: print(f'  {i}', flush=True)
hd = {p: v[0] for p, v in NB.items()}; mc = {p: delta_class(v[1]) for p, v in NB.items()}
modded = lambda p: hd.get(p) == 1 and bool(mc[p] - {'other'})
C = pd.DataFrame({'seq': ctrl_seqs}); C['hd'] = C.seq.map(hd); C['len'] = C.seq.str.len()
C['lbin'] = pd.cut(C.len, [0, 8, 10, 13, 17, 60]); C['mod'] = C.seq.map(modded)
nE = (T.e == 1).sum(); nH = (T.e == 0).sum(); frac_mod = C['mod'].mean()
print(f'control: Hamming-1 {np.mean(C.hd == 1):.4f}, of which modification-equivalent {C["mod"].sum() / max((C.hd == 1).sum(), 1):.3f} '
      f'(library fraction mod-equivalent {frac_mod:.4f})', flush=True)
PR = pd.concat(pair_rows); PR['hd'] = PR.useq.map(hd); PR['mod'] = PR.useq.map(modded)
for name in dict.fromkeys(n for n, _ in runs):
    rr = [t for n, t in runs if n == name]
    acc = pd.concat([t[(t.is_decoy == 0) & (t.q_value <= 0.01)] for t in rr]).join(Pi[['e', 'useq']])
    N, E = len(acc), int(acc.e.sum()); ae = acc[acc.e == 1].copy()
    ae['hd'] = ae.useq.map(hd); ae['len'] = ae.useq.str.len(); ae['lbin'] = pd.cut(ae.len, [0, 8, 10, 13, 17, 60]); ae['mod'] = ae.useq.map(modded)
    print(f'== {name}: runs {len(rr)} N {N} E {E} lower {E / N:.4f} combined {E * (1 + nH / nE) / N:.4f}', flush=True)
    # (a) length-controlled acceptance rate per library member (relative to class 3+ in the same length bin)
    lines = []
    for lb, g in ae.groupby('lbin', observed=True):
        cg = C[C.lbin == lb]
        rate = {h: (g.hd == h).sum() / max((cg.hd == h).mean(), 1e-9) for h in (1, 2, 3)}
        lines.append(f'{lb}: rel rate H1 {rate[1] / max(rate[3], 1e-9):.2f} H2 {rate[2] / max(rate[3], 1e-9):.2f} (E {len(g)})')
    print('  (a) ' + ' | '.join(lines), flush=True)
    # (b) pair fairness by class
    pr = PR[PR.arm == name]
    print('  (b) P(entrapment beats own decoy): ' + ' | '.join(
        f'H{h}: {pr[pr.hd == h].win.mean():.3f} (n {int((pr.hd == h).sum())})' for h in (1, 2, 3)) +
        f' | H1 mod-equivalent {pr[pr["mod"] == True].win.mean():.3f} (n {int(pr["mod"].sum())}) | H1 other {pr[(pr.hd == 1) & (pr["mod"] == False)].win.mean():.3f}', flush=True)
    # (c) substitution class enrichment among accepted Hamming-1
    a1 = ae[ae.hd == 1]; c1 = C[C.hd == 1]
    def cnt(df):
        z = {}
        for p in df.seq if 'seq' in df else df.useq:
            for m in mc[p]: z[m] = z.get(m, 0) + 1
        return z
    za, zc = cnt(a1), cnt(c1)
    print('  (c) accepted H1 by substitution class (share accepted / share library): ' + ', '.join(
        f'{m} {za.get(m, 0)} ({za.get(m, 0) / max(len(a1), 1) / max(zc.get(m, 0) / max(len(c1), 1), 1e-9):.2f}x)'
        for m in list(MODS) + ['other'] if za.get(m, 0) or zc.get(m, 0)), flush=True)
    # (d) set aside modification-equivalent entrapments (pre-listed classes), correct r
    am = int(ae['mod'].sum()); r2 = nE * (1 - frac_mod) / nH
    print(f'  (d) set aside {am} modification-equivalent accepted entrapments: lower {(E - am) / (N - am):.4f} '
          f'combined {(E - am) * (1 + 1 / r2) / (N - am):.4f}', flush=True)
