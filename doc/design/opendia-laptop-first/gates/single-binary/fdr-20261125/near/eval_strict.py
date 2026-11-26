#!/usr/bin/env python3
"""STRICT evaluation (diagnosis): the step A relabel evaluator, plus a strict homolog-free entrapment definition.
Entrapment peptides in the STRICT list (min Hamming <= 2 to a same-length human window, I=L; strict_list.py, library only)
are AMBIGUOUS: they are dropped from accepted N and E, and from the library entrapment count used for r. Everything else
is as eval_pdec_relabel.py (tool q as written, no re-competition). Arms: layout <runs dir>/<SET>_<arm>_s<seed>/.
Prints per arm: pooled lower / combined (90 % Poisson) under (a) relabel only and (b) strict; median N@FDP_c 1 % (b).
Usage: eval_strict.py <lib.pqp> <runs dir> <SET>   env RELABEL=<audit tsv> STRICT=<strict list tsv>"""
import sqlite3, sys, glob, os, numpy as np, pandas as pd
from scipy.stats import chi2
L, R, S = sys.argv[1:4]
RELABEL = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq) if os.environ.get('RELABEL') else set()
STRICT = set(pd.read_csv(os.environ['STRICT'], sep='\t').seq)
# env HUMNEAR=<tsv>: human library peptides with a near human neighbour (hd <= 2, strict_list.py HUMAN=1). Adds a
# STRATIFIED estimate: near-homolog entrapments count as errors, each stratum s (near / far) uses its own
# r_s = entrapment_s / human_s, FDP = sum_s E_s (1 + 1 / r_s) / N over all accepted (nothing removed).
HUMNEAR = set(pd.read_csv(os.environ['HUMNEAR'], sep='\t').seq) if os.environ.get('HUMNEAR') else None
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
pq.loc[(pq.d == 0) & pq.seq.isin(RELABEL), 'e'] = 0
tg = pq[pq.d == 0].copy(); tg['amb'] = (tg.e == 1) & tg.seq.isin(STRICT)
info = tg.set_index('precursor_id')[['e', 'amb']]
if HUMNEAR is not None:
    tg['near'] = np.where(tg.e == 1, tg.amb, tg.seq.isin(HUMNEAR))
    rs = {s: ((tg.e == 1) & (tg.near == s)).sum() / ((tg.e == 0) & (tg.near == s)).sum() for s in (True, False)}
    print(f'stratified: r near {rs[True]:.4f} (human near {int(((tg.e == 0) & tg.near).sum())}), r far {rs[False]:.4f}')
    pnear = tg.set_index('precursor_id').near
r1 = (tg.e == 1).sum() / (tg.e == 0).sum(); r2 = ((tg.e == 1) & ~tg.amb).sum() / (tg.e == 0).sum()
k1, k2 = 1 + 1 / r1, 1 + 1 / r2
print(f'{S}: library target precursors {len(tg)}; entrapment {int((tg.e == 1).sum())} (ambiguous {int(tg.amb.sum())}); '
      f'r relabel {r1:.4f}, r strict {r2:.4f}')
rows = []
for d in sorted(glob.glob(f'{R}/{S}_*_s*')):
    if not os.path.isdir(d): continue
    arm, seed = os.path.basename(d).split('_')[1], int(d.split('_s')[-1])
    f = [x for x in glob.glob(f'{d}/*.precursors.tsv*') if '/precursors.tsv' not in x]
    t = pd.read_csv(f[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
    t = t[(t.score > -1e29) & (t.is_decoy == 0)].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    t = t.join(info, on='precursor_id'); t['e'] = t.e.fillna(0); t['amb'] = t.amb.fillna(False).astype(bool)
    acc = t[t.q_value <= 0.01]
    N1, E1 = len(acc), int(acc.e.sum()); a2 = acc[~acc.amb]; N2, E2 = len(a2), int(a2.e.sum())
    u = t[~t.amb].sort_values(['q_value', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(u) * k2 / np.arange(1, len(u) + 1); ok = np.where(fdp <= 0.01)[0]
    rows.append(dict(arm=arm, seed=seed, N1=N1, E1=E1, N2=N2, E2=E2, amb=int(acc.amb.sum()), N_at=int(ok.max() + 1) if len(ok) else 0,
                     En=int((acc.e * acc.amb).sum()), Ef=int((acc.e * ~acc.amb).sum()),
                     Hn=int(((acc.e == 0) & acc.precursor_id.map(pnear).fillna(False).astype(bool)).sum()) if HUMNEAR is not None else 0))
df = pd.DataFrame(rows).sort_values(['arm', 'seed'])
print(df.to_string(index=False))
ci = lambda E: (chi2.ppf(0.05, 2 * E) / 2 if E else 0.0, chi2.ppf(0.95, 2 * E + 2) / 2)
for arm, g in df.groupby('arm'):
    E1, N1, E2, N2 = g.E1.sum(), g.N1.sum(), g.E2.sum(), g.N2.sum(); l1, h1 = ci(E1); l2, h2 = ci(E2)
    print(f'{arm}: relabel lower {E1 / N1:.4f} combined {E1 * k1 / N1:.4f} ({l1 * k1 / N1:.4f}-{h1 * k1 / N1:.4f}) | '
          f'strict lower {E2 / N2:.4f} combined {E2 * k2 / N2:.4f} ({l2 * k2 / N2:.4f}-{h2 * k2 / N2:.4f}); '
          f'ambiguous accepted {g.amb.sum()}; strict median N@FDP_c 1 % {g.N_at.median():.0f}')
    if HUMNEAR is not None:
        st = (g.En.sum() * (1 + 1 / rs[True]) + g.Ef.sum() * (1 + 1 / rs[False])) / N1
        Nn = g.En.sum() + g.Hn.sum(); Nf = N1 - Nn
        print(f'{arm}: stratified combined {st:.4f} (accepted E near {g.En.sum()}, far {g.Ef.sum()}); '
              f'near stratum N {Nn} FDP_c {g.En.sum() * (1 + 1 / rs[True]) / Nn:.4f} | far stratum N {Nf} FDP_c {g.Ef.sum() * (1 + 1 / rs[False]) / Nf:.4f}')
