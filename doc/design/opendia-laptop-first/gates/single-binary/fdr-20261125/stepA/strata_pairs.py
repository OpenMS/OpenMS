#!/usr/bin/env python3
"""Step A, HeLa: where does the entrapment-vs-own-decoy asymmetry live? Diagnosis only.
For each (entrapment target, its own DECOY_ partner) pair where the better member has q <= 0.05 in a run table,
P = fraction where the entrapment scores higher (0.5 under exchangeability). Pooled over the given run tables, split by
library strata of the entrapment: charge, length, precursor m/z, library RT, IM, #transitions, C-term residue,
missed cleavages, Met, Pro, and the entrapment/decoy near-isobar exposure (same charge, 10 ppm, >= 2 shared fragments).
Also: own q <= 0.01 entrapments per stratum vs that stratum's entrapment share -> class bound E_c (1 + 1/r_c) / N.
Usage: strata_pairs.py <lib.pqp> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT rt, p.LIBRARY_DRIFT_TIME im,
  pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
nt = pd.read_sql("select PRECURSOR_ID id, count(*) nt from TRANSITION_PRECURSOR_MAPPING group by PRECURSOR_ID", c)
P = P.merge(nt, on='id', how='left')
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
T = P[P.d == 0]; D = P[P.d == 1]
pairs = T.set_index('key').join(D.set_index('key')[['id']].rename(columns={'id': 'did'}), how='inner')
pairs = pairs[pairs.e == 1].reset_index()
r = (T.e == 1).sum() / (T.e == 0).sum()
s = pairs
s['len'] = s.seq.str.len()
strata = {
    'charge': s.z.clip(upper=4).astype(str),
    'length': pd.cut(s['len'], [0, 9, 12, 16, 22, 60]).astype(str),
    'mz': pd.qcut(s.mz, 4).astype(str),
    'lib_rt': pd.qcut(s.rt, 4).astype(str),
    'im': pd.qcut(s.im, 4).astype(str),
    'n_trans': s.nt.clip(upper=8).astype(str),
    'cterm': s.seq.str[-1].where(s.seq.str[-1].isin(['K', 'R']), 'other'),
    'missed': s.seq.str[:-1].str.count('[KR]').clip(upper=2).astype(str),
    'Met': s.seq.str.contains('M').map({True: 'M', False: 'noM'}),
    'Pro': s.seq.str.contains('P').map({True: 'P', False: 'noP'}),
}
share = {k: v.value_counts(normalize=True) for k, v in strata.items()}
eid = pd.Series(np.arange(len(s)), index=s.id)
print(f'library: entrapment pairs {len(s)}, r {r:.4f}', flush=True)
for spec in sys.argv[2:]:
    arm, g = spec.split('=', 1)
    files = sorted(glob.glob(g))
    W = np.zeros(len(s)); Nn = np.zeros(len(s)); acc = np.zeros(len(s)); N = 0
    for f in files:
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        st, sd = t.score.reindex(s.id).values, t.score.reindex(s.did).values
        qt, qd = t.q_value.reindex(s.id).values, t.q_value.reindex(s.did).values
        ok = ~np.isnan(st) & ~np.isnan(sd) & (np.fmin(qt, qd) <= 0.05)
        Nn += ok; W += ok & (st > sd)
        a = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]
        N += len(a)
        acc += np.isin(s.id.values, a.index.values)
    E = acc.sum()
    print(f'== {arm}: runs {len(files)} N {N} E {int(E)} lower {E / N:.4f} combined {E * (1 + 1 / r) / N:.4f}; '
          f'P {W.sum() / Nn.sum():.3f} (n {int(Nn.sum())})', flush=True)
    for k, v in strata.items():
        rows = []
        for lev in sorted(v.unique()):
            m = (v == lev).values
            n, w, e = Nn[m].sum(), W[m].sum(), acc[m].sum(); rc = r * share[k][lev]
            rows.append(f'{lev}: share {share[k][lev]:.3f} P {w / max(n, 1):.3f} (n {int(n)}) E {int(e)} bound {e * (1 + 1 / rc) / N * 1 / 1:.4f}')
        print(f'  {k:8s} ' + ' | '.join(rows), flush=True)
