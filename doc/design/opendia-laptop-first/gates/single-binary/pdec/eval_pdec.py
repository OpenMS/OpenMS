"""PLAN_PDEC evaluation (entrapment labels read here only). Usage: eval_pdec.py <library.pqp> <runs dir> <D|H>.
Runs are <runs dir>/<D|H>_<cur|new>_s<seed>/<run>.precursors.tsv[.gz] (tool output as written, no re-competition).
Per run, at the tool's q <= 0.01: N, entrapments E, lower bound E/N, combined bound E(1+1/r)/N.
Per arm, pooled over seeds: lower and combined bound with a 90 % Poisson interval on E, median N@FDP_c 1 %.
Pair fairness: P(entrapment target beats its own predicted decoy) among pairs whose better member has q <= 0.05."""
import sqlite3, sys, glob, os, numpy as np, pandas as pd
from scipy.stats import chi2
L, R, P = sys.argv[1:4]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
tg = pq[pq.d == 0]; ent = tg.set_index('precursor_id').e
r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
key = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq = pq.assign(key=key)
pairs = pq[pq.d == 0].set_index('key')[['precursor_id', 'e']].join(
    pq[pq.d == 1].set_index('key')[['precursor_id']].rename(columns={'precursor_id': 'dec'}), how='inner')
print(f'library: target precursors {len(tg)}, r = {r:.4f}, targets with a paired decoy {len(pairs)} '
      f'({len(pairs) / len(tg):.4f})')
rows, pw = [], []
for d in sorted(glob.glob(f'{R}/{P}_*_s*')):
    if not os.path.isdir(d): continue
    arm, seed = os.path.basename(d).split('_')[1], int(d.split('_s')[-1])
    f = [x for x in glob.glob(f'{d}/*.precursors.tsv*') if '/precursors.tsv' not in x]
    t = pd.read_csv(f[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    tt = t[t.is_decoy == 0].copy(); tt['e'] = tt.precursor_id.map(ent).fillna(0).values
    acc = tt[tt.q_value <= 0.01]; N, E = len(acc), int(acc.e.sum())
    e = tt.sort_values(['q_value', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    rows.append(dict(arm=arm, seed=seed, N_q01=N, E_q01=E, lower=E / max(N, 1), combined=E * k / max(N, 1),
                     N_at_FDPc1=int(ok.max() + 1) if len(ok) else 0))
    s = t.set_index('precursor_id')
    x = pairs.assign(st=s.score.reindex(pairs.precursor_id).values, sd=s.score.reindex(pairs.dec).values,
                     qt=s.q_value.reindex(pairs.precursor_id).values, qd=s.q_value.reindex(pairs.dec).values)
    x = x.dropna(subset=['st', 'sd']); x = x[(np.minimum(x.qt, x.qd) <= 0.05) & (x.e == 1)]
    pw.append(dict(arm=arm, n=len(x), wins=int((x.st > x.sd).sum())))
df = pd.DataFrame(rows).sort_values(['arm', 'seed']); pw = pd.DataFrame(pw)
print(df.round(4).to_string(index=False))
for arm, g in df.groupby('arm'):
    E, N = g.E_q01.sum(), g.N_q01.sum()
    lo, hi = chi2.ppf(0.05, 2 * E) / 2 if E else 0.0, chi2.ppf(0.95, 2 * E + 2) / 2
    w = pw[pw.arm == arm]; n, wins = w.n.sum(), w.wins.sum()
    print(f'{arm}: pooled q<=0.01 lower {E / N:.4f} (90 % {lo / N:.4f}-{hi / N:.4f}), combined {E * k / N:.4f} '
          f'(90 % {lo * k / N:.4f}-{hi * k / N:.4f}); E {E}, N {N}; median N@FDP_c 1 % {g.N_at_FDPc1.median():.0f} '
          f'({g.N_at_FDPc1.min()}-{g.N_at_FDPc1.max()}); P(entrapment beats own decoy | best q<=0.05) '
          f'{wins / max(n, 1):.3f} (n {n}, 90 % +-{1.645 * np.sqrt(0.25 / max(n, 1)):.3f})')
