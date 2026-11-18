"""PLAN_PGENT evaluation (entrapment labels read here only). Usage: eval_pgent.py <library.pqp> <runs dir> <D|H>.
Per run: N at q <= 0.01, entrapments, FDP_c = N_e (1 + 1/r) / N; N@FDP_c 1 % from the precursor-winner ranking
(the tool's own q order); then per arm the pooled FDP_c over seeds and the median N@FDP_c 1 %."""
import sqlite3, sys, glob, os, numpy as np, pandas as pd
L, R, P = sys.argv[1:4]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
tg = pd.read_sql('select ID precursor_id from PRECURSOR where DECOY=0', c).precursor_id
r = (ent.reindex(tg) == 1).sum() / (ent.reindex(tg) == 0).sum(); k = 1 + 1 / r
print(f'library: target precursors {len(tg)}, r = {r:.4f}')
rows = []
for d in sorted(glob.glob(f'{R}/{P}_*_s*')):
    if not os.path.isdir(d): continue
    arm, seed = os.path.basename(d).split('_')[1], int(d.split('_s')[-1])
    f = glob.glob(f'{d}/*.precursors.tsv*'); f = [x for x in f if '/precursors.tsv' not in x]
    t = pd.read_csv(f[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    tt = t[t.is_decoy == 0].copy(); tt['e'] = tt.precursor_id.map(ent).fillna(0).values
    acc = tt[tt.q_value <= 0.01]; N, E = len(acc), int(acc.e.sum())
    # N@FDP_c 1 %: walk the tool's ranking of target winners, largest prefix with FDP_c <= 0.01
    e = tt.sort_values('q_value', kind='stable').e.values  # q order (ties keep score order)
    tt2 = tt.sort_values(['q_value', 'score'], ascending=[True, False]); e = tt2.e.values
    fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1)
    ok = np.where(fdp <= 0.01)[0]; n1 = int(ok.max() + 1) if len(ok) else 0
    rows.append(dict(arm=arm, seed=seed, N_q01=N, E_q01=E, FDPc_q01=E * k / max(N, 1), N_at_FDPc1=n1))
df = pd.DataFrame(rows).sort_values(['arm', 'seed'])
print(df.round(4).to_string(index=False))
for arm, g in df.groupby('arm'):
    print(f'{arm}: pooled FDP_c at q<=0.01 {g.E_q01.sum() * k / g.N_q01.sum():.4f} (E {g.E_q01.sum()}, N {g.N_q01.sum()}); '
          f'median N@FDP_c 1 % {g.N_at_FDPc1.median():.0f} ({g.N_at_FDPc1.min()}-{g.N_at_FDPc1.max()}); median N at q<=0.01 {g.N_q01.median():.0f}')
