"""PRECHECK_GAP2 part B evaluation (diagnostic, 22 Nov 2026).
Per run (D, H) and level (r0, r7): arms def / V7M6 / V8M8 on the sibling default (opendia-pd10-sibdef.sif, seed 0).
HeLa IDs at own q<=0.01 (best target row per precursor in <run>.precursors.tsv), at r7 also N@FDP_c 1 % and bounds,
wall time and peak RSS from the /usr/bin/time log.
Usage: precheck_gap2_b.py <gap2 dir>   (expects <dir>/<R>_<L>_<arm>/ and <dir>/<R>_<L>_<arm>.log)"""
import sys, glob, re, sqlite3, numpy as np, pandas as pd
W = sys.argv[1]
C = '/nfs/research/juan/singj/opendia'
LIB = {'r0': f'{C}/libsize-20261121/r0/lib.pqp', 'r7': f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'}
R7 = 7.0539; k = 1 + 1 / R7
def wall(log):
    t = open(log).read(); m = re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([0-9:.]+)', t)
    rss = re.findall(r'Maximum resident set size \(kbytes\): (\d+)', t)
    p = [float(x) for x in m[-1].split(':')]; s = p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) > 2 else 0)
    return s, int(rss[-1]) / 1e6
ents = {}
for L, lib in LIB.items():
    c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
    ents[L] = pd.read_sql("""select p.ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
      join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
      where p.DECOY=0 group by p.ID""", c).set_index('pid').e
rows = []
for R in ('D', 'H'):
    for L in ('r0', 'r7'):
        for arm in ('def', 'V7M6', 'V8M8'):
            rd, log = f'{W}/{R}_{L}_{arm}', f'{W}/{R}_{L}_{arm}.log'
            try: txt = open(log).read()
            except OSError: txt = ''
            if 'exit 0' not in txt: print(f'{R} {L} {arm}: missing or failed'); continue
            pt = [f for f in glob.glob(f'{rd}/*.precursors.tsv') if not f.endswith('/precursors.tsv')][0]
            t = pd.read_csv(pt, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
            t = t[(t.is_decoy == 0) & (t.score > -1e29)].sort_values('score', ascending=False).drop_duplicates('precursor_id')
            t['e'] = t.precursor_id.map(ents[L]).fillna(0).values
            acc = t[t.q_value <= 0.01]; N, E = len(acc), int(acc.e.sum())
            row = dict(run=R, lvl=L, arm=arm, N=N, hela=N - E)
            if L == 'r7':
                e = t.sort_values(['q_value', 'score'], ascending=[True, False]).e.values
                fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
                row.update(lower=E / N, combined=E * k / N, n_fdp1=int(ok.max() + 1))
            row['wall_s'], row['rss_gb'] = wall(log)
            rows.append(row)
df = pd.DataFrame(rows)
for (R, L), g in df.groupby(['run', 'lvl']):
    b = g[g.arm == 'def']
    if b.empty: continue
    b = b.iloc[0]
    df.loc[g.index, 'd_hela'] = g.hela / b.hela - 1; df.loc[g.index, 'wall_x'] = g.wall_s / b.wall_s
    if L == 'r7': df.loc[g.index, 'd_fdp1'] = g.n_fdp1 / b.n_fdp1 - 1; df.loc[g.index, 'd_lower_pp'] = (g.lower - b.lower) * 100
pd.set_option('display.width', 250)
print(df.round(4).to_string(index=False))
for arm in ('V7M6', 'V8M8'):
    a = df[df.arm == arm]
    if len(a) < 4: print(f'RULE B {arm}: incomplete'); continue
    r0 = a[a.lvl == 'r0']; r7 = a[a.lvl == 'r7']
    ok = [bool((r0.d_hela >= 0.03).all()), bool((r7.d_fdp1 >= 0.02).all()), bool((r7.d_lower_pp <= 0.2).all()), bool((a.wall_x <= 1.5).all())]
    print(f'RULE B {arm}: r0 HeLa>=+3% {ok[0]}, r7 FDP_c>=+2% {ok[1]}, lower<=+0.2pp {ok[2]}, wall<=1.5x {ok[3]} -> '
          f'{"PASS" if all(ok) else "FAIL"}; mean r7 gain {r7.d_fdp1.mean():+.4f}')
