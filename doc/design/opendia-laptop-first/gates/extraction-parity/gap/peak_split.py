"""PRECHECK_PEAK (diagnostic; DIA-NN 2.2 as a reference annotation only): break the 'wrong peak' losses of the current
C++ scorer (arm-E arrays + C++ featI, seeds 0-4) into causes, and test label-free alternatives to the W rule offline.
Per DIA-NN human q<=0.01 ID in our library (the 'reference'):
  A no candidate at all in the search table;  B candidates, none within 6 s of DIA-NN's apex in the search table
  (B1 DIA-NN apex outside our RT window, B2 inside: picker/cap miss);  C a candidate within 6 s exists but the pair
  population drops it (support rule);  D in the population, a within-6 s window exists, but another window is kept;
  E kept = DIA-NN's peak.
What-ifs (same scores, only the per-precursor window choice changes): keep the window with the highest full score
instead of the largest linear part."""
import json, sqlite3, sys, numpy as np, pandas as pd
RAW = sys.argv[1]
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'; B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = f'{G}/armE_{RAW}'
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,p.LIBRARY_RT lrt,pe.MODIFIED_SEQUENCE ms,p.DECOY decoy from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2)
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key')
REF = {'raw01': B + '/diann-pairlib-20261113/report_raw01.parquet', 'raw02': B + '/diann-decoys-20261114/own/raw02/report.parquet'}[RAW]
ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid'); ref['rt_s'] = ref.RT * 60
H = set(ref.index)
sf = pd.read_csv(f'{D}/search.features.tsv', sep='\t', usecols=['candidate_row', 'precursor_id', 'rt_apex', 'transition_support', 'peak_rank', 'peaks_in_window'])
pop = pd.read_csv(f'{D}/population.tsv', sep='\t'); pop['row'] = np.arange(len(pop))
cal = json.load(open(f'{D}/cal/cal.json'))['rt']; T = cal['tolerance_half_width_s']
kx, ky = np.array(cal['knots_library_rt']), np.array(cal['knots_run_rt_s'])
lib_rt = pre.set_index('precursor_id').lrt
def m(x):
    x = np.asarray(x, float); p = np.interp(x, kx, ky); lo = x < kx[0]; hi = x > kx[-1]
    p[lo] = ky[0] + (x[lo] - kx[0]) * cal['low_slope_s']; p[hi] = ky[-1] + (x[hi] - kx[-1]) * cal['high_slope_s']; return p
sfh = sf[sf.precursor_id.isin(H)].copy(); sfh['d'] = (sfh.rt_apex.values - ref.rt_s.reindex(sfh.precursor_id).values)
in_search = set(sfh.precursor_id); near_search = set(sfh[sfh.d.abs() <= 6].precursor_id)
poph = pop[pop.precursor_id.isin(H)].copy(); poph['d'] = poph.rt_apex.values - ref.rt_s.reindex(poph.precursor_id).values
in_pop = set(poph.precursor_id); near_pop = set(poph[poph.d.abs() <= 6].precursor_id)
A = H - in_search
Bset = in_search - near_search
win_err = np.abs(ref.loc[list(Bset), 'rt_s'].values - m(lib_rt.loc[list(Bset)].values)) if Bset else np.array([])
B1 = {p for p, e in zip(Bset, win_err) if e > T}; B2 = Bset - B1
C = near_search - near_pop
print(f'{RAW}: reference human IDs {len(H)}; window half width T = {T:.1f} s')
print(f'  A no search candidate {len(A)} | B no candidate within 6 s {len(Bset)} (B1 DIA-NN apex outside window {len(B1)}, B2 inside {len(B2)}) | C dropped by population {len(C)}')
b2 = sfh[sfh.precursor_id.isin(B2)].groupby('precursor_id').agg(n=('rt_apex', 'size'), cap=('peaks_in_window', 'max'), dmin=('d', lambda v: np.abs(v).min()))
print(f'  B2 detail: candidates per precursor {b2.n.value_counts().sort_index().to_dict()} | peaks found in window (cap 5) median {b2.cap.median():.0f}, share at cap {np.mean(b2.cap >= 5):.2f} | nearest candidate |d| p50 {b2.dmin.median():.1f} s')
# D vs E per seed and the what-if.
rows = []
for s in range(5):
    sc = np.load(f'{G}/cppscore/{RAW}_s{s}.npy').astype(float)
    kept = sc > -1e29
    w, best = tdc(pop, sc)
    tw = w[w.is_decoy == 0].sort_values('s', ascending=False); e = tw.precursor_id.map(ent).values
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(tw) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1); acc = set(tw.precursor_id.iloc[:n1])
    W = w.set_index('precursor_id')
    keptd = {p: abs(pop.rt_apex.values[int(W.row[p])] - ref.rt_s[p]) for p in near_pop if p in W.index}
    E = {p for p, d in keptd.items() if d <= 6}; Dw = near_pop - E
    rows.append(dict(seed=s, D_other_window_kept=len(Dw), E_right=len(E), right_accepted=len(E & acc), N1=n1))
x = pd.DataFrame(rows); print(x.to_string(index=False)); print('  median', x.median().to_dict())
