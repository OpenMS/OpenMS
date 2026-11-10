"""PLAN7.md reading: P1 (refined peak picking for the boundary group) on raws 01/02, ID rule + bootstrap, quant rule
(naive vs refined on fixed base windows), references and mechanism (reported only). q-derived counts are decoy-count estimates."""
import os, re, sqlite3, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
sys.path.insert(0, B + '/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled
O = B + '/astral-single-scorer-20261109'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
RULE = {'raw01': (1050, 0.015, 0.070), 'raw02': (1324, 0.015, 0.079)}
BASE = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
rng = np.random.default_rng(7)
def winners(d):
    cr = pd.read_csv(d + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'precursor_best', 'in_rt_window', 'stage2_score'])
    return cr, cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id')
def n_at(score, e, t):
    o = np.argsort(-score, kind='stable'); e = e[o]
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(fdp <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0
def accepted(w, t):
    w = w.sort_values('stage2_score', ascending=False)
    fdp = np.cumsum(w.e.values) * (1 + 1 / R) / np.arange(1, len(w) + 1); ok = np.where(fdp <= t)[0]
    return w.iloc[:ok.max() + 1 if len(ok) else 0]
def usage(path):
    s = open(path).read() if os.path.exists(path) else ''
    wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', s); rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', s)
    parts = [float(x) for x in wall.group(1).split(':')] if wall else [float('nan')]
    return sum(v * 60 ** (len(parts) - 1 - i) for i, v in enumerate(parts)) / 60, (int(rss.group(1)) / 1024 ** 2 if rss else float('nan'))
def auc(e, d): return mannwhitneyu(e, d).statistic / (len(e) * len(d))
cmp_txt = open(O + '/CMP7.txt').read() if os.path.exists(O + '/CMP7.txt') else ''
id_ok, base_acc = [], {}
print('== IDs')
for raw in ('raw01', 'raw02'):
    d = f'{O}/{raw}/P1'
    crb, wb = winners(BASE[raw]); base_acc[raw] = accepted(wb, 0.01)
    if not os.path.exists(d + '/crossfit.tsv'):
        print(f'{raw} P1 missing'); id_ok.append(False); continue
    cr, w = winners(d)
    curve = [n_at(w.stage2_score.values, w.e.values, t) for t in (0.005, 0.01, 0.015, 0.02)]
    best, _ = paired_pooled(cr.merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left'), 'stage2_score')
    t = best[best.is_decoy == 0].merge(ent, on='precursor_id')
    fdp = {q: t[t.q <= q].e.sum() * (1 + 1 / R) / max((t.q <= q).sum(), 1) for q in (0.01, 0.05)}
    nq = {q: int((t.q <= q).sum()) for q in (0.01, 0.05)}
    minutes, gb = usage(d + '/stderr.log')
    internal = f'{raw}/P1: crossfit cols 1-6,9 IDENTICAL' in cmp_txt and f'{raw}/P1: features cols 1-9,15-16 IDENTICAL' in cmp_txt
    n1, f1, f5 = RULE[raw]
    ok = curve[1] >= n1 and fdp[0.01] <= f1 and fdp[0.05] <= f5 and gb <= 30 and minutes <= 35 and internal
    id_ok.append(ok)
    print(f'{raw} P1 N at FDP_c 0.5/1/1.5/2%: {curve} | q<=0.01: {nq[0.01]} FDP_c {fdp[0.01]:.3f} | q<=0.05: {nq[0.05]} FDP_c {fdp[0.05]:.3f}'
          f' | {minutes:.1f} min {gb:.1f} GB | internal {internal} | rule {"PASS" if ok else "fail"}', flush=True)
    # Paired bootstrap against base.
    tab = wb[['precursor_id', 'stage2_score', 'e']].rename(columns={'stage2_score': 'b'}).merge(
        w[['precursor_id', 'stage2_score']].rename(columns={'stage2_score': 'a'}), on='precursor_id', how='outer')
    tab['e'] = tab.precursor_id.map(ent.set_index('precursor_id').e); tab = tab.fillna({'a': -np.inf, 'b': -np.inf})
    e = tab.e.values.astype(int); diffs = []
    for _ in range(1000):
        idx = rng.integers(0, len(tab), len(tab))
        diffs.append([n_at(tab.a.values[idx], e[idx], q) - n_at(tab.b.values[idx], e[idx], q) for q in (0.01, 0.02)])
    diffs = np.array(diffs)
    lo1, hi1 = np.percentile(diffs[:, 0], [5, 95]); lo2, hi2 = np.percentile(diffs[:, 1], [5, 95])
    print(f'   bootstrap P1 - base: @1% 90% {lo1:.0f}..{hi1:.0f} | @2% 90% {lo2:.0f}..{hi2:.0f}{" | CLEAR LOSS" if hi2 < 0 else ""}')
    # Mechanism.
    f = pd.read_csv(d + '/features.tsv', sep='\t').merge(cr[['sample_id', 'in_rt_window']], on='sample_id')
    fb = pd.read_csv(BASE[raw] + '/features.tsv', sep='\t').merge(crb[['sample_id', 'in_rt_window']], on='sample_id')
    for name, g in (('base', fb), ('P1', f)):
        g = g[g.in_rt_window == 1].sort_values(['stage2_score', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').merge(ent, on='precursor_id', how='left')
        ee = g[(g.is_decoy == 0) & (g.e == 1)]; dd = g[g.is_decoy == 1]
        print(f'   AUC(entrapment > decoy) {name}: ' + ' '.join(f'{x}={auc(ee[x].values, dd[x].values):.3f}' for x in ('boundary_log_width', 'boundary_coelution', 'stage2_score')))
print(f'ID rule: {"PASS" if all(id_ok) else "FAIL"} ({id_ok})')

print('== Quant')
qt = {}
for raw in ('raw01', 'raw02'):
    path = f'{O}/{raw}/P1/quant.tsv'
    if not os.path.exists(path): print(f'{raw} quant missing'); sys.exit(0)
    q = pd.read_csv(path, sep='\t')
    a = base_acc[raw]; a = a[a.e == 0][['precursor_id', 'sample_id']]
    qt[raw] = a.merge(q.drop(columns='precursor_id'), on='sample_id')
    changed = (qt[raw].refined_left != qt[raw].naive_left) | (qt[raw].refined_right != qt[raw].naive_right)
    wn = (qt[raw].naive_right - qt[raw].naive_left).median(); wr = (qt[raw].refined_right - qt[raw].refined_left).median()
    print(f'{raw}: base 1% non-entrapment IDs {len(a)}, quantified rows {len(qt[raw])} | boundaries changed {changed.mean():.3f} | median width naive {wn:.2f}s refined {wr:.2f}s'
          f' | quantitative fragments hist ' + ' '.join(f'{k}:{(qt[raw].quantitative == k).mean():.2f}' for k in range(7)))
m = qt['raw01'].merge(qt['raw02'], on='precursor_id', suffixes=('_1', '_2'))
n_set = len(m)
both = np.ones(len(m), bool); num = np.zeros(len(m)); den = np.zeros(len(m)); kq = np.zeros(len(m), int)
for k in range(6):
    sel = (m[f'r{k}_1'] >= 0.9) & (m[f'r{k}_2'] >= 0.9) & (m[f'area{k}_1'] > 0) & (m[f'area{k}_2'] > 0)
    num += np.where(sel, m[f'area{k}_1'], 0); den += np.where(sel, m[f'area{k}_2'], 0); kq += sel.values
m['ref1'], m['ref2'], m['kq'] = num, den, kq
use = (m.kq >= 2) & (m.naive_area_1 > 0) & (m.naive_area_2 > 0)
u = m[use]
def rsd(x):
    x = x - np.median(x); return 1.4826 * np.median(np.abs(x))
ln = np.log2(u.naive_area_1 / u.naive_area_2).values; lr = np.log2(u.ref1 / u.ref2).values
sd_n, sd_r = rsd(ln), rsd(lr)
boot = []
for _ in range(1000):
    idx = rng.integers(0, len(u), len(u)); boot.append(rsd(lr[idx]) - rsd(ln[idx]))
lo, hi = np.percentile(boot, [5, 95])
cov = use.mean()
ok_q = sd_r <= 0.95 * sd_n and hi < 0 and cov >= 0.8
print(f'precursor set (base 1% non-entrapment in both raws) {n_set}; refined quantifies {use.sum()} ({cov:.3f})')
print(f'robust SD log2(01/02): naive {sd_n:.4f} refined {sd_r:.4f} (ratio {sd_r / sd_n:.3f}); bootstrap 90% of refined-naive {lo:.4f}..{hi:.4f}')
print(f'Quant rule: {"PASS" if ok_q else "FAIL"}')
# References (reported only).
d = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                    columns=['Run', 'Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Precursor.Quantity', 'Decoy'])
d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0) & (d.Run.str.endswith('_01'))]
pep = pd.read_sql("""select p.ID precursor_id, pe.MODIFIED_SEQUENCE seq, p.CHARGE z from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID""", c)
d['seq'] = d['Modified.Sequence'].str.replace(r'\(UniMod:(\d+)\)', r'(UniMod:\1)', regex=True)
x = u.merge(pep, on='precursor_id').merge(d, left_on=['seq', 'z'], right_on=['seq', 'Precursor.Charge'])
if len(x) > 50:
    print(f'reference: Spearman with DIA-NN Precursor.Quantity (raw 01, n={len(x)}): naive {spearmanr(x.naive_area_1, x["Precursor.Quantity"])[0]:.3f} '
          f'refined {spearmanr(x.ref1, x["Precursor.Quantity"])[0]:.3f}; DIA-NN 01/02 robust SD 0.335 (72k)')
else:
    print(f'reference: DIA-NN match too small (n={len(x)}); sequence format differs')
