"""Suite evaluation, old (develop merge-base) vs new (total peptide order).
Per file: PSMs at q <= 0.01 (mean of seeds 1/42/137, normal database), delta, z = delta / sqrt(var_new/3 + var_old/3)
(up/down at |z| >= 2), PIN identical, PIN rows that differ; entrapment PSMs and FDP (r = 1 shuffled: 2 N_E / N,
all proteins ENTRAP_; Velos: frozen natural entrapment, N_E (1 + 1/r) / N). Groups: N and N_E summed over files.
usage: evaluate.py <out dir>"""
import collections, json, math, statistics, subprocess, sys
from pathlib import Path
S = Path('/home/user/bench/suite'); R = S/'runs'; OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
HY = '/home/user/bench/frozen/z/prose-andes-reproduction/broad_benchmark/data/hy_td.fasta'
DS = json.loads((S/'cfg/datasets.json').read_text())
group = lambda d: next(g for g in ('velos', 'hfx', 'astral', 'lumos_lfq', 'lumos_tmt', 'eclipse_tmtpro', 'tims') if d.startswith(g))
def diff_rows(a, b):
    A = open(a).read().splitlines(); B = open(b).read().splitlines()
    return sum(x != y for x, y in zip(A, B)) + abs(len(A) - len(B))
def velos(arm, d):
    files = [str(R/d/f'{arm}_normal/s{s}/target.tsv') for s in (1, 42, 137)]
    out = subprocess.run(['python3', '-I', str(S/'scripts/frozen_eval.py'), 'velos', HY] + files, capture_output=True, text=True, check=True).stdout
    r = float(out.split('r = ')[1].split(';')[0])
    rows = [l.split('\t') for l in out.splitlines() if l.startswith('/')]
    return [int(x[1]) for x in rows], [int(x[2]) for x in rows], 1 + 1 / r
rows = []; G = collections.defaultdict(lambda: collections.defaultdict(float))
for d in DS:
    m = {k: json.loads((R/d/f'{k}/metrics.json').read_text()) for k in ('old_normal', 'new_normal', 'old_entrap', 'new_entrap') if (R/d/k/'metrics.json').exists()}
    if 'old_normal' not in m or 'new_normal' not in m: continue
    ps = {a: [v['psms'] for v in m[f'{a}_normal']['seeds'].values()] for a in ('old', 'new')}
    mean = {a: statistics.mean(ps[a]) for a in ps}
    delta = mean['new'] - mean['old']; se = math.sqrt(statistics.variance(ps['new']) / 3 + statistics.variance(ps['old']) / 3)
    z = delta / se if se else (0.0 if delta == 0 else math.copysign(math.inf, delta))
    call = 'up' if z >= 2 else 'down' if z <= -2 else 'same'
    same_pin = m['old_normal']['pin_sha256'] == m['new_normal']['pin_sha256']
    nd = 0 if same_pin else diff_rows(R/d/'old_normal/input.pin', R/d/'new_normal/input.pin')
    row = dict(file=d, group=group(d), old=round(mean['old'], 1), new=round(mean['new'], 1), delta=round(delta, 1),
               delta_pct=round(100 * delta / mean['old'], 3), z=round(z, 2), call=call, seeds_old=ps['old'], seeds_new=ps['new'],
               pin_identical=same_pin, pin_rows=m['old_normal']['pin']['rows_target'] + m['old_normal']['pin']['rows_decoy'], pin_rows_differ=nd)
    g = G[group(d)]; g['files'] += 1; g['old'] += mean['old']; g['new'] += mean['new']; g['up'] += call == 'up'; g['down'] += call == 'down'
    g['pin_identical'] += same_pin
    for a in ('old', 'new'):
        if group(d) == 'velos':
            n, ne, f = velos(a, d)
        elif f'{a}_entrap' in m:
            sd = m[f'{a}_entrap']['seeds'].values(); n = [v['psms'] for v in sd]; ne = [v['entrapment_psms'] for v in sd]; f = 2.0
        else:
            continue
        row[f'fdp_{a}_pct'] = round(100 * f * sum(ne) / sum(n), 2); row[f'N_{a}'] = round(statistics.mean(n), 1)
        g[f'N_{a}'] += sum(n); g[f'NE_{a}'] += sum(ne); g['f'] = f
    rows.append(row)
cols = list(dict.fromkeys(k for r in rows for k in r))
with open(OUT/'per_file.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
tot = collections.defaultdict(float)
with open(OUT/'per_group.tsv', 'w') as f:
    f.write('group\tfiles\told\tnew\tdelta_pct\tup\tdown\tpin_identical\tfdp_old_pct\tfdp_new_pct\n')
    for k, g in G.items():
        fd = lambda a: f"{100 * g['f'] * g[f'NE_{a}'] / g[f'N_{a}']:.2f}" if g[f'N_{a}'] else ''
        f.write(f"{k}\t{int(g['files'])}\t{g['old']:.1f}\t{g['new']:.1f}\t{100 * (g['new'] - g['old']) / g['old']:+.3f}\t{int(g['up'])}\t{int(g['down'])}\t{int(g['pin_identical'])}\t{fd('old')}\t{fd('new')}\n")
        for x in ('files', 'old', 'new', 'up', 'down', 'pin_identical'): tot[x] += g[x]
    f.write(f"total\t{int(tot['files'])}\t{tot['old']:.1f}\t{tot['new']:.1f}\t{100 * (tot['new'] - tot['old']) / tot['old']:+.3f}\t{int(tot['up'])}\t{int(tot['down'])}\t{int(tot['pin_identical'])}\t\t\n")
print(open(OUT/'per_group.tsv').read())
