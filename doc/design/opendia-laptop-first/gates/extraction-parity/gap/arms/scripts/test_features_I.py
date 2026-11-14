"""Synthetic check of features_I.py (no real data): two precursors in one isolation window, overlapping peaks.
Precursor 1 (the claimant b) has 12 fragments; its fragment 0 has the same m/z as fragment 0 of precursor 2 (a).
b's unshared fragments show it at abundance 100 x lib. a has no own signal on fragment 0 beyond b's contribution, and
its other fragments match its library at abundance 10 x lib. Expected:
  - a's explained share = b's expected contribution on the shared fragment / a's total central area;
  - a's residual cosine to its library rises above the observed cosine (i_dlibcos > 0);
  - a third precursor c in another isolation window gets explained 0 and nclaim 0;
  - shifting b's apex by more than the 32-bin window removes the overlap (no link): explained 0.
Usage: test_features_I.py (writes a temporary root under $JCODE_SCRATCH_DIR or /tmp)."""
import json, os, sqlite3, subprocess, sys, tempfile, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
K = 12


def build(root, b_apex):
    os.makedirs(f'{root}/raw01/cal', exist_ok=True)
    json.dump({'ms2_mz': {'tolerance_half_width_ppm': 7.0}}, open(f'{root}/raw01/cal/cal.json', 'w'))
    pd.DataFrame({'lo': [400.0, 402.0], 'hi': [402.0, 404.0]}).to_csv(f'{root}/iso.tsv', sep='\t', index=False)
    pop = pd.DataFrame({'sample_id': [10, 20, 30], 'precursor_id': [1, 2, 3], 'rt_apex': [b_apex, 100.0, 100.0]})
    pop.to_csv(f'{root}/raw01/population.tsv', sep='\t', index=False)
    pd.DataFrame({'candidate_row': [10, 20, 30], 'rt_start': [b_apex - 3, 97.0, 97.0], 'rt_end': [b_apex + 3, 103.0, 103.0]}).to_csv(
        f'{root}/feat.tsv', sep='\t', index=False)
    db = f'{root}/lib12.pqp'
    if os.path.exists(db): os.remove(db)
    c = sqlite3.connect(db)
    c.execute('create table PRECURSOR(ID int, PRECURSOR_MZ real)')
    c.execute('create table TRANSITION(ID int, PRODUCT_MZ real)')
    c.execute('create table TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID int, PRECURSOR_ID int)')
    c.executemany('insert into PRECURSOR values (?,?)', [(1, 401.0), (2, 401.5), (3, 403.0)])
    tid = 0
    for p, base in ((1, 300.0), (2, 1300.0), (3, 2300.0)):
        for s in range(K):
            mz = 500.0 if s == 0 and p in (1, 2) else base + 17.3 * s
            c.execute('insert into TRANSITION values (?,?)', (tid, mz)); c.execute('insert into TRANSITION_PRECURSOR_MAPPING values (?,?)', (tid, p)); tid += 1
    c.commit(); c.close()
    lib = np.tile(np.linspace(1.0, 0.3, K), (3, 1)).astype(np.float32)
    prof = np.exp(-0.5 * ((np.arange(32) - 15.5) / 2.0) ** 2)
    x = np.zeros((3, K, 32), np.float32)
    x[0] = 100 * lib[0][:, None] * prof[None, :]      # claimant b: abundance 100
    x[1] = 10 * lib[1][:, None] * prof[None, :]       # a: abundance 10 ...
    x[1, 0] = 100 * lib[0][0] * prof                  # ... but its fragment 0 carries only b's signal
    x[2] = 10 * lib[2][:, None] * prof[None, :]
    np.save(f'{root}/raw01/xic12.npy', x); np.save(f'{root}/raw01/lib12.npy', lib)
    return x, lib, prof


def run(root):
    r = subprocess.run([PY, f'{HERE}/features_I.py', 'raw01', root, f'{root}/lib12.pqp', f'{root}/iso.tsv', f'{root}/feat.tsv'],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return np.load(f'{root}/raw01/featI.npy')


root = tempfile.mkdtemp(dir=os.environ.get('JCODE_SCRATCH_DIR'))
x, lib, prof = build(root, 100.0)
f = run(root)
area_a = x[1][:, 12:20].sum(1)
expect = min(100 * lib[0][0] * prof[12:20].sum(), area_a[0]) / area_a.sum()
print('a: explained', f[1, 0], 'expected', round(float(expect), 4), '| dlibcos', f[1, 2], '| nclaim', f[1, 3])
assert abs(f[1, 0] - expect) < 1e-4, 'explained share'
assert f[1, 2] > 0, 'residual cosine should rise'
assert f[2, 0] == 0 and f[2, 3] == 0, 'other isolation window must be untouched'
# b's shared fragment is explained by a at a's abundance (10 x lib), not b's own 100 x: small explained share for b.
assert 0 < f[0, 0] < f[1, 0], 'claimant gets a smaller explained share'
x, lib, prof = build(root, 200.0)  # no RT overlap
f = run(root)
assert f[1, 0] == 0 and f[1, 3] == 0, 'no overlap -> no link'
print('features_I synthetic checks: PASS')
