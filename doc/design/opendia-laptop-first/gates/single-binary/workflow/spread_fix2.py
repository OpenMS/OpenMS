"""Offline check, variant 2: fit in window coordinates. Every stage-3 winner's residual to its window centre lies in
[-hw, hw], targets and decoys alike, so the decoy background has the targets' window edges exactly. A true target peak
sits at shift_m + N(mu, s), shift_m = final map - window centre at the target's library RT. EM over pi, mu, s.
Usage: spread_fix2.py <library.pqp> <cal dir> ..."""
import json, sys, sqlite3, numpy as np, pandas as pd
LIB, DIRS = sys.argv[1], sys.argv[2:]
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pre = pd.read_sql("select ID precursor_id, TRAML_ID traml, DECOY decoy, LIBRARY_RT lrt from PRECURSOR", c)
tl = pre[pre.decoy == 0].set_index('traml').lrt
pre['target_lrt'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tl), pre.lrt)
lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])
for d in DIRS:
    cal = json.load(open(f'{d}/cal.json')); kx, ky = np.array(cal['rt']['knots_library_rt']), np.array(cal['rt']['knots_run_rt_s'])
    hw = cal['stage3']['half_window_s']
    t = pd.read_csv(f'{d}/diag.tsv', sep='\t', comment='#'); t = t[t.stage3_score.notna()].merge(pre[['precursor_id', 'target_lrt']], on='precursor_id')
    res = (t.stage3_rt - t.stage3_centre_rt).values; tg = (t.decoy == 0).values
    shift = (np.interp(t.target_lrt, kx, ky) - t.stage3_centre_rt).values
    keep = np.abs(res) <= hw
    x, sh = res[tg & keep], shift[tg & keep]
    e = np.linspace(-hw, hw, 41); h, _ = np.histogram(res[~tg & keep], e); h = (h + 1) / (h.sum() + len(h)) / np.diff(e)
    bg = h[np.clip(np.searchsorted(e, x) - 1, 0, len(h) - 1)]
    pi, mu, s = 0.3, 0.0, hw / 4
    for _ in range(300):
        z = x - sh
        g = np.exp(-0.5 * ((z - mu) / s) ** 2) / (s * np.sqrt(2 * np.pi)); q = pi * g / (pi * g + (1 - pi) * bg)
        pi, mu = q.mean(), (q * z).sum() / q.sum(); s = np.sqrt((q * (z - mu) ** 2).sum() / q.sum())
    raw = '01' if 'raw01' in d else '02'
    rr = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    err = rr.RT.values * 60 - np.interp(rr.lrt, kx, ky)
    rs = 1.4826 * np.median(np.abs(err - np.median(err)))
    print(f'{d.split("/")[-1]:22s} hw {hw:5.1f} shift p50 {np.median(sh):+6.1f} | written {cal["rt"].get("spread_s", float("nan")):5.2f} | window-frame {s:5.2f} '
          f'(pi {pi:.3f} mu {mu:+.1f}) | DIA-NN robust {rs:5.2f} | cover 4*s {np.mean(np.abs(err) <= 4 * s):.3f} vs written {np.mean(np.abs(err) <= 4 * cal["rt"].get("spread_s", 0)):.3f}', flush=True)
