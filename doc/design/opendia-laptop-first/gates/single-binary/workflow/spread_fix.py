"""Offline check of the RT spread fit (evaluation against DIA-NN only for coverage). The C++ fit takes the target
residuals to the final map but the decoy background residuals to the stage-3 window centre (stage-1/2 map). When the two
maps differ, the background's window edges no longer line up with the targets' and the Gaussian absorbs the mismatch.
Fix: decoy residuals to the final map at their target's library RT, cut to the same |r| <= half window.
Usage: spread_fix.py <library.pqp> <cal dir> ... (cal dir with cal.json and diag.tsv)"""
import json, sys, sqlite3, numpy as np, pandas as pd
LIB, DIRS = sys.argv[1], sys.argv[2:]
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]


def mixture(x, bgx, hw):
    e = np.linspace(-hw, hw, 41); h, _ = np.histogram(bgx, e); h = (h + 1) / (h.sum() + len(h)) / np.diff(e)
    bg = h[np.clip(np.searchsorted(e, x) - 1, 0, len(h) - 1)]; pi, mu, s = 0.3, 0.0, hw / 4
    for _ in range(200):
        g = np.exp(-0.5 * ((x - mu) / s) ** 2) / (s * np.sqrt(2 * np.pi)); q = pi * g / (pi * g + (1 - pi) * bg)
        pi, mu = q.mean(), (q * x).sum() / q.sum(); s = np.sqrt((q * (x - mu) ** 2).sum() / q.sum())
    return pi, mu, s


for d in DIRS:
    lib_path = LIB if not d.startswith('lib=') else None
    c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
    pre = pd.read_sql("select ID precursor_id, TRAML_ID traml, DECOY decoy, LIBRARY_RT lrt from PRECURSOR", c)
    tl = pre[pre.decoy == 0].set_index('traml').lrt
    pre['target_lrt'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tl), pre.lrt)
    lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
      join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
    lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])
    cal = json.load(open(f'{d}/cal.json')); kx, ky = np.array(cal['rt']['knots_library_rt']), np.array(cal['rt']['knots_run_rt_s'])
    hw = cal['stage3']['half_window_s']
    t = pd.read_csv(f'{d}/diag.tsv', sep='\t', comment='#'); t = t[t.stage3_score.notna()]
    t = t.merge(pre[['precursor_id', 'target_lrt']], on='precursor_id', how='left')
    r = (t.stage3_rt - np.interp(t.target_lrt, kx, ky)).values; tg = (t.decoy == 0).values
    xt = r[tg]; xt = xt[np.abs(xt) <= hw]
    old = mixture(xt, (t.stage3_rt - t.stage3_centre_rt).values[~tg], hw)
    xd = r[~tg]; xd = xd[np.abs(xd) <= hw]
    new = mixture(xt, xd, hw)
    raw = '01' if 'raw01' in d else '02'
    rr = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    e = rr.RT.values * 60 - np.interp(rr.lrt, kx, ky)
    rs = 1.4826 * np.median(np.abs(e - np.median(e)))
    print(f'{d.split("/")[-1]:22s} windowed {cal["stage3"].get("windowed_second_pass")} hw {hw:5.1f} | written spread {cal["rt"].get("spread_s", float("nan")):5.2f} '
          f'| offline old {old[2]:5.2f} new {new[2]:5.2f} (pi {new[0]:.3f}) | DIA-NN robust {rs:5.2f} | cover 4*new {np.mean(np.abs(e) <= 4 * new[2]):.3f} '
          f'4*old {np.mean(np.abs(e) <= 4 * old[2]):.3f}', flush=True)
