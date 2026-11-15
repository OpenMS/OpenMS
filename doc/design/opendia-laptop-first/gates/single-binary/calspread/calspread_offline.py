"""PLAN_CALSPREAD offline evaluation (frozen design; DIA-NN used only for criteria 1-2).
Window-frame mixture: target x_m = stage3_rt - stage3_centre_rt in [-hw, hw]; true peaks at s_m + N(mu, sigma) with
s_m = final map - centre at the target's library RT; background = decoy winners' window-frame histogram (40 bins, +1);
anchors (targets with stage3_q <= anchor FDR) excluded; start sigma = stage-1 sigma, mu 0, pi 0.3; EM to convergence;
accepted if pi n >= 20 and sigma < hw, else the map sigma. Only members inside the stage-1 map support are used, as in
C++ (approximated here by excluding targets whose window was extended to the run edge: |x| > hw).
Usage: calspread_offline.py <library.pqp> <cal dir> ..."""
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
    hw = cal['stage3']['half_window_s']; s1 = cal['stage1']['sigma_s']; afdr = cal['options']['anchor_fdr']
    t = pd.read_csv(f'{d}/diag.tsv', sep='\t', comment='#'); t = t[t.stage3_score.notna()].merge(pre[['precursor_id', 'target_lrt']], on='precursor_id')
    x = (t.stage3_rt - t.stage3_centre_rt).values; tg = (t.decoy == 0).values
    s = (np.interp(t.target_lrt, kx, ky) - t.stage3_centre_rt).values
    inside = np.abs(x) <= hw
    anchor = tg & (t.stage3_q.values <= afdr)
    use = tg & inside & ~anchor
    e = np.linspace(-hw, hw, 41); h, _ = np.histogram(x[~tg & inside], e); h = (h + 1) / (h.sum() + len(h)) / np.diff(e)
    xt, st = x[use], s[use]; bg = h[np.clip(np.searchsorted(e, xt) - 1, 0, len(h) - 1)]
    pi, mu, sig = 0.3, 0.0, s1
    for it in range(1000):
        z = xt - st
        g = np.exp(-0.5 * ((z - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi)); q = pi * g / (pi * g + (1 - pi) * bg)
        npi, nmu = q.mean(), (q * z).sum() / q.sum(); nsig = np.sqrt((q * (z - nmu) ** 2).sum() / q.sum())
        done = abs(nsig - sig) < 1e-4 and abs(nmu - mu) < 1e-4
        pi, mu, sig = npi, nmu, nsig
        if done: break
    ok = pi * len(xt) >= 20 and 0 < sig < hw
    spread = sig if ok else cal['rt']['sigma_s']
    raw = '01' if 'raw01' in d else '02'
    rr = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    err = rr.RT.values * 60 - np.interp(rr.lrt, kx, ky)
    rs = 1.4826 * np.median(np.abs(err - np.median(err)))
    cover = np.mean(np.abs(err) <= 4 * spread)
    c1 = abs(spread / rs - 1) <= 0.15; c2 = cover >= 0.95
    print(f'{d.split("/")[-1]:20s} n {len(xt):6d} anchors_excl {int((anchor & inside).sum()):4d} it {it:4d} | sigma {spread:5.2f} (pi {pi:.3f} mu {mu:+5.1f}, {"fit" if ok else "fallback"}) '
          f'| written {cal["rt"].get("spread_s", float("nan")):5.2f} | DIA-NN {rs:5.2f} ratio {spread / rs:.3f} {"ok" if c1 else "FAIL"} | cover {cover:.3f} {"ok" if c2 else "FAIL"}', flush=True)
