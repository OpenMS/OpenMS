"""Offline check: RT sigma from anchors whose acceptance does not depend on their RT residual.
An anchor is RT-independent if its score evidence plus the RT log-likelihood ratio at the window edge (the least
favourable residual) still reaches the weakest accepted evidence. Its residual is then not selected for, so its spread
is not truncated. Reports: sigma as written (all anchors), sigma on RT-independent anchors (1.4826 x median |residual|
to the written map), their count, and the share of DIA-NN 1 % IDs within 4 sigma for both (evaluation only).
Usage: sigma_check.py <cal dir> ..."""
import json, sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]
c = sqlite3.connect(f'file:{B}/pairlib-20261113/lib/library.pqp?mode=ro', uri=True)
lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])
for d in sys.argv[1:]:
    cal = json.load(open(f'{d}/cal.json')); raw = '01' if 'raw01' in d else '02'
    s1, hw = cal['stage1']['sigma_s'], cal['stage3']['half_window_s']
    kx, ky = np.array(cal['rt']['knots_library_rt']), np.array(cal['rt']['knots_run_rt_s'])
    t = pd.read_csv(f'{d}/diag.tsv', sep='\t', comment='#')
    t = t[t.stage3_score.notna()]
    dec = np.sort(t[t.decoy == 1].stage3_score.values)
    sev = -np.log((len(dec) - np.searchsorted(dec, t.stage3_score.values, side='left') + 1) / (len(dec) + 1))
    llr = lambda r: -0.5 * (r / s1) ** 2 - np.log(s1 * np.sqrt(2 * np.pi)) + np.log(2 * hw)
    ev = sev + llr(t.stage3_rt.values - t.stage3_centre_rt.values)
    acc = (t.decoy.values == 0) & (t.stage3_q.values <= 0.01)
    estar = ev[acc].min()
    indep = acc & (sev + llr(hw) >= estar)
    res = t.stage3_rt.values - np.interp(t.library_rt.values, kx, ky)
    sig_all = cal['rt']['sigma_s']
    sig_ind = 1.4826 * np.median(np.abs(res[indep])) if indep.sum() >= 20 else np.nan
    r = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    e = np.abs(np.interp(r.lrt, kx, ky) - r.RT.values * 60)
    print(f'{d:20s} anchors {acc.sum():4d} RT-independent {indep.sum():4d} | sigma written {sig_all:5.2f} cover {np.mean(e <= 4 * sig_all):.3f} | '
          f'sigma indep {sig_ind:5.2f} cover {np.mean(e <= 4 * sig_ind):.3f}')
