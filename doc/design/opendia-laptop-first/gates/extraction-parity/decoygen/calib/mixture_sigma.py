"""Offline: selection-free RT spread and offset. All stage-3 target winners (not only anchors) give residuals to the
final map; their density over the stage-3 window is a mixture pi N(mu, s) + (1 - pi) background, the background shape
taken from the decoy winners' residuals (histogram, same window) so it carries the window and edge effects. EM fits
pi, mu, s. Compare s and mu to DIA-NN's robust sigma and bias (evaluation only). Usage: mixture_sigma.py <cal dir> ..."""
import json, sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]
c = sqlite3.connect(f'file:{B}/pairlib-20261113/lib/library.pqp?mode=ro', uri=True)
lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])


def mixture(x, bg_x, hw, iters=200):
    edges = np.linspace(-hw, hw, 41)
    h, _ = np.histogram(bg_x, edges); h = (h + 1) / (h.sum() + len(h)) / np.diff(edges)
    bg = h[np.clip(np.searchsorted(edges, x) - 1, 0, len(h) - 1)]
    pi, mu, s = 0.3, 0.0, hw / 4
    for _ in range(iters):
        g = np.exp(-0.5 * ((x - mu) / s) ** 2) / (s * np.sqrt(2 * np.pi))
        r = pi * g / (pi * g + (1 - pi) * bg)
        pi, mu = r.mean(), (r * x).sum() / r.sum()
        s = np.sqrt((r * (x - mu) ** 2).sum() / r.sum())
    return pi, mu, s


for d in sys.argv[1:]:
    raw = '01' if 'raw01' in d else '02'
    cal = json.load(open(f'{d}/cal.json')); kx, ky = cal['rt']['knots_library_rt'], cal['rt']['knots_run_rt_s']
    hw = cal['stage3']['half_window_s']
    t = pd.read_csv(f'{d}/diag.tsv', sep='\t', comment='#'); t = t[t.stage3_score.notna()]
    lr = t.library_rt.values
    # Decoys' library_rt is their own; use the target's (the window centre) via the stored centre and map difference.
    res = t.stage3_rt.values - t.stage3_centre_rt.values  # residual to the window centre (stage-2 / stage-1 map)
    shift = np.interp(lr, kx, ky) - t.stage3_centre_rt.values  # final map minus centre, at the member's own library RT
    tg = t.decoy.values == 0
    x = (res - shift)[tg]; x = x[np.abs(x) <= hw]
    bgx = res[~tg]
    pi, mu, s = mixture(x, bgx, hw)
    r = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    e = r.RT.values * 60 - np.interp(r.lrt, kx, ky)
    rs = 1.4826 * np.median(np.abs(e - np.median(e)))
    print(f'{d:20s} written sigma {cal["rt"]["sigma_s"]:5.2f} | mixture pi {pi:.3f} mu {mu:+5.2f} s {s:5.2f} | DIA-NN sigma {rs:5.2f} bias {np.median(e):+5.2f} '
          f'| cover 4*s {np.mean(np.abs(e) <= 4 * s):.3f}, 4*s about mu {np.mean(np.abs(e - mu) <= 4 * s):.3f}')
