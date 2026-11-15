"""PLAN_CALREFINE criteria 1-2 (DIA-NN used only here): spread_s vs DIA-NN robust spread around the written map, and the
share of DIA-NN 1 % IDs within the written tolerance, for the refined calibration and the current-code baseline.
Usage: calrefine_eval.py <library.pqp> <raw 01|02> <baseline cal.json> <refined cal.json>"""
import json, sys, sqlite3, numpy as np, pandas as pd
LIB, RAW, BASE, NEW = sys.argv[1:5]
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])
rr = ref[ref.Run.str.endswith('_' + RAW)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
out = {}
for name, path in (('baseline', BASE), ('refined', NEW)):
    j = json.load(open(path)); rt = j['rt']
    err = rr.RT.values * 60 - np.interp(rr.lrt, rt['knots_library_rt'], rt['knots_run_rt_s'])
    rs = 1.4826 * np.median(np.abs(err - np.median(err)))
    cover = np.mean(np.abs(err) <= rt['tolerance_half_width_s'])
    out[name] = cover
    print(f'{NEW.split("/")[-2]:14s} {name:8s} spread {rt.get("spread_s", float("nan")):6.2f} DIA-NN {rs:5.2f} ratio {rt.get("spread_s", float("nan")) / rs:.3f} '
          f'tol {rt["tolerance_half_width_s"]:6.1f} cover {cover:.4f} |err| p50 {np.median(np.abs(err)):5.1f} anchors {rt["anchors_used"]} '
          f'refined {j["stage3"].get("refined_on_final_map", False)}', flush=True)
j = json.load(open(NEW)); rt = j['rt']
err = rr.RT.values * 60 - np.interp(rr.lrt, rt['knots_library_rt'], rt['knots_run_rt_s'])
rs = 1.4826 * np.median(np.abs(err - np.median(err)))
c1 = abs(rt['spread_s'] / rs - 1) <= 0.15
c2 = out['refined'] >= 0.95 and out['refined'] >= out['baseline'] - 0.01
print(f'{NEW.split("/")[-2]:14s} criterion 1 {"pass" if c1 else "FAIL"} | criterion 2 {"pass" if c2 else "FAIL"}')
