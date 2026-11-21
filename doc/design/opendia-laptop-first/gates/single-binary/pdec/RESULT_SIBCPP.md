# RESULT_SIBCPP: C++ charge-state sibling rescoring passes on r ≈ 1.56 (21 Nov 2026)

The gate is `PRECHECK_SIBCPP.md`, committed before any job ran.

- **Code:** OpenMS `a3eea6d9ee` adds `-sibling_features`, on top of `f3005a4b14`. I reviewed it, and 26/26 tests pass.
- **Container:** `opendia-pd10-sibf.sif`, sha256 `c67d7d9fa8cf0f127f8326357ee02ec8d34b14d6922e881250331ca8a23e28f6`.
- **Runs:** `sibcpp-20261121/runs`. All 8 exited 0.
- **Evaluator:** `eval_pdec.py`, sha256 `bb03387d…`. Full output is in `EVAL_SIBCPP.txt`.

**Default unchanged.** With the flag off, `precursors.tsv` is byte-identical to `bt156-20261121/<R>_s<seed>` on all 4 runs (D and H, seeds 0 and 1).

## Numbers (r ≈ 1.1225 for this library)

| Run | Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Gain |
|---|---|---|---:|---:|---:|---:|---:|---:|
| D | cur | 0 | 17,839 | 154 | 0.86 % | 1.63 % | 16,678 | |
| D | cur | 1 | 18,083 | 158 | 0.87 % | 1.65 % | 17,038 | |
| D | sib | 0 | 24,110 | 179 | 0.74 % | 1.40 % | 23,428 | +40.5 % |
| D | sib | 1 | 24,413 | 182 | 0.75 % | 1.41 % | 23,392 | +37.3 % |
| H | cur | 0 | 32,719 | 216 | 0.66 % | 1.25 % | 32,310 | |
| H | cur | 1 | 32,844 | 222 | 0.68 % | 1.28 % | 32,415 | |
| H | sib | 0 | 37,320 | 245 | 0.66 % | 1.24 % | 36,845 | +14.0 % |
| H | sib | 1 | 37,055 | 243 | 0.66 % | 1.24 % | 36,610 | +12.9 % |

**Pooled bounds and pair fairness** (P = probability the entrapment target beats its own decoy):

| Run | Arm | Lower | Combined | P |
|---|---|---:|---:|---:|
| D | cur | 0.87 % | 1.64 % | 0.523 |
| D | sib | 0.74 % | 1.41 % | 0.499 |
| H | cur | 0.67 % | 1.26 % | 0.511 |
| H | sib | 0.66 % | 1.24 % | 0.507 |

## Rule

| Criterion | Bar | Observed | |
|---|---|---|---|
| H mean gain | ≥ +5 % | +13.5 % | pass |
| D mean gain | ≥ 0 | +38.9 % | pass |
| Every H seed gains | yes | +14.0 %, +12.9 % | pass |
| Lower bound change | ≤ +0.2 pp | D −0.13 pp, H −0.01 pp | pass |

**PASS.** Freeze PLAN_SIB, the r ≈ 7 run on D seeds 0-4.

## Notes

- **Symmetry check.** Library-wide sibling counts are the same for targets and decoys. The share with 0, 1 or 2 siblings is 13.8 / 49.7 / 36.5 % for targets and 13.8 / 49.7 / 36.5 % for decoys. In D_sib_s0 windows, mean `s_n` (1.297 both) and `s_co` (0.470 vs 0.461) are equal. Only `s_best` differs (0.364 vs 0.124), which is the evidence the feature is meant to carry. Entrapment calibration improves, so the gain is not a target-only artifact.
- **Bigger than the Python estimate.** The C++ gain on D (+38.9 %) is far above the Python siblings-only arm (−0.3 / +3.5 %) and the Bts arm (+21-26 %). Two likely reasons:
  - the C++ version also rescores the RT-refine first pass, so the ridge RT fit trains on more IDs;
  - the network scorer uses the new features differently from the GBM.

  This is not decomposed here. The entrapment bounds are the arbiter.
- **Cost.** Wall time is within node noise of cur. D: sib 5:21 and 5:35, cur 7:41 and 4:39. H: sib 14:16 and 8:03, cur 11:25 and 7:28. Peak RSS is unchanged at 9.5-9.7 GB on D and 14.4 GB on H. The second scorer call is cheap next to extraction.
- **Fold validity.** Siblings share the family fold (0 of 245,711 multi-precursor peptides span two folds), so each first-pass score feeding `s_best` is out-of-fold for every precursor in the group.
