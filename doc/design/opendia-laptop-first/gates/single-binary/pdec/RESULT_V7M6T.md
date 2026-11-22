# RESULT_V7M6T: V7M6 on matched nodes (22 Nov 2026)

Gate: PRECHECK_V7M6T (`68e37bd097`, amended `01a6479896` before any job). Container `opendia-pd10-best.sif` (sha256 `1bfd94e4`). Four jobs, each holding a whole Intel Xeon Gold 6336Y node (`hl-codon-132-01` to `-04`), run in ABBA order after a cache warm-up. Evaluator `precheck_v7m6t.py`, output below.

**Checks.**
- **Sanity check:** H r7 def rep 1 is byte-identical to `entlib-gbt-20261122/runs/H_gbt_s0`.
- **Determinism:** rep 2 is byte-identical to rep 1 for all 8 arm-cells.
- **Timing noise:** within a cell, the two reps of an arm differ by at most 3 % (D r7 def: 1,921 s and 1,978 s).

## Result

| Cell | HeLa IDs def → V7M6 | r7 N @ FDP_c 1 % | Lower bound change | Wall (matched, 2 reps each) |
|---|---:|---:|---:|---:|
| D r0 | 28,174 → 32,539 (+15.5 %) | | | 1.18× |
| H r0 | 42,289 → 44,067 (+4.2 %) | | | 1.08× |
| D r7 | 21,950 → 24,683 | 21,468 → 24,324 (**+13.3 %**) | −0.05 pp | **1.504×** |
| H r7 | 38,009 → 39,570 | 37,339 → 38,888 (**+4.2 %**) | −0.05 pp | 1.16× |

**Rule B: FAIL.** The r0 HeLa, r7 N and lower-bound bars all pass, but the wall bar fails on D r7: 1.504× against ≤ 1.5×. By the rules this is a fail, and the bar is not moved.

**Reading.**
- **The GAP2 wall result was not node noise.** On matched nodes the D r7 ratio is 1.504×, essentially the 1.52× seen in GAP2. H r7 is 1.16× here, not 0.79×, so GAP2's H number was the noisy one.
- **The gain holds on the GBM default.** D r7 is +13.3 % and H r7 +4.2 %, which matches or beats GAP2 (+13.3 % and +3.6 %). The tool's own q is slightly tighter (lower bound −0.05 pp on both runs).
- **Where the extra time goes (D r7, rep 1, per pass):**

  | Stage | def | V7M6 | Ratio |
  |---|---:|---:|---:|
  | Candidate windows | 10.1 M | 14.1 M | 1.39× |
  | Extraction | 337 s | 494 s | 1.47× |
  | Interference features | 230 s | 561 s | **2.44×** |
  | Scoring | 286 s | 353 s | 1.23× |

  The interference step grows with the square of crowding, not linearly with windows. Over both passes it adds about 610 s of the roughly 980 s that V7M6 adds per run (62 %). If interference grew only linearly (1.39×), the D r7 wall would be near 1.27×. That is a diagnostic projection, not a measurement.
- **Peak memory:** D r7 rises from 34.4 to 46.6 GB (+36 %).

**What this means.** V7M6 is not adopted as it stands. The direct path is making the interference resolver near-linear in the number of candidate windows, with byte-identical output, and then re-running this same gate. That is a code change to the default path, so it needs a decision.

## Evaluator output

```
run lvl  arm     N  hela  wall1_s  wall2_s  rss_gb  rep_same  lower  combined  n_fdp1  d_hela  wall_x  d_fdp1  d_lower_pp
  D  r0  def 28174 28174   372.67   371.92  5.3017      True    NaN       NaN     NaN  0.0000  1.0000     NaN         NaN
  D  r0 V7M6 32539 32539   435.80   439.94  7.0111      True    NaN       NaN     NaN  0.1549  1.1761     NaN         NaN
  D  r7  def 22204 21950  1921.28  1978.03 34.3817      True 0.0114    0.0131 21468.0  0.0000  1.0000  0.0000      0.0000
  D  r7 V7M6 24956 24683  2937.77  2928.32 46.6054      True 0.0109    0.0125 24324.0  0.1245  1.5044  0.1330     -0.0500
  H  r0  def 42289 42289   593.00   593.71 13.3106      True    NaN       NaN     NaN  0.0000  1.0000     NaN         NaN
  H  r0 V7M6 44067 44067   641.86   637.35 13.3076      True    NaN       NaN     NaN  0.0420  1.0779     NaN         NaN
  H  r7  def 38507 38009  2385.83  2341.52 27.9989      True 0.0129    0.0148 37339.0  0.0000  1.0000  0.0000      0.0000
  H  r7 V7M6 40069 39570  2737.67  2743.47 31.0293      True 0.0125    0.0142 38888.0  0.0411  1.1595  0.0415     -0.0479
SANITY H r7 def vs entlib-gbt H_gbt_s0: SAME
RULE B V7M6: r0 HeLa>=+3% True, r7 FDP_c>=+2% True, lower<=+0.2pp True, matched wall<=1.5x False -> FAIL; max wall_x 1.504; mean r7 gain +0.0873
```
