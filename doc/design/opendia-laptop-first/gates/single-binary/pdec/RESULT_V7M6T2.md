# RESULT_V7M6T2: V7M6 passes rule B with the faster interference resolver (23 Nov 2026)

Gate: PRECHECK_V7M6T2 (`756063e1de`, pushed before any cell job). Container `opendia-pd10-fastint.sif` (sha256 `8de3a968`, OpenMS `955516d84c`). Same nodes as RESULT_V7M6T (`hl-codon-132-01` to `-04`, Xeon Gold 6336Y, whole node, ABBA). Evaluator `precheck_v7m6t.py`.

**Checks.**
- **Byte-identity precondition** (`fastint-20261122`): the new container is byte-identical on all 9 files of D r≈1.56 default, H r7 default and D r7 V7M6.
- **Sanity check:** H r7 def rep 1 is byte-identical to `entlib-gbt-20261122/runs/H_gbt_s0`.
- **Determinism:** rep 2 is byte-identical to rep 1 on all 8 arm-cells.
- **IDs:** identical to RESULT_V7M6T, as the byte-identity requires.

## Result

| Cell | HeLa IDs def → V7M6 | r7 N @ FDP_c 1 % | Lower bound change | Wall (V7M6T, old resolver) | Wall (now) |
|---|---:|---:|---:|---:|---:|
| D r0 | 28,174 → 32,539 (+15.5 %) | | | 1.18× | 1.18× |
| H r0 | 42,289 → 44,067 (+4.2 %) | | | 1.08× | 1.07× |
| D r7 | 21,950 → 24,683 | 21,468 → 24,324 (+13.3 %) | −0.05 pp | 1.504× | **1.375×** |
| H r7 | 38,009 → 39,570 | 37,339 → 38,888 (+4.2 %) | −0.05 pp | 1.16× | 1.13× |

**Rule B: PASS.** All four bars hold. The largest wall ratio is 1.375× (D r7), and the mean r7 gain is +8.7 %.

**Absolute wall times also fell.** On D r7, def went from 1,950 s to 1,894 s and V7M6 from 2,933 s to 2,604 s per run. The resolver per pass (rep 1, def / V7M6) went from 230 / 561 s to 181 / 402 s on the first pass. The remaining V7M6 cost is mostly the larger candidate set itself: extraction, more interference links, and scoring. Peak memory on D r7 is 47.0 GB against 34.2 GB (+37 %).

**Next.** Per the gate, a pass goes to a frozen 5-seed r7 plan, PLAN_V7M6, with D seeds 0-4 first and H after the D verdict.

## Evaluator output

```
run lvl  arm     N  hela  wall1_s  wall2_s  rss_gb  rep_same  lower  combined  n_fdp1  d_hela  wall_x  d_fdp1  d_lower_pp
  D  r0  def 28174 28174   372.58   371.72  5.2108      True    NaN       NaN     NaN  0.0000  1.0000     NaN         NaN
  D  r0 V7M6 32539 32539   439.99   438.60  7.1155      True    NaN       NaN     NaN  0.1549  1.1804     NaN         NaN
  D  r7  def 22204 21950  1894.52  1893.67 34.1705      True 0.0114    0.0131 21468.0  0.0000  1.0000  0.0000      0.0000
  D  r7 V7M6 24956 24683  2578.71  2628.40 47.0008      True 0.0109    0.0125 24324.0  0.1245  1.3746  0.1330     -0.0500
  H  r0  def 42289 42289   592.70   598.09 13.3095      True    NaN       NaN     NaN  0.0000  1.0000     NaN         NaN
  H  r0 V7M6 44067 44067   633.59   638.31 13.3052      True    NaN       NaN     NaN  0.0420  1.0681     NaN         NaN
  H  r7  def 38507 38009  2357.21  2406.86 28.3351      True 0.0129    0.0148 37339.0  0.0000  1.0000  0.0000      0.0000
  H  r7 V7M6 40069 39570  2677.87  2683.50 31.9036      True 0.0125    0.0142 38888.0  0.0411  1.1254  0.0415     -0.0479
SANITY H r7 def vs entlib-gbt H_gbt_s0: SAME
RULE B V7M6: r0 HeLa>=+3% True, r7 FDP_c>=+2% True, lower<=+0.2pp True, matched wall<=1.5x True -> PASS; max wall_x 1.375; mean r7 gain +0.0873
D_r7_def_1 interference_seconds=180.541983 interference_seconds=141.733659
D_r7_V7M6_1 interference_seconds=402.272655 interference_seconds=334.580408
```
