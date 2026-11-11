# Voxel v4 phase 2 result: 12-fragment evidence features (F) and linear-part window choice (W)

Plan: `PLAN_V4B.md` (sha256 `77bd7201e6…`, written 2026-11-11T11:35:23Z, sha check OK after the run). Data: ProteoBench Astral raws 01 and 02. Entrapment labels are used for evaluation only. DIA-NN is a reference only (funnel and window choice).

**Library and inputs:**
- 12-fragment library: DIALibGen 0.13.1 in generate mode (Astral, NCE 25, 14 threads), `astral-lib12-20261111/`.
- It covers all 300,223 original targets. 92 % of them have 12 fragments.
- Decoys get their target's predicted intensities, with m/z on the decoy sequence. The m/z code reproduces all 3.6 M original fragment m/z within 9e-5.
- Extraction: `-out_xipm` for the existing population windows only. The search and the candidate windows are unchanged.
- Check: the 6 original fragments in the new export reproduce `xic.npy` (median relative difference 0, p99 0.037 / 0.034).

## Verdict: D11F and D11FW pass; D11W fails

| Arm | Rule 1 (PLAN_V2 1–4 against BASE) | Rule 2 (N@2 % above D11, both raws) | Rule 3 (null-winner AUC ≤ D11 + 0.005) | Result |
|---|---|---|---|---|
| D11W (D11, linear-part window choice) | FAIL: FDP_c at q ≤ 0.01 on raw 01 is 0.0173 | FAIL (1560 / 1568) | PASS (0.4947 / 0.4951) | FAIL |
| D11F (+6 12-fragment features) | PASS: ΔN@2 % +288..+500 / +21..+415; ΔN@1 % −203..+651 / −738..+327; FDP_c 0.0136 / 0.0149; entrapment AUC 0.516 / 0.515 (BASE 0.531 / 0.535) | PASS (1688 / 1635) | PASS (0.5042 / 0.5040, limit 0.5065 / 0.5072) | **PASS** |
| D11FW (D11F, linear-part window choice) | PASS: ΔN@2 % +295..+484 / +18..+405; ΔN@1 % −214..+643 / −474..+305; FDP_c 0.0136 / 0.0150; entrapment AUC 0.491 / 0.496 | PASS (1714 / 1620) | PASS (0.4936 / 0.4934) | **PASS** |

FDP_c values at q ≤ 0.01 are exact, not rounded:

| Arm | raw 01 | raw 02 |
|---|---|---|
| D11F | 0.01360 (11 of 1630) | 0.01491 (11 of 1487) |
| D11FW | 0.01362 (11 of 1628) | 0.01497 (11 of 1481) |

Both arms are close to the 0.015 limit on raw 02.

## Numbers (precursor target winners after pooled paired TDC; FDP_c = N_e (1 + 1/r) / N; decoy-count estimates)

| Raw | Arm | N@0.5 % | N@1 % | N@1.5 % | N@2 % | q ≤ 0.01 N (FDP_c) | Entrapment AUC | Null-winner AUC | Right window (of 2392 / 2400) |
|---|---|---|---|---|---|---|---|---|---|
| 01 | BASE | 900 | 1024 | 1267 | 1314 | 1081 (0.013) | 0.531 | 0.518 | 2278 |
| 01 | D11 | 1256 | 1348 | 1496 | 1572 | 1504 (0.0147) | 0.511 | 0.5015 | 1971 |
| 01 | D11W | 1214 | 1325 | 1481 | 1560 | 1519 (0.0173) | 0.489 | 0.4947 | 2255 |
| 01 | D11F | 920 | 1525 | 1640 | 1688 | 1630 (0.0136) | 0.516 | 0.5042 | 2030 |
| 01 | D11FW | 924 | 1521 | 1637 | 1714 | 1628 (0.0136) | 0.491 | 0.4936 | 2279 |
| 02 | BASE | 1215 | 1288 | 1333 | 1376 | 1241 (0.008) | 0.535 | 0.523 | 2306 |
| 02 | D11 | 980 | 1334 | 1417 | 1572 | 1353 (0.0119) | 0.514 | 0.5022 | 2017 |
| 02 | D11W | 1061 | 1323 | 1483 | 1568 | 1355 (0.0104) | 0.497 | 0.4951 | 2286 |
| 02 | D11F | 880 | 1149 | 1539 | 1635 | 1487 (0.0149) | 0.515 | 0.5040 | 2093 |
| 02 | D11FW | 885 | 1115 | 1524 | 1620 | 1481 (0.0150) | 0.496 | 0.4934 | 2291 |

Paired bootstrap against D11 (1,000 precursor resamples, seed 5, 90 % interval; `VS_D11_B.txt`):

| Raw | Arm − D11 | ΔN@1 % | ΔN@2 % |
|---|---|---|---|
| 01 | D11W | −23 (−48..−4) | −12 (−25..+19) |
| 01 | D11F | +177 (−441..+293) | +116 (+30..+221) |
| 01 | D11FW | +173 (−450..+301) | +142 (+23..+222) |
| 02 | D11W | −11 (−13..+107) | −4 (−8..+17) |
| 02 | D11F | −185 (−440..+261) | +63 (−48..+221) |
| 02 | D11FW | −219 (−427..+220) | +48 (−57..+213) |

Funnel against DIA-NN 1 % (`funnel_<arm>.txt`):

| Raw | Arm | C wrong window | D right window, below cut | E accepted |
|---|---|---|---|---|
| 01 | D11 | 421 | 681 | 1290 |
| 01 | D11W | 137 | 985 | 1270 |
| 01 | D11F | 362 | 595 | 1435 |
| 01 | D11FW | 113 | 847 | 1432 |
| 02 | D11 | 383 | 744 | 1273 |
| 02 | D11W | 114 | 1024 | 1262 |
| 02 | D11F | 307 | 982 | 1111 |
| 02 | D11FW | 109 | 1213 | 1078 |

## Limits (read with the verdict)

1. **The gain is at 2 %, not at 1 %.** ΔN@2 % against D11 is above 0 on raw 01 for D11F and D11FW. On raw 02 the interval includes 0. At 1 % the D11F and D11FW curves are unstable: raw 02 N@1 % is 1149 / 1115 (D11 1334) while N@1.5 % is 1539 / 1524. The ΔN@1 % intervals are wide (about −450..+300). This is the same tail-shape instability already seen for D11 at 0.5 % on raw 02. The frozen rule does not test N@1 % against D11.
2. **FDP_c on raw 02 is at the limit** (0.0149 / 0.0150, 11 entrapment winners). One more entrapment winner at q ≤ 0.01 would fail rule 1.
3. **W moves windows but not IDs.** Linear-part window choice cuts stage C from about 400 to about 110 per raw. The windows it picks are right ones that score below the cut, so they move to stage D and E does not rise. W lowers entrapment AUC by about 0.02 (to 0.49–0.50), and null-winner AUC falls below 0.5. On D11 alone, W raises FDP_c at q ≤ 0.01 on raw 01 (0.0173).
4. **What the F features add.** `f12_libcos` (the library cosine over 12 fragments) separates stage D from the top 1 % of decoy winners with AUC 0.90 on raw 01 (pre-plan check, DIA-NN reference). D11F raises E by 145 on raw 01. On raw 02 it moves fewer precursors above the 1 % cut (1111 against 1273), although N@1.5 % and N@2 % are higher.
5. Only two raws. Raw 03 is not used (needs approval).

## Next

- D11F (or D11FW) is the first single network to beat D11 under a frozen rule. By PLAN_V4B it is the new base and the candidate for the C++ `OpenDIAScoringModel`.
- Moving it to C++ needs the 12-fragment library in the OpenDIA input, plus extraction of fragments that the search does not use. That is a design decision for the user (todo `sm-cpp`).
- The 1 % instability (limit 1) is the main open problem. One candidate is a per-fold tail calibration; it needs a new plan.
