# Voxel v4 phase 1 result: window ranking term (R) and MS1 isotope features (M) for the single decoy-whitened network

Plan: `PLAN_V4.md` (sha256 `d0439a2bf6…`, written 2026-11-11T10:09:55Z, sha check OK after the run). Notes: `PLAN_V4.notes` (P1, written before training). Data: ProteoBench Astral raws 01 and 02. Entrapment labels are used for evaluation only. Training: one process, 10 torch threads, 14 min (raw 01) and 10 min (raw 02) for all three arms.

## Verdict: no arm passes

| Arm | Rule 1 (PLAN_V2 1–4 against BASE) | Rule 2 (N@2 % above D11, both raws) | Rule 3 (null-winner AUC ≤ D11 + 0.005) | Result |
|---|---|---|---|---|
| D11R (ranking term) | FAIL: ΔN@2 % on raw 02 is −85..+258 | FAIL (1523 / 1451 against 1572 / 1572) | PASS (0.5022 / 0.5037) | FAIL |
| D11M (MS1 features) | FAIL: FDP_c at q ≤ 0.01 on raw 01 is 0.016 | FAIL (1592 / 1550) | PASS (0.5052 / 0.5029, limit 0.5065 / 0.5072) | FAIL |
| D11MR (both) | FAIL: FDP_c at q ≤ 0.01 is 0.017 on raw 01 | FAIL on raw 01 (1573 against 1572 is above by 1, raw 02 1637); see note | PASS (0.5051 / 0.5036) | FAIL |

Note on D11MR rule 2: 1573 > 1572 and 1637 > 1572 literally pass rule 2. The arm fails on rule 1 (FDP_c 0.017 at q ≤ 0.01 on raw 01, limit 0.015). Paired intervals against D11 include 0 for every arm and raw (table below). D11 stays the best single network.

## Numbers (precursor target winners after pooled paired TDC; FDP_c = N_e (1 + 1/r) / N; decoy-count estimates)

| Raw | Arm | N@0.5 % | N@1 % | N@1.5 % | N@2 % | q ≤ 0.01 N (FDP_c) | Entrapment AUC | Null-winner AUC |
|---|---|---|---|---|---|---|---|---|
| 01 | BASE | 900 | 1024 | 1267 | 1314 | 1081 (0.013) | 0.531 | 0.518 |
| 01 | D11 | 1256 | 1348 | 1496 | 1572 | 1504 (0.015) | 0.511 | 0.5015 |
| 01 | D11R | 1091 | 1370 | 1448 | 1523 | 1463 (0.015) | 0.515 | 0.5022 |
| 01 | D11M | 1241 | 1418 | 1534 | 1592 | 1544 (0.016) | 0.518 | 0.5052 |
| 01 | D11MR | 1217 | 1439 | 1494 | 1573 | 1535 (0.017) | 0.517 | 0.5051 |
| 02 | BASE | 1215 | 1288 | 1333 | 1376 | 1241 (0.008) | 0.535 | 0.523 |
| 02 | D11 | 980 | 1334 | 1417 | 1572 | 1353 (0.012) | 0.514 | 0.5022 |
| 02 | D11R | 862 | 1267 | 1370 | 1451 | 1345 (0.013) | 0.518 | 0.5037 |
| 02 | D11M | 1285 | 1338 | 1440 | 1550 | 1408 (0.013) | 0.516 | 0.5029 |
| 02 | D11MR | 1244 | 1309 | 1513 | 1637 | 1459 (0.014) | 0.517 | 0.5036 |

Paired bootstrap against D11 (1,000 precursor resamples, seed 5, 90 % interval; `VS_D11.txt`):

| Raw | Arm − D11 | ΔN@1 % | ΔN@2 % |
|---|---|---|---|
| 01 | D11R | +22 (−132..+55) | −49 (−108..+33) |
| 01 | D11M | +70 (−113..+169) | +20 (−56..+135) |
| 01 | D11MR | +91 (−98..+186) | +1 (−87..+107) |
| 02 | D11R | −67 (−261..+270) | −121 (−182..+27) |
| 02 | D11M | +4 (−113..+446) | −22 (−134..+197) |
| 02 | D11MR | −25 (−147..+419) | +65 (−75..+191) |

Funnel against DIA-NN 1 % (human, exact library match, reference only; `funnel_<arm>.txt`, stages as in `../v3-seq/MISSES_D11.md`):

| Raw | Arm | C wrong window | D right window, below cut | E accepted |
|---|---|---|---|---|
| 01 | D11 | 421 | 681 | 1290 |
| 01 | D11R | 441 | 641 | 1310 |
| 01 | D11M | 389 | 644 | 1359 |
| 01 | D11MR | 416 | 613 | 1363 |
| 02 | D11 | 383 | 744 | 1273 |
| 02 | D11R | 384 | 799 | 1217 |
| 02 | D11M | 372 | 746 | 1282 |
| 02 | D11MR | 395 | 744 | 1261 |

## Observations (not part of the verdict)

1. **The ranking term had almost nothing to learn.** It acts on positives only. In round 0 the positives are the anchors, and D11 already ranks the anchor window first for 99.2 % / 99.5 % of anchor precursors. Only 1 / 0 stage-C precursors are anchors. In rounds 1–2 the positives are the model's own q-accepted winners, which are its top windows by construction. The term therefore cannot reach the stage-C precursors, which sit below the accepted set. Stage C did not drop (441 / 384 against 421 / 383).
2. **D11 chooses windows worse than a linear model on the same features** (`PLAN_V4.notes` P1, before training). For precursors with a window within 10 s of the DIA-NN RT, the winner is the right window for D11 1971 / 2017, the linear decoy-whitened LDA (L11d) 2266 / 2301, and BASE 2278 / 2306. D11 already has the linear features in its embedding (h = [z, tanh(g(z))]), yet it chooses the right window less often than L11d. Why is not shown. One candidate is that the readout weights the tanh units, whose values saturate below the top 1 % (`MISSES_D11.md`). The readout weights were not inspected. A test (for example the window-choice rate of the z part of the D11 readout alone) would need a new plan.
3. **MS1 features help a little in the funnel but not at the cut.** D11M lowers stage C by 32 / 11 and raises E by 69 / 9, but FDP_c at q ≤ 0.01 rises on raw 01 (0.016) and N@2 % on raw 02 falls by 22. Entrapment AUC rises by 0.002–0.007 for every arm, which is inside the base limit but in the wrong direction.

## Next

- Phase 2 (planned in `PLAN_V4.md`): 12-fragment features from the DIALibGen 0.13.1 library, extracted for the existing windows. The library and the export are in progress (`astral-lib12-20261111/`). Phase 2 is frozen in its own plan before its runs. Its base stays D11, because no phase 1 arm passed.
