# R-WP2: Identification-free run calibration

[Master plan](../README.md)

**Objective.** Map library RT, fragment m/z and precursor 1/K0 to each run and estimate extraction tolerances, without IDs or iRT peptides. Modern runs often have no iRT kit.

**Delivered:** `OpenDIACalibrate` (`f5ffc1f578`, frozen design), with the CNN and search wiring in `63d354a6a8`, `49c08f097c` and `3671d01f94`.
- **Stage 1** (after DIATraceExtractor): whole-run co-elution anchors chosen by paired target-decoy competition, then an RT-mode filter.
- **Stage 2:** the v0.13 latent peak/no-peak mixture over 11 monotone knots, with fixed π, β and σ.
- **Stage 3:** anchors at paired q ≤ 0.01 within the stage-2 window.
- **Stage 4:** the RT map, the MS2 ppm offset and a robust 1/K0 line, with tolerances of 4 σ, 5 σ and 4 σ.

| HeLa (design run) | Value |
|---|---|
| Anchors | 1278 |
| RT error vs 19132 teacher IDs | p50 1.8 s, p95 7.3 s; 96.3% within ±8.8 s |
| 1/K0 offset | +0.0216 (teacher +0.0203) |
| MS2 | +0.35 ppm, ±10.7 ppm |

**Known failures and differences:**
- On the independent 500 pg run it failed closed with 0 anchors. Stage 1 kept only 36 anchors, and that run's library RTs predict poorly (teacher σ 46.8 s).
- Stage 2 uses penalised least squares plus an isotonic fit, not v0.13's SLSQP, so numerical parity is not shown.
- The stage-2 prior is the stage-1 map, not v0.13's fixed linear prior.

**Open gates:**
1. Recover a calibration on low-input runs without tuning on an evaluation run. This needs a new held-out acquisition.
2. A v0.13 parity fixture for stage 2.
3. MS1 calibration, which has not been started.

## Panel scaling and escalation, 7 November 2026 (HeLa)

Anchors grow linearly with the panel, at about 0.13 per target. Panel 2500 gives 316, 10000 gives 1278 and 40000 gives 5479, and 95.8 to 96.9% of teacher IDs fall within the RT tolerance.
Panels of 500 and 1000 fail on the (D+1)/T floor, since q ≤ 0.01 needs at least 100 winning targets.
Memory grows by about 75 KB per pair: 4.0 GB and 57 s at 40000.

`-max_panel_targets <n>` (`da5b13983b`, default off, byte-identical) retries a failed calibration with a 4× panel.
The gate passed: 1000 → 4000, 477 anchors, 97.6% coverage.
This fixes failures caused by too few anchors. It does not fix weak target-decoy separation (500 pg: 26 decoys in the top 100), which has not been retried, by design.
A full-library panel needs streaming or a two-pass panel to bound memory.

The stage-2 EM already runs over every candidate peak of every panel target, but it is seeded by the stage-1 map.
A seed-free EM over many local starts would remove the stage-1 dependency. It is not yet designed, and needs a low-signal design run that is not 500 pg.
