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

## Independent run 2: ProteoBench Astral entrapment, raw 01 (7 November 2026)

The frozen teacher-free design (`4cbd7701bb`, panel 40000) was applied once, under a gate written before the run: `astral-entrapment-20261107/GATE.md`, with its sha256 recorded.
- The Thermo `.raw` streamed natively; OpenDIACalibrate took 2:04 at 10.4 GB.
- Calibration **failed closed with 0 anchors**: stage 1 kept 95, with σ 18.8 s. The stage-3 window was ±75 s, and 19 decoys were in the top 100.
- No IDs were produced, so the entrapment FDP could not be measured. The gate is not passed.

This is the second independent run, after 500 pg, on which the calibration finds too few separable anchors. The HeLa design does not transfer.
Both runs have low identifiable populations: about 1500 OpenSWATH result rows per raw on Astral, and a low-load 500 pg sample. The Astral run also has no 1/K0 and uses predicted, normalised library RTs.
Calibration robustness on weak runs is now the blocking problem for the whole teacher-free path. It needs a new design, on a design run that is not an evaluation run.

### Option B: raw 01 as the design run (user decision, 7 November 2026)

Raw 01 is now the design run. Entrapment labels are never joined to its outputs. Raw 02 is the single evaluation and raw 03 is held back. See `astral-design-raw01-20261107/PROTOCOL.md` and `DIAGNOSIS.md`.
- **Diagnosis** (new `-out_diagnostics`): a power problem.
  - Only 11264 of 40000 panel targets have a peak, and no stage reaches q ≤ 0.01.
  - The stage-2 EM map is nearly flat: 109 s median distance to the stage-1 anchors, against 20 s for the stage-1 map.
  - Centring stage 3 on the stage-1 map alone gives 105 targets at q ≤ 0.05 but still 0 at q ≤ 0.01.
- **Fix** (`290051546d`, opt-in `-window_panel_targets <n>`): when stage 3 fails after a successful stage 1, one more raw pass collects up to `n` targets.
  - It stores only spectra inside the stage-3 windows around the stage-1 map.
  - Stages 3 and 4 then repeat unchanged.
  - On raw 01 at n = 300000: 512 anchors, 4 decoys above the cut, RT σ 21 s, 3:37 and 16.5 GB.
  - HeLa p40000 output is byte-identical with and without the option.
- **Rest of the raw 01 design run** (no entrapment labels, `pipeline/RESULT.md`):
  - Seeds: 509 positive and 509 decoy.
  - Per-fold paired TDC: 122 targets at q ≤ 0.01, all from one fold. The other folds sit on the (D+1)/T floor.
  - `-pooled_tdc` (`0eccb5f039`, opt-in) standardises each fold's logits by its decoy winners' median and MAD, then runs one competition over all folds: 321 targets.
  - `-ss_iterations 1` gives 0, with 321 or 634 positives (overfit).
  - Frozen design: `-ss_iterations 0 -pooled_tdc` at `46be19510d`.
- **Raw 02 single evaluation** (`astral-entrapment-raw02-20261108/GATE.md` and `RESULT.md`, gate hashed and pushed in `914ad5e14e` before the run): **FAIL**.
  - Calibration transfers: 604 anchors, 0.99% of them entrapment.
  - 360 targets at q ≤ 0.01, of which 5 are entrapment.
  - FDP_c = 0.028 > 0.015 (95% interval 0.009 to 0.065). The lower bound is 0.014.
  - At q ≤ 0.05: 635 targets and FDP_c 0.044.
  - All 5 entrapment hits score high, not at the margin. The experimental q underestimates the entrapment FDP on this run.
  - The failure is in score separation and error control (R-WP4), not in calibration.
  - Raw 02 is spent, and raw 03 is still held back.

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
