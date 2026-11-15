# PLAN_CALREFINE2: five-seed end-to-end re-test of the refined calibration

Written 15 November 2026, before any run. PLAN_CALREFINE (frozen) failed on a single-seed end-to-end criterion. Its
calibration criteria passed on all 4 runs. This plan re-tests only that end-to-end criterion with five seeds. The code
under test is unchanged (`calrefine/calrefine.patch`).

## Arms (arm-E library, both raws, OpenDIAWorkflow, seeds 0-4)
- **CUR:** the current-code calibration. The calibration is deterministic: the self-calibrated JSON is byte-equal to a
  standalone OpenDIACalibrate run. It is reused through `-calibration/-anchors`, from OpenDIACalibrate output of the
  current code.
- **REF:** the refined calibration, from the 4-run check of PLAN_CALREFINE; also deterministic and byte-equal to the
  workflow run. It is reused the same way.

The calibration is reused rather than recomputed, so seeds change only the scorer. The calibration step itself was
shown to be identical.

## Criteria (labels only here; evaluated with eval_run.py)
1. **Median N@FDP_c 1 %:** REF ≥ CUR − 1 % on each raw.
2. **FDP_c at q ≤ 0.01:** ≤ 0.010 for every REF seed on both raws.

The refined calibration is adopted if both hold. Otherwise nothing changes. There is no retuning.
