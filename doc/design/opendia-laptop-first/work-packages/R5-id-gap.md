# R-WP5: Close the gap to OpenSWATH + PyProphet

[Master plan](../README.md)

**Objective.** Explain and close the gap between 11830 rewrite precursors and 19132 SVM / 19725 XGBoost precursors at 1% on HeLa.

**Method:** a loss funnel over the teacher's q ≤ 0.01 IDs. At each step, count how many teacher IDs survive:
1. Library and decoy pairing eligibility.
2. A calibrated-search candidate exists.
3. A candidate apex lies inside the teacher peak bounds.
4. The precursor is in the CNN scoring manifest.
5. The CNN-selected window is the teacher peak.
6. The precursor wins its paired competition and passes q ≤ 0.01.

For each loss, record the size, an RT/IM/intensity/support stratification, and whether the loss is in candidates, choice or separation.
Fix the cause with a predeclared HeLa gate, then reapply unchanged to an independent run.

**State:** active.

## First funnel, 7 November 2026 (HeLa design run, diagnostic only)

Scripts and logs: `OpenMS-r1-benchmark-20261107/raw-trace-cnn-hela/gap-funnel-20261107/`.

| Where SVM q ≤ 0.01 IDs (19132) end up | Count | Share |
|---|---|---|
| Identified by the rewrite | 9912 | 51.8% |
| Not in library or no decoy pair | 343 | 1.8% |
| No candidate apex inside the teacher peak bounds | 1612 | 8.4% |
| Wrong window chosen | 153 | 0.8% |
| **Right window, score below the 1% threshold** | **6600** | **34.5%** |
| Right window, lost to its paired decoy | 511 | 2.7% |

XGBoost shows the same pattern. The main loss is score separation, not candidates, windows or calibration.
The below-threshold group has weaker teacher scores and about half the teacher area, but its RT error matches the identified group.
Its 1/K0 error is a little larger: median 0.0146 against 0.0126.

**Predeclared rescoring diagnostic** (`GATE_rescore.md`, rule: more than 1000 extra targets means information is missing from the inputs).
A cross-fitted, decoy-negative logistic model was fitted on the CNN logit plus seven window features the CNN does not see: |RT error|, |1/K0 error|, apex intensity, area, transition support, matched scans and width.
It gave 14695 targets at paired q ≤ 0.01. The logit-only control gave 11830. That is +2865, so the rule is met.
Next step: add these features to the native scorer, opt-in and behind a gate.
Not yet checked: baseline overlap of the gained IDs and an entrapment check. This is a same-run, teacher-seeded diagnostic, not a result.
