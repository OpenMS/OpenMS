# R-WP3: Raw-trace CNN scorer

[Master plan](../README.md)

**Objective.** Score candidate windows from raw fragment XIC/XIM tensors, in place of the handcrafted OpenSWATH feature suite.

**Delivered:**
- Python prototype and parity checks (`a7ce24802e`).
- Native LibTorch training and scoring in `OpenDIACNN` (`595807253b`).
- `OpenDIATraces` builds the tensors, with a configurable memory budget (`bf11435231`).
- Coordinate shifts from the calibration, which gave +339 precursors on HeLa (`49c08f097c`).

**Open gates:**
1. Seeds come from the same-run OpenSWATH teacher (SVM q ≤ 0.001, 3000 targets plus 3000 decoys). Replace them with a teacher-free or cross-run seed source.
2. Trace geometry (±15 s, ±0.03 1/K0, 0.02 Da) does not use the calibrated tolerances.
3. No MS1 or precursor features yet.

## Per-run training from calibration anchors, 7 November 2026 (HeLa design run)

Proposal: train the CNN for each run from the ID-free calibration, with no pretrained model and no teacher.
Each run was a single variable change against the 11830 teacher-seeded crossfit. The rules were in `raw-trace-cnn-hela/anchor-seeds-20261107/GATE.md` before running, and the pass mark was ≥ 10647 (90%).

| Seeds | Seed pairs | Targets at paired q ≤ 0.01 | In SVM |
|---|---|---|---|
| Teacher (OpenSWATH SVM q ≤ 0.001) | 3000 | 11830 | 85.9% |
| Calibration anchors (q ≤ 0.01) | 1186 | 9627 (fail) | 88.3% |
| Round 2: self-training on round-1 q ≤ 0.001 | 2251 | 7093 (fail) | n/a |

Crossfit took 1.5 min at 2.0 GB.

Interpretation:
- Teacher-free per-run training works, at 81% of the teacher-seeded count.
- Naive self-training makes it worse. The model's most confident targets are the easy peaks.
- The anchors are fewer than the teacher seeds and biased toward strong peaks.

Next candidates, each needing its own gate:
- More anchors, from a looser stage-3 cut or a larger panel.
- A PyProphet-style semi-supervised loop on the score itself rather than on reseeding.
- Running anchor extraction and CNN training inside OpenDIACalibrate.
