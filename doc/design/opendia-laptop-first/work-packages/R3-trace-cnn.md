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
