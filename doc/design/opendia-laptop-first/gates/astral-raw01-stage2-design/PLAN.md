# Stage-2 design run on raw 01 (written before the run)

Design run only, per `../PROTOCOL.md`. Only target and decoy labels are read; entrapment labels are never joined.

- **What runs:** the frozen raw 01 CNN step (`../pipeline/cnn_it0.sh` flags plus `-pooled_tdc`), with `-stage2` at its defaults:
  - 3 LDA rounds;
  - initial positive cut 0.15, later rounds 0.05;
  - core window ±4 s, MS1 window ±20 ppm.
- **What is read:**
  - target precursors at pooled paired q ≤ 0.01, against 321 for the same run without `-stage2`;
  - the per-fold LDA weights;
  - target and decoy counts above score cuts;
  - the per-window feature table.
- **Decision rule, fixed now:**
  - `-stage2` defaults are frozen as they are if all three hold:
    1. the run completes;
    2. every fold gives finite weights with a positive weight on `cnn_z`;
    3. the q ≤ 0.01 count is at least 321.
  - No parameter is tuned on this run. If any condition fails, the failure is recorded and the design is revisited, one variable at a time.
- **No ID claim:** a higher raw 01 count is not an ID claim. Whether the extra targets are real is only testable with entrapment, on raw 03 under a new hashed gate, which needs the user's approval.
