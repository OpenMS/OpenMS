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
| **Round 2b: PyProphet-style, round-1 q ≤ 0.01 targets vs best decoy windows** | 6000 | **11332 (pass)** | 86.6% |
| **Native `OpenDIACNN -ss_iterations 1`** (`6d02e94ac6`) | 6000 | **11118 (pass)** | n/a |
| **Panel-40000 anchors (3000 pairs) + `-ss_iterations 1`** (`da5b13983b`) | 3000 then 6000 | **13611 (pass)** | 84.7% |

Crossfit took 1.5 min at 2.0 GB.

Interpretation:
- Anchors alone reach 81% of the teacher-seeded count.
- Round 2 retrained on only the model's most confident targets and used easy decoy windows as negatives; it got worse.
- Round 2b iterates the way PyProphet does: targets at q ≤ 0.01 against each decoy's best-scoring window. That reaches 95.8% of the teacher-seeded count with no OpenSWATH input.
- Rounds 1, 2 and 2b were native crossfit runs with seeds rebuilt by a script in between.
- The loop is now in the tool as opt-in `-ss_iterations` (default off, byte-identical). Its subset order differs from the script's, giving 11118 against 11332.
- Only this design run was used; an independent run is open.

Next, each with its own gate:
- Entrapment check for the teacher-free path; 13611 exceeds the teacher-seeded count, which needs independent evidence.
- Move anchor training into OpenDIACalibrate.
- Cross-run global model (R-WP4).

## TOPAZ-style modality-flexible scorer, 8 November 2026 (raw 01 design, target/decoy labels only)

Code: `07ed9580e4` and `edc775d6bc`. Every option is opt-in, and the HeLa default `crossfit.tsv` is byte-identical after each rebuild. Gate record: [`gates/astral-raw01-cnn-topaz/`](../gates/astral-raw01-cnn-topaz/RESULT.md).

**The design** follows redeem TOPAZ (`add/topaz`: `conv_encoder.rs`, `coelution.rs`, `trace_input.rs`, `topaz.rs`).
- **Inputs:** 1D conv branches over a dual input view (max-normalised plus log1p) of:
  - the fragment XIC, with IM summed;
  - MS1 M+0, M+1 and M+2;
  - the fragment mobilogram over the core RT bins.
- **Co-elution head per branch:** lag-tolerant pairwise cosines, soft-apex spread, entropy and peakiness, plus a learned similarity embedding.
- **Scorer inputs:** MS1-vs-MS2 agreement, library fit, and an MLP.
- **Absent modalities:** MS1 and IM branch outputs are multiplied by `has_ms1` and `has_im`, so absent modalities add exact zeros (unit tested). One network serves Astral (no IM) and timsTOF.
- **Options:**
  - `-cnn_arch topaz`;
  - `-topaz_width` (32 is the TOPAZ layout);
  - `-weight_decay` (AdamW);
  - `-cnn_training mil`: masked-max MIL over each precursor's in-RT-window candidates, bag label target/decoy, with sampled bags per epoch.

**Raw 01, under the rule fixed in advance:** at least 732 targets at stage-2 q ≤ 0.01 (the reference is 717), plus a null-symmetry guard. Round 2 also had a compute budget of CNN time ≤ 1.5× grid2d.

| Arm | Stage-2 q ≤ 0.01 / 0.05 | CNN s | Outcome |
|---|---|---|---|
| grid2d reference | 717 / 1412 | 243 | frozen |
| A: topaz w32, seeds (about 100K params) | 0 / 853 | 228 | fail (overfits 600 seeds, heavy decoy tail) |
| B: A + MIL, 40000 bags × 40 epochs | not run | est. ≥ 1 h | stopped: too slow for OpenDIA |
| C: topaz w8, weight decay 0.01, seeds | 473 / 1177 | 146 | fail |
| D: C + MIL, 4000 bags × 10 epochs | 653 / **1539** | 253 | fail at 1%; best at 5%; decoy tail z 4 vs 7 |

**HeLa timsTOF diagnostic** (topaz w32, `-ss_iterations 1`; not a gate): **15119** targets at q ≤ 0.01, against 13253 for grid2d (+14%). CNN time was 342 s against about 145 s.

**Interpretation**
- The architecture is not slower than grid2d and does learn. On a strong IM run with 6000 seeds it clearly beats grid2d.
- On the weak Astral run, a few hundred seeds is too little supervision for it.
- MIL over all precursors fixes the overconfident decoy tail. Its lower R at the 0.05 cut (0.71 against a π0 estimate of 0.86) is a warning: with bag labels, part of the gain may be a target/decoy shift rather than peptide evidence.

**Decision:** grid2d plus stage 2 `search` stays frozen, and TOPAZ stays experimental. Raw 01 has now been used for four TOPAZ arms, so the next step should not be another raw 01 tweak. Options:
1. A second design run, or a fresh Astral replicate, for MIL with a symmetry guard at the 0.05 cut.
2. Feeding MIL-scale positives into seeds training (PyProphet-style iteration with a compact model).

**Port fidelity** ([`FIDELITY.md`](../gates/astral-raw01-cnn-topaz/FIDELITY.md))
- **Architecture: faithful.** A Rust-vs-C++ harness on identical inputs reproduces TOPAZ's co-elution head and MS1/MS2 agreement features to within 1e-6, with identical parameter counts per block.
- **Training regime: not faithful.** TOPAZ max-normalises each trace block before the dual view (`normalize_max: true`). It also uses cosine LR warm-up, gradient clipping and auxiliary pairwise and in-bag ranking losses. The OpenDIACNN arms had none of these.
- **Reading:** the failures on raw 01 tested this port's training regime, not TOPAZ as specified. The next TOPAZ variant should adopt these defaults as one declared change, on an unspent design run.
