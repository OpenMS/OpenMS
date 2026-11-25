# DIAG_FDR: why the GBM rescore makes own q loose on Astral (25 Nov 2026, diagnostic)

This is a diagnosis only. No plan was frozen, and nothing here changes a default.

- **Data:** local laptop runs of Astral raw 01 and raw 02, seeds 0-4 (10 runs, the dev set), 16 or 10 threads.
- **Library and inputs:** the arm-E library (`search.pqp`, r = 1.0) and its Parquet inputs.
- **Binary:** OpenMS `b7f19ff751`, sha256 in `binary.sha256`. `-features_out` adds the pre-GBM `network_score` and `network_q_value`. With no flag, raw01_s0 is byte-identical to the baseline (3/3 files).
- **Scripts and outputs:** in this directory. The main tables are `*_all10.txt`.
- **Held out:** Astral raw 03 and HeLa D/H were not used.

## Result 1: the GBM over-admits false targets only near q ≤ 0.01, not deeper

Pooled over the 10 dev runs (`diag_sweep_all10.txt`). D is the number of decoys above the threshold, and 2E is the entrapment estimate of false targets at r = 1.

| q cut | network N / D / 2E | GBM N / D / 2E |
|---|---|---|
| 0.005 | 6,569 / 21 / 40 | 16,882 / 84 / 100 |
| 0.01 | 18,674 / 179 / 136 | 20,121 / 221 / 236 |
| 0.02 | 23,295 / 458 / 446 | 23,888 / 535 / 486 |
| 0.05 | 27,353 / 1,379 / 1,460 | 27,793 / 1,529 / 1,362 |
| 0.1 | 31,001 / 3,178 / 3,166 | 31,216 / 3,524 / 3,054 |

- **Network:** 2E/D is about 1 at every depth.
- **GBM:** 2E/D is above 1 at q ≤ 0.01 and below 1 from q 0.02 on.
- **Reading:** the GBM does not create a global decoy/false-target mismatch. Its excess is confined to the top of the list, where only about 20 to 25 decoys per run set the q.

## Result 2: decoys and false targets are exchangeable in every bulk and tail test

All of these were measured on raw01_s0 and raw02_s0:

- **Separability:** a cross-fitted classifier separates decoys from entrapments with AUC 0.51 to 0.52, on all windows and on best windows (`diag_artifact.txt`). In the network tail (top 3 % of decoys) the AUC is 0.50 (`diag_shift.txt`).
- **GBM shift:** at matched network score, the GBM moves decoys and entrapments by the same amount (`diag_shift.txt`).
- **Ruled out as the cause** (`diag_windows.txt`, `diag_leak.txt`, `diag_relative.txt`, `diag_twin.txt`):
  - window counts per precursor;
  - cross-fold peak sharing;
  - cross-fold sequence relatives of the hard decoy negatives;
  - co-elution of the entrapment with its parent's kept window.

## Result 3: the excess sits in entrapments that share fragments with a human isobar, a class decoys cannot represent

- **Construction:** each entrapment is a shuffle of a human peptide that keeps K/R and the C-terminus, so it has the same precursor m/z as its parent.
- **Exposure:** share of each class with an isobaric human target that shares library fragments (`diag_exposure.txt`, 40,000 sampled per class).

| class | has isobaric human | shares ≥ 1 fragment | ≥ 2 | ≥ 3 |
|---|---:|---:|---:|---:|
| human targets | 0.096 | 0.080 | 0.070 | 0.059 |
| entrapment targets | 1.000 | 0.176 | 0.078 | 0.022 |
| decoys of human | 0.008 | 0.004 | 0.002 | 0.001 |
| decoys of entrapment | 0.008 | 0.004 | 0.003 | 0.001 |

- **Bound by class:** own q ≤ 0.01 split by the entrapment's fragment overlap with its parent (`summary_ov0_all10.txt`). Each class bound is E_c (1 + 1/r_c) / N, with r_c = r × the class share.

| score | N | E (ov 0 / 1 / 2+) | all-class bound | ov 0 bound | ov 2+ bound |
|---|---:|---|---:|---:|---:|
| network | 18,674 | 23 / 19 / 26 | 0.73 % | 0.27 % | 2.04 % |
| GBM | 20,121 | 49 / 24 / 45 | 1.17 % | 0.54 % | 3.28 % |

**Reading:**
- About 8 % of entrapments (and about 7 % of real human targets) have a near-isomer in the library that shares at least 2 of 6 fragments. Pseudo-reverse decoys almost never do (0.2 %). So no decoy stands in for a false target that borrows a real co-isolated peptide's signal.
- The network already over-admits this class (2.0 % on 7 % of entrapments). The GBM learns "shared fragment evidence present" from its positives and over-admits it further (3.3 %). It also doubles the ov-0 entrapments, from 23 to 49.
- On the remaining 83 % of entrapments, both scores are well inside 1 %.
- This is a limit of target-decoy competition with these decoys, not a GBM bug. The GBM amplifies it because it scores evidence that a decoy cannot have.
- A shadow rule tested offline was too weak to fix it (`sim_shadow_v1.txt`, −3 E of 118). It removes a precursor whose apex, m/z and ≥ 2 fragments coincide with a higher-scored target.

## Result 4: an offline per-precursor blend keeps the GBM's IDs at the network's calibration (dev data, simulation)

The blend is b = w × network + (1 − w) × GBM, using both calibrated scores of each precursor. It is then recalibrated with the tool's `calibrateTail` and recompeted with pooled paired TDC (`sim_combo.py`, `sim_combo_sweep_all10.txt`).

| arm | own q ≤ 0.01 N | E | pooled combined | worst run | runs > 1 % | median N @ FDP_c 1 % | sum N @ FDP_c 1 % |
|---|---:|---:|---:|---:|---:|---:|---:|
| network only | 18,674 | 68 | 0.73 % | 1.10 % | 3 | 2,101 | 20,877 |
| GBM (current default) | 20,121 | 118 | 1.17 % | 1.77 % | 6 | 2,136 | 20,474 |
| blend w = 0.5 | 20,291 | 76 | 0.75 % | 1.13 % | 3 | 2,165 | 21,422 |
| blend w = 0.35 | 20,566 | 82 | 0.80 % | 1.11 % | 4 | 2,154 | 21,330 |
| blend w = 0.65 | 19,758 | 75 | 0.76 % | 1.16 % | 3 | 2,140 | 21,269 |
| min(network, GBM) | 20,816 | 125 | 1.20 % | 1.70 % | 7 | 2,138 | 21,167 |

**Reading:**
- w = 0.5 keeps the GBM's own-q N (+0.8 %) and N at true FDP (median +1.4 %). The pooled bound is back at the network's level (0.75 % against 0.73 %).
- **Caveats:**
  - The weights were compared on the dev set itself, so w = 0.5 is the a-priori equal weight, not a tuned value.
  - The worst run still sits slightly above 1 % (raw01 seeds 0 and 1, 1.07 / 1.13 %), as it does for the network alone.
- **Why it works:** it is consistent with Result 3. A precursor must look good to both models, which were trained on different labels, so the excess of either one is diluted. That is an interpretation, not a proof.
- **Not yet known:** whether it holds on HeLa, where the GBM adds 7-9 % IDs. A blend may give some of that back.

## Next

1. **Done:** the C++ opt-in `-gbm_blend <w>` (`blendScores` in `OpenDIAGBM`, OpenMS `b3391c92e8`), default 0 = off.
   - **Tests:** 19/19 pass, including new `OpenDIAGBM_test` cases. These check that the kept rows are the GBM's, that w = 0 keeps the GBM order and w = 1 the network order, that the result is deterministic, and that w outside [0, 1] is rejected.
   - **Flag off:** raw01_s0 is byte-identical to the diagnosis run, 3/3 files.
   - **`-gbm_blend 0.5` on raw01_s0:** N 2047, E 11, combined 1.07 %, N @ FDP_c 1 % 2035. This equals the offline simulation (2047 / 11 / 2035). The tool's default on the same run is 2257 / 20 / 1.77 % / 2136. The record is `run_blend_check.sh` and `blend_binary.sha256`.
2. A frozen precheck before any remote run:
   - **arms:** w = 0.5 only, fixed now from this note;
   - **dev set:** Astral raw 01 and raw 02, seeds 0-4, on codon;
   - **held out:** raw 03 seeds 0-4, and HeLa D and H at r ≈ 7, seeds 0-4;
   - **criteria:** set in the precheck, covering the own-q bound on both instruments and N @ FDP_c 1 % at or above the current default per set.
