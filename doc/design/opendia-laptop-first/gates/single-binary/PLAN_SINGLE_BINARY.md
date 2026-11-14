# PLAN_SINGLE_BINARY: one C++ OpenDIA run, raw files to report

Written 14 November 2026. Goal: one `OpenDIA` binary that analyses all raw files end to end like DIA-NN, with the
identification performance of the best validated chain (PLAN_GAP D11FDI with the C++ interference features, arm D11FDC).

## Where we are

The best result needs 4 C++ tools and about 12 Python scripts, and reads each raw file 6 or more times:

| Step | Now | Target |
|---|---|---|
| Calibration | OpenDIACalibrate (1-3 raw passes) | OpenDIA, same code |
| Candidate search, 11 heuristic features | OpenDIA (pass), OpenDIACNN (boundary, byion; LibTorch, CNN unused by D11FDI) | one raw pass after calibration |
| Seed windows, population, training rows | prepare_trace_windows.py, make_population.py | in memory |
| 12-fragment features, deficit | OpenDIA 2 more passes, build_pqp12, build12, features12, deficit_features (Python) | same raw pass |
| Interference | OpenDIAInterference (C++, own raw pass) | same raw pass |
| Scorer | train_gap.py (torch MLP + decoy-whitened LDA readout), make_w, calib_tail | C++ library, no LibTorch |
| q-values, report | evaluate (Python, labels) | pooled paired TDC, precursor report |

The fully C++ scorer today (CNN + LDA stage2_score) gives N@FDP_c 1 % 1170 / 267 (raw01 / raw02), against 1837-1857 /
1903-1976 for D11FDC over five seeds.

## Steps and acceptance

1. **C++ scorer** (`OpenDIA::CrossFitScorer`): the D11 network h = [z, tanh(g(z))], g = d-64-16, the decoy-whitened
   readout w = (S_- + 0.01 I)^-1 (mu_+ - mu_-), loss -dmu'w, Adam 1e-3, 30 epochs x 100 steps of 256 + 1024, early stop
   on the validation objective (patience 5), 3 folds, rounds 0 (anchor rows) and 1, 2 (self-labels at pooled q <= 0.15,
   0.05), the W rule (window by the linear part), per-fold decoy-tail calibration, pooled paired TDC. Same validation
   split (FNV-1a of 'v' + family / pair key). Exact reproduction of torch is not possible (RNG, float32), so acceptance
   is statistical, on the arm-E arrays (heuristics, feat12, featD, featC), over five C++ seeds:
   - median N@FDP_c 1 % at least 98 % of the D11FDC median (raw01 1848, raw02 1925), each raw;
   - FDP_c at q <= 0.01 <= 0.015 for every seed;
   - deterministic for a seed and independent of the thread count.
2. **Feature extraction in OpenDIA**: the 11 heuristics (search 4, boundary 5, byion 2), feat12 (6), featD (3) and featC
   (4) from one raw pass after calibration, without the xipm export. Acceptance per block: equal to the current arrays on
   the same inputs (tolerance per block, documented), then step 1's test on the C++ arrays.
3. **End-to-end mode**: `OpenDIA -in <raws> -lib <pqp> -out_dir <dir>` runs calibration, search, features, population
   and training rows, scorer and TDC per raw, and writes a precursor report (run, precursor, peptide, charge, RT, q,
   score, quantity).
4. **Several raws**: per-run q and a global precursor q over runs (best score per precursor over runs, pooled TDC).
   Match-between-runs later.
5. **Acceptance on raw01 / raw02**: N@FDP_c 1 % and FDP_c at q <= 0.01 as step 1 against the Python chain; wall time and
   peak memory per raw reported.

The Python scripts stay for evaluation against entrapment labels only. OpenDIACNN stays optional (OPENDIA_WITH_TORCH)
and leaves the default path.
