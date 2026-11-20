# PRECHECK_FAINT: do per-scan trace features add information for faint peaks? (20 Nov 2026, diagnostic)

This is a diagnostic, not a plan. DIA-NN labels are used for evaluation only. It follows PRECHECK_NOCAND: with the
looser search (`-valley_fraction 0.8 -max_peaks 8`) 98 % of DIA-NN IDs have a right candidate, but the faint ones score
below the decoy tail. The question here is whether evidence the 36 features do not use can separate them.

**Data.**
- r ≈ 1.56 HeLa D, seed 0, the valley 0.8 / max 8 arm (`ev/v08m8`).
- A mini library of 6,000 DIA-NN non-entrapment q ≤ 0.01 targets (about 1,500 per quantity quartile) and 30,000 random decoy precursors.
- Searched with the dev tool and the arm's calibration and picker settings. The 266,207 candidate windows join 1:1 to the workflow's 36 features. 5,909 of the 6,000 targets have a right window, meaning an apex within 6 s of DIA-NN's RT.
- Scripts: `faint_mini.py`, `precheck_faint.py`, `precheck_faint_boot.py`.

**New features.** Each candidate peak gets 7 new features, computed from the search's own trace matrix (6 fragments × native scans over the peak):

| Feature | What it measures |
|---|---|
| `ps_cos_w` | Mean cosine of each scan's fragment vector to the library pattern (square-root intensities), weighted by consensus |
| `ps_cos_apex` | The same cosine on the summed apex ±1 scans |
| `ps_n4` | Share of peak scans with at least 4 fragments present |
| `apex_spread` | Spread of the per-fragment apex positions |
| `frag_run` | Longest run of consecutive nonzero scans, relative to the peak width |
| `local_sn` | Apex height against the median of the trace outside the peak |
| `rank_corr` | Spearman correlation of the apex fragment sums with the library intensities |

**Test.** A GBM trained on DIA-NN labels and cross-fitted on the workflow's folds, with the 36 features vs 36 + new. This
is a ceiling comparison: a gain means the new features carry new information, not that a method is ready. The measure
is the TPR of each target's right window against the decoy precursors' best windows, at a fixed decoy FPR.

| Scorer | TPR @ FPR 1e-3 (Q1 / Q2 / Q3 / Q4 / all) | TPR @ FPR 1e-2 (all) |
|---|---|---|
| Tool score | 0.215 / 0.414 / 0.600 / 0.871 / 0.519 | 0.590 |
| Ceiling, 36 features | 0.197 / 0.407 / 0.592 / 0.865 / 0.509 | 0.658 |
| Ceiling, 36 + 7 new | 0.218 / 0.432 / 0.633 / 0.879 / **0.535** | **0.693** |
| 36 + `ps_cos_apex` only | 0.527 | 0.673 |
| 36 + `local_sn` only | 0.502 | 0.673 |
| 36 + any other single new feature | 0.500-0.516 | 0.657-0.665 |

Paired bootstrap (300 resamples of decoys and targets) on the TPR change from 36 to 43 features, for two GBM seeds:

| Decoy FPR | Seed 0 | Seed 1 |
|---|---|---|
| 1e-3 | +0.026 (90 % +0.008 to +0.064) | +0.037 (+0.018 to +0.053) |
| 3e-3 | +0.029 (+0.014 to +0.042) | +0.029 (+0.019 to +0.041) |
| 1e-2 | +0.036 (+0.028 to +0.044) | +0.032 (+0.024 to +0.040) |

Univariate AUC, Q1-Q2 right windows against decoy best windows:

| Feature | AUC |
|---|---|
| ps_cos_apex | 0.784 |
| ps_cos_w | 0.758 |
| f12_libcos (existing) | 0.731 |
| frag_run | 0.707 |
| boundary_coelution (existing) | 0.676 |
| ps_n4 | 0.673 |
| local_sn | 0.662 |
| rank_corr | 0.647 |
| x_library_correlation (existing) | 0.608 |
| apex_spread | 0.504 |

## Reading

- **The information is real but small.** The ceiling gains about 3 TPR points at the FPRs where the 1 % cut falls, which is +5-7 % relative. This is consistent across seeds and FPRs, and the bootstrap intervals exclude zero.
- Most of the gain comes from **library-pattern agreement at the apex** (`ps_cos_apex`) and **local signal-to-background** (`local_sn`). Both are cheap to compute from the 6-fragment tensor the workflow already extracts.
- The scan-by-scan shape features add little: fragment-count share, apex spread and run length.
- **The gain is concentrated in Q1-Q3** but is small in absolute terms: Q1 goes from 0.20 to 0.22 at FPR 1e-3. Even with these features, about 78 % of the faintest DIA-NN IDs remain below the decoy tail.
- **Caveats.** This is a DIA-NN-labelled ceiling, so a label-free scorer may realise less of the gain. It is one run and one seed of the search, with 30,000 decoys, so FPR 1e-3 rests on 30 decoys. The bootstrap includes that noise.

## Decision

- Worth one cheap implementation: `ps_cos_apex` and `local_sn` (plus `ps_cos_w`) as additional stage-2 features computed from the existing trace tensor. Expect at most about +5 % IDs, and less without labels.
- Not a route to DIA-NN parity. A gap of this size (OpenDIA at 44-62 % of DIA-NN) needs much stronger per-peak evidence than these pattern statistics.
- Measuring the remaining distance would need a different kind of test: DIA-NN's own per-precursor scores on these windows, or a learned trace model (the earlier OpenDIACNN line).
- No defaults changed.
