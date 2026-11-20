# PLAN_FAINT: faint trace evidence features on the r ≈ 7 entrapment (frozen 20 Nov 2026)

**Question.** Do two cheap, pooled trace-shape features for faint peaks recover IDs at equal true-FDP without materially loosening the tool's q values?

## Diagnostic evidence before this freeze

The fresh current-default r ≈ 7 D decomposition found 8,787 DIA-NN IDs whose right peak was present in OpenDIA's best window but scored too low. On this B1 class, the OpenDIA MS2-area rank correlation with DIA-NN quantity was 0.543 versus 0.776 for already accepted A IDs, and fragment coelution fell toward the decoy distribution. RT and calibrated ion mobility were already aligned, so the proposed change targets trace evidence rather than candidate windows.

The C++ implementation adds two inputs to the existing 36-feature scorer, without another extraction pass:

- `p_cos_apex`: cosine agreement of the pooled fragment trace around the candidate apex.
- `p_local_sn`: candidate-apex signal divided by the median outside the picked peak group.

The block is opt-in with `-faint_features`. The default path remains unchanged. The code is OpenMS commit `11117425cc` on `feature/opendia-clean-rewrite`; the affected build passes 25/25 OpenDIA and calibration tests.

The r ≈ 1.56 two-seed precheck used the same current and faint binaries, library, raw file, and seeds 0 and 1:

| Arm | Median N at FDPc 1 % | Pooled lower | Pooled combined | Pair win |
|---|---:|---:|---:|---:|
| CUR | 14,908 | 0.82 % | 1.56 % | 0.523 |
| FAINT | 15,849 | 0.85 % | 1.62 % | 0.512 |

The faint arm therefore gained 6.3 % in median N at FDPc 1 %, increased the lower bound by 0.03 percentage points and the combined bound by 0.06 percentage points, and moved pair fairness toward 0.5. The arm's wall time was 3:46 versus 2:25 to 2:59 for CUR in these runs, with peak RSS about 9.8 to 10.0 GB. The choice to proceed was made after seeing this precheck and is disclosed here.

## Setup (fixed)

- **Binary and container:** source commit `11117425cc`. The tested container is `/nfs/research/juan/singj/opendia/opendia-pd10-2257c55c7a.sif`, sha256 `3739e33e9952b9245184a8e8b91d554fb47f990bcf06044d13fe21c6fe64e6f2`.
- **Library:** the existing r ≈ 7 entrapment library used by the current default, with its decoy construction and cache unchanged.
- **Arms:** CUR is the current default. NEW is CUR plus `-faint_features`.
- **Runs:** D seeds 0 through 4, 8 threads, with the same raw file, library, output layout, decoy mode, and calibration settings as the current-default r ≈ 7 benchmark. H is submitted only after the D verdict is committed and is reported either way.
- **Evaluation:** the frozen `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`, on each tool output as written. Entrapment labels are read only by that evaluator.

## Criteria on D

1. **Sensitivity:** NEW median N at FDPc 1 % is at least 1.05 times CUR median.
2. **Calibration is not materially worse:** NEW pooled lower bound at own q ≤ 0.01 is no more than CUR plus 0.0020, and NEW pooled combined bound is no more than CUR plus 0.0030.
3. **Pair fairness:** NEW entrapment versus own-decoy win probability is no more than 0.56.
4. **Cost:** NEW median wall time is no more than 1.3 times CUR, and peak RSS is no more than 64 GB.

All criteria and thresholds are frozen before the r ≈ 7 D runs. No scorer retuning, feature removal, seed removal, or post-hoc filtering is allowed.

## Outcome rule

- If all four criteria pass, report H and then decide separately whether to make the feature block the default.
- If sensitivity passes but criterion 2 fails, report the gain together with the calibration cost and keep the block opt-in. Do not retune.
- Any other failure is reported as measured, with the block remaining opt-in.
