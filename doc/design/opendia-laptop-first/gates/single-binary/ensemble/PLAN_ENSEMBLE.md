# PLAN_ENSEMBLE: seed-ensemble scoring in OpenDIAWorkflow (option 1, scorer capacity)

Written 16 November 2026, before the implementation and before any run.

## Evidence, offline on the arm-E arrays with the C++ scorer (`gap/OPTION1_PROBES.txt`)

| probe | raw01 N@1 % | raw02 N@1 % | verdict |
|---|---|---|---|
| single seed, seeds 0-4 (median) | 1822 | 1905 | baseline |
| 5 signal-level interaction inputs, seeds 0-2 | 1658-1759 | 1658-1695 | worse; FDP_c up to 0.0116; dropped |
| self-trained gradient-boosted trees on the network's winners | unstable: 616-1754 | 592-1827 | dropped |
| mean of 2 seeds' decoy-z scores (5 disjoint pairs, median) | 1850 | 1910 | +28 / +5 |
| mean of 3 seeds (5 sets, median) | 1848 (1826-1855) | 1939 (1889-1952) | +26 / +34; spread 29 / 63 against 60 / 195 single |

- At q ≤ 0.01 the 3-seed ensembles have FDP_c 0.0059-0.0071 (raw01) and 0.0025-0.0049 (raw02).
- The labelled ceiling on the same inputs (`CEILING_current.txt`) is about +100 per raw. The ensemble reaches a quarter
  to a third of it.
- Note: both probes use offline arrays. The interaction inputs were constructed after seeing PRECHECK_JOINT.

## Design (`-ensemble <m>` in OpenDIAWorkflow; plan value m = 3)
1. **Train.** The CrossFitScorer is trained m times, with seeds s, s+1000, …, s+1000(m−1), on the same input. Each
   training is unchanged.
2. **Combine.** Each model's raw window scores are standardised per fold by the median and 1.4826 × MAD of the decoy
   rows of that fold. A row's ensemble score is the mean over the models.
3. **Then as now.** The W rule keeps each precursor's window with the largest mean linear part (standardised the same
   way). Then come tail calibration and pooled paired TDC.

The combination needs no new constants: it uses the scorer's own decoy units. It costs m × the scoring time, which is
about 20 s per model on 8 threads. Memory stays the same, because the models run one after another.

## Arms, runs, criteria
- **Arms:** CUR = m = 1 (today), NEW = m = 3. Both are separate OpenDIAWorkflow runs: self-calibrating, arm-E library
  plus `library12_full.pqp`, seeds 0-4.
- **Runs:** raws 01/02 are the design runs. Raw 03 is held out and evaluated once, after the verdict is committed.
- **Criteria**, on raws 01 and 02 separately:
  1. NEW median N@FDP_c 1 % ≥ CUR median + 1 %.
  2. NEW FDP_c at q ≤ 0.01 ≤ 0.010 on every seed.
  3. NEW seed range (max − min of N@1 %) ≤ CUR seed range.
- **Outcome:** pass makes m = 3 the default. Fail removes the code. Nothing is retuned.
