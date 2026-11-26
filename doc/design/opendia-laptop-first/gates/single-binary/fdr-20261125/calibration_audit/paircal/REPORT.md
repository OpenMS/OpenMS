# Fixed pair-stratified TDC screen results

2026-11-26. **Both fixed rules fail the success criterion. Do not promote either as a calibration fix.**

## Protocol

Exact coordinator-specified rules are frozen in RULES.md before this screen. A separates pairs whose NET and G0 winner identities agree/disagree. B separates pairs whose max NET score is above/below the 99th percentile of all NET-kept decoy scores. A tiny missing-NET class is kept explicitly in both rules rather than dropping G0 pairs. No biological annotation enters pairing, strata, scores, thresholds, or selection. No tuning after results.

Data: all D and H cur/gbt seed pairs 0..4, ten complete runs, about 1.2 million competing pairs per run. Same G0 pair-winner universe for baseline/A/B within each run. Plus-one TDC recomputed identically within fixed strata, equal-score groups together, exact member ties to decoy. This recomputed baseline is the comparator, not a promise to reproduce the tool's original fold-standardized q-values.

## Corrected-label evaluation

Combined bound = E*(1+1/r)/N at estimated q<=1%. Exact-label correction is evaluation only. Intervals are 90% Poisson count intervals, not a validation guarantee.

| Data | Rule | Pooled N | Pooled E | Combined bound | 90% interval | Median N@FDP_c 1% | Delta vs baseline |
|---|---|---:|---:|---:|---:|---:|---:|
| D | baseline | 109910 | 1032 | 1.0747% | 1.0203–1.1314% | 21805 | 0 |
| D | A agreement | 110027 | 1045 | 1.0871% | 1.0323–1.1440% | 21767 | -38 |
| D | B net99 | 109742 | 1024 | 1.0680% | 1.0137–1.1245% | 21779 | -26 |
| H | baseline | 192239 | 2144 | 1.2765% | 1.2315–1.3228% | 37897 | 0 |
| H | A agreement | 191083 | 2106 | 1.2615% | 1.2166–1.3076% | 37699 | -198 |
| H | B net99 | 190646 | 2031 | 1.2193% | 1.1752–1.2648% | 37747 | -150 |

Corrected per-seed N@FDP_c 1%, seeds 0 through 4:

- D baseline: 21990, 21541, 22009, 21776, 21805.
- D A: 21935, 21546, 22023, 21755, 21767.
- D B: 21905, 21470, 21997, 21718, 21779.
- H baseline: 37928, 37876, 37897, 37959, 37719.
- H A: 37660, 37816, 37704, 37699, 37543.
- H B: 37659, 37780, 37754, 37747, 37573.

Both candidates remain above the 1% corrected pooled bound and both lose median N@FDP_c1% in D and H. B modestly improves the H calibration diagnostic but does not satisfy the fixed objective. A slightly worsens the D bound. These are exploratory reused-data results, not held-out acceptance.

## Original-label evaluation

| Data | Rule | Combined bound | Median N@FDP_c 1% |
|---|---|---:|---:|
| D | baseline | 1.2487% | 21495 |
| D | A | 1.2609% | 21470 |
| D | B | 1.2381% | 21452 |
| H | baseline | 1.4558% | 37376 |
| H | A | 1.4419% | 37188 |
| H | B | 1.3991% | 37246 |

## Artifacts and observed checks

- `RULES.md`: frozen formulas and handling decisions, before evaluation.
- `paircal_helpers.py`: annotation-free SQL and TDC unit checks. All helper checks passed, including all-target +1 correction, all-decoy q=1, tie-group handling, and row-permutation invariance.
- `run_selection.py`: actual fixed selector. Only numeric IDs, pair keys, target/decoy flags and scores enter selection. Existing cached library columns other than precursor_id/key/decoy are discarded before processing.
- `selected.tsv.gz`: every pair winner, stratum, score and q-value from all 30 run-rule combinations, compressed after evaluation to conserve disk.
- `strata.tsv`: per-stratum pair counts, target/decoy winner counts, frozen NET-decoy 99th-percentile cutoff.
- `selection.log`: all seeds and all three rule outputs, exit 0 after 266 s.
- `evaluate.py`: separate evaluation-only process that reads protein annotations and the preexisting exact-label map. These annotations cannot enter the already completed selector process.
- `evaluation.tsv`: all 60 policy/data/seed/rule rows. `evaluation.txt`: full table plus pooled bounds and intervals. Evaluation exit 0 after 55 s.
- `provenance.sha256`: hashes of frozen rules, selector, helper and exact-label map.
- `run.sh`: rerun selection and isolated evaluation. Uses preexisting local Python environment and annotation-free metadata cache.

No repository edits, source commits, refs changes, remote jobs or pushes. No GBM model was trained here. Original source symmetry diagnostic remains under sibling symmetry-20261126 directory.
