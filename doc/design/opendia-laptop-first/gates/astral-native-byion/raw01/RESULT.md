# Stage-2 byion group: raw 01 native design result

`PLAN.md` was written before the run (`PLAN.written_at`, sha256 in `PLAN.md.sha256`). Binary is commit `8c49a5bbb2`.

## Equivalence
- Columns 1–6 and 9 of `crossfit.tsv` (sample, precursor, group, decoy, fold, CNN logit, window) are byte-identical to `../stage2-boundary-design/boundary`.
- Feature columns 1–14 (ids, cnn_z, search, boundary) are identical.
- Only `byion_b`/`byion_y`, the stage-2 score, `precursor_best` and q differ.
- `run.sh` compared columns 1–7 and 1–15. Those include `precursor_best` and the first byion column, so its `differ` lines are expected. See `CMP.txt` for the corrected comparison.

## Result (entrapment-matched N, the `matched_fdp.py` metric; FDP_c combined)

| | FDP_c ≤ 0.5% | **1%** | 1.5% | 2% | q≤0.01 N / FDP_c | q≤0.05 N / FDP_c |
|---|---|---|---|---|---|---|
| search+boundary | 641 | 757 | 923 | 950 | 746 / 0.008 | 1441 / 0.070 |
| **+byion** | 899 | **1029 (+35.9%)** | 1270 | 1321 | 1068 / 0.013 | 1638 / 0.065 |
| OSW+Percolator (same metric) | 1015 | 1151 | 1251 | 1286 | | |

- **PASS on N:** 1029 ≥ 795.
- **Plain pooled q fails the 0.06 bar:** q ≤ 0.05 has FDP_c 0.065.
- **Bulk-null weighted q passes (`bulknull_native.py`, 3 seeds):**
  - q ≤ 0.01: 939 / FDP_c 0.006;
  - q ≤ 0.05: 1520–1537 / FDP_c 0.048–0.049.
- So the predeclared rule holds with bulk-null as the primary q, as recommended in CALIB section 9.
- **Fold weights:** `byion_b` (+2.8 to +2.9) and `byion_y` (+3.6 to +3.9) carry the second- and third-largest standardised weights, after the search terms. They are stable across folds.
- **Cost:** wall time 12:58 against 12:52, and peak RSS 24.9 GB against 24.7 GB.

## Reading
- The native stage-2 LDA uses the counts at least as well as the offline MLP (+36% here against +15% there). The native baseline is weaker, so there was more to gain.
- At 0.5–1.5% we are now at 89% of OSW+Percolator at FDP_c 1%, and above it at 1.5% and 2%.
- Raw 02 replication: `../../astral-raw02-byion-design-20261109`.
