# RESULT_XIC_IDS (PLAN_XIC_IDS.md, amendment 1)

Run 2026-11-12, 20:02-23:03 UTC, binaries in `binary.expected.sha256` (OpenDIACNN 40af886a4a, results byte-identical to
the frozen 0298804e on a 70k-window subset). Chain status `done`; the funnel step first failed on a missing output
directory (`v3seq/`), it was created and only the two descriptive scripts (misses.py, report_xic.py) were rerun.

## Verdict: FAIL (rule 1, FDP_c at q <= 0.01)

| Rule | raw01 | raw02 | |
|---|---|---|---|
| 1a N@2 % - BASE, 90 % interval > 0 | +432..+789 | +610..+1044 | pass |
| 1b N@1 % - BASE not entirely < 0 | +390..+743 | +546..+1051 | pass |
| 1c FDP_c at q <= 0.01 <= 0.015 | 0.016 | 0.017 | **fail** |
| 1d entrapment>decoy AUC <= BASE | 0.4907 <= 0.5190 | 0.4999 <= 0.5152 | pass |
| 2 N@FDP_c 2 % > 1716 / 1736 | 1862 | 1751 | pass |
| 3 null winner AUC <= 0.4983 / 0.4986 | 0.4966 | 0.4981 | pass |

Per the plan: recorded, nothing reverted, next step is the user's decision.

## IDs, D11FWT

| | old base | rtlin | xic (this) |
|---|---|---|---|
| raw01 N@FDP_c 0.5 / 1 / 1.5 / 2 % | 971 / 1496 / 1616 / 1716 | 1142 / 1519 / 1694 / 1839 | 1414 / 1567 / 1716 / 1862 |
| raw02 N@FDP_c 0.5 / 1 / 1.5 / 2 % | 1127 / 1426 / 1625 / 1736 | 1250 / 1420 / 1644 / 1736 | 1302 / 1558 / 1651 / 1751 |
| raw01 q <= 0.01 (FDP_c) | 1614 (0.014) | 1728 (0.016) | 1767 (0.016) |
| raw02 q <= 0.01 (FDP_c) | 1536 (0.013) | 1716 (0.019) | 1653 (0.017) |
| raw01 q <= 0.01 entrap / paired decoys | 11 / 15 | 14 / 16 | 14 / 16 |
| raw02 q <= 0.01 entrap / paired decoys | 10 / 14 | 16 / 16 | 14 / 15 |
| entrapments in N@FDP_c 1 % | 7 / 7 | 7 / 7 | 7 / 7 |

- At a fixed entrapment FDP (the ID count that matters), xic gains on both raws: +71 / +132 at 1 %, +443 / +175 at
  0.5 %, +146 / +15 at 2 % against the old base. raw02 at 2 % is only just above the reference.
- The q-value is miscalibrated: at q <= 0.01 the entrapment FDP is 1.6-1.7 %, so rule 1c fails. The old base was at
  1.3-1.4 %. The extra IDs at q <= 0.01 cost 3-4 entrapments and 1 paired decoy per raw.
- raw01 against DIA-NN's 2357 q <= 1 % IDs: in the population 2250 -> 2345, accepted at q <= 0.01 1490 -> 1616
  (entrapments 3 -> 3); accepted but not in DIA-NN 124 -> 151 (entrapments 8 -> 11).
- Population 1.27 M -> 1.43 M windows, 187 k -> 289 k precursors (no 30 s gap filter, run-edge windows).

## Where DIA-NN IDs are lost (raw01 funnel, human exact-match reference 2743)
A not in population 42, B no window within 10 s 55, C wrong window won 82, D right window below cut 1074,
E accepted 1490. Extraction stages A-C are 179 (6.5 %); the bulk (1074) has the right window and loses on score.
raw02: 36 / 71 / 83 / 1055 / 1486.

## Cost (raw01 / raw02)
OpenDIA search 3:04 / 2:50, 4.4 GB. BASE (OpenDIACNN crossfit) 14:53 / 12:33, 9.95 GB (was 23.1 GB before the
compact tensor). export6 20:56 / 20:43, 8.7 GB. Whole chain per raw about 1.5 h; DIA-NN raw01 11:03.

## Cleanup after the run (OpenMS a8be5bed3c, not part of the frozen plan)
The run-mode extraction is removed; native-scan XIC picking is the only mode. Option renames for scripts written
against the plan binary: `-peak_picking xic` is gone, `-xic_max_peaks` -> `-max_peaks`, `-xic_valley_fraction` ->
`-valley_fraction`, `-xic_boundary_half_widths` -> `-boundary_half_widths`, `-xic_min_profile_scans` ->
`-min_profile_scans`, `-xic_traces_out/_precursors` -> `-traces_out/_precursors`. The feature table drops the two
always-zero columns `has_im_match_evidence`, `mean_abs_im_error_vssc` (19 columns). raw01 check: peak table
byte-identical to `raw01/pipeline/search.tsv`, feature table identical after cutting those two columns, identical at
1 and 10 threads. Search wall time on raw01: 2:07 -> 0:51 at `-threads 10` (1:10 at 1 thread); 32 s of that is
streaming the RAW file. Peak RSS 4.4 GB (unchanged).
