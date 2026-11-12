# PLAN_XIC_IDS: identifications after the extraction-parity changes (D11FWT recipe unchanged)

Written 2026-11-12, before any search, population, model or score of this rebuild exists. Root `xicids-20261112/`.

## Question
The extraction now matches DIA-NN much more closely on raw01. The changes are: a monotone RT map with run-edge
windows, a fitted fragment tolerance, and the xic picker with profile-gated co-elution and the estimated boundary cap
(commits 92cc9a2469, 8d2771bcec, 24dad4f842). Do these changes give more identifications, and do they keep
entrapments and decoys low? Recipe, seeds and thresholds stay as they are.

## What changes against the rtlin rebuild (`rtlin-20261111`)
1. Binaries: `OpenMS-r1-full-ninja-build/bin` at 24dad4f842 (`binary.expected.sha256`; the chain stops on mismatch).
2. Every OpenDIA search (the pipeline search and the 6- and 12-fragment exports) uses
   `-min_transitions 1 -peak_picking xic -xic_max_peaks 5` in place of
   `-min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10`. The boundary cap comes from
   cal.json, and `-xic_min_profile_scans` is at its default of 2.
   - A candidate window is now one of up to 5 xic peaks per precursor, inside the calibrated RT window.
   - The old candidate was a streaming signal run, split at local apexes, which ended only after a 30 s gap.
   - The 30 s gap rule therefore no longer forms candidates.
   - The search features keep their names. In xic mode they mean:
     - `transition_support`: fragments seen within 1 scan of the apex.
     - `matched_scans`: signal scans in the capped peak.
     - `peak_width_seconds`: the capped width.
3. Everything else is unchanged: calibration options, `prepare_trace_windows.py`, the BASE OpenDIACNN crossfit
   (`-stage2 search+boundary+byion`, same seeds, folds and epochs), the population rule (pairs with max
   transition_support ≥ 3), the 12-fragment library, the arrays, the features, and `train_v4.py D11F` →
   `make_w.py` → `calib_tail.py`. The scripts are copies of `rtlin-20261111/scripts` with only the root path changed.
   One script is new, the evaluation-only `report_xic.py`.

Pre-look (no labels, apart from one DIA-NN-membership count, done before this plan; raw01 picker2 search of the same
binaries):

| | Pairs | Precursors | Windows | DIA-NN q ≤ 1 % in population |
|---|---|---|---|---|
| xic | 144,601 | 288,722 | 1.43 M | 2345 / 2357 |
| rtlin | 85,422 | 169,169 | 1.16 M | 2301 / 2357 |

## Chain
`chain.sh`: calibrate, search, seeds, BASE and population for each raw, then pqp12, then per raw the exports, voxels,
build12 and features12, then training and evaluation.
- Each step takes under about 25 min. The total is about 3 h, run locally because the cluster is not reachable from
  this session.
- A memory guard waits before BASE until there is enough free memory. It does not change results.
- One run. No retuning after a failure.
- Precheck before freezing: the xic search was rerun on raw01 with the picker2 calibration, to confirm it is
  deterministic. The chain's export step needs a byte-identical `search.features.tsv`. The precheck output is in
  `precheck/`.

## Pass rule (unchanged: PLAN_V2 rules 1–4, PLAN_V8/V19 additions), on both raws
1. Against the rebuilt BASE:
   - The bootstrap 90 % interval of N@2 % − BASE lies above 0.
   - The interval of N@1 % − BASE is not entirely below 0.
   - FDP_c at q ≤ 0.01 is ≤ 0.015.
   - The entrapment-over-decoy winner AUC is ≤ the rebuilt BASE AUC.
2. N@FDP_c 2 % is above the current base D11FWT: 1716 (raw01) and 1736 (raw02).
3. The null-window winner AUC is ≤ D11FWT + 0.005: 0.4983 (raw01) and 0.4986 (raw02). The references are the base
   D11FWT values 0.4933 / 0.4936 from NULL_V19.txt.

The N@2 % reference 1716 / 1736 and the null references are those of the current base D11FWT, which is the old
population (`astral-voxel-pilot-20261110/v2`). The rtlin D11FWT failed PLAN_V8 and is not the base.

## Reported (descriptive, labels only in evaluation)
- N at FDP_c 0.5 / 1 / 1.5 / 2 %.
- q ≤ 0.01 / 0.05 N and FDP_c.
- `report_xic.py` against the old (`astral-voxel-pilot-20261110/v2`) and the rtlin D11FWT:
  - q ≤ 0.01 human / entrapment / paired decoys;
  - the N@FDP_c 1 % sets;
  - gained and lost precursors;
  - raw01 overlap with the 2357 DIA-NN q ≤ 1 % IDs.
- The A–E miss funnel (`misses.py`).
- Null AUC for BASE and D11FWT.
- Wall time and peak memory of each step, next to DIA-NN's 11:03 on raw01.

Absolute increases in entrapments or decoys are reported as costs even when the rule passes.

## Known limits
- The populations differ between builds, so comparisons with old or rtlin are not paired. Precursor sets are compared.
- D11F was designed on run-level search features. In xic mode three of its inputs change meaning (see above). This
  measures the recipe as it is. A failure is not evidence against the extraction.
- The boundary cap k and the tolerance are estimated per raw by the calibration. Nothing is set from raw01.
- Raw02 is replication, not a pristine holdout.

## Outcome
- Pass: the xic build becomes the base for later work.
- Fail: the result is recorded, nothing is reverted, and the next step is the user's decision.

## Amendment 1 (2026-11-12T20:02Z, before any BASE output; no result had been seen)
- Reason: raw01 BASE needed about 27 GB (8.9 KiB per window x 2.75 M windows) and the machine had 22-25 GB free.
  The user said to go ahead and run it.
- Change: OpenDIACNN binary replaced by commit 40af886a4a. Without library ion mobility the trace tensor keeps only the
  central mobility bin, where every point already landed; the network still sees the same [6,32,8] input. Results are
  meant to be unchanged.
- Evidence: raw01 70k-window crossfit subset (all seeds plus every 40th window, 4 epochs), old binary 0298804e vs new
  70cc963c: byte-identical crossfit.tsv, features.tsv and the 3 fold models (`~/.jcode/scratch/compactim/eq`).
  Peak RSS 5.06 -> 4.69 GB.
- Chain edits (no effect on results): the raw01 cal, search and seeds already produced under this plan are reused;
  the memory guard is 3.9 KiB x windows + 1 GB. OpenDIA and OpenDIACalibrate are unchanged.
- Nothing else in the plan, the scripts or the pass rule changes.
