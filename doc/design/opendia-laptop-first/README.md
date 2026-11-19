# OpenDIA rewrite: master plan

**Status, 18 November 2026.** This plan replaced the laptop-first WP-00 to WP-10 plan on 7 November 2026. The old plan is closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
All of its code is under `src/topp/OpenDIA/`, apart from the Thermo read-ahead change, which is split out as OpenMS#10925 (approved). The branch tip is `963543c7f0`.
It is an experimental branch, not a pull request. It must not merge into `develop` until the acceptance gates of R-WP6 pass.
The branch's [`src/topp/OpenDIA/README.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README.md) describes the tools. The evidence is in [`gates/`](gates/).

## Current state, 18 November 2026

**Summary.**
- **Astral (entrapment, true FDR):** OpenDIA reaches 81-84 % of DIA-NN at FDP_c 1 %, unchanged since 16 Nov (the table
  below). Every lever tried there has failed or shown no room.
- **HeLa diaPASEF (own q ≤ 0.01, no entrapment):** three fixes on 17 Nov raised the design run from 11,653 to 16,116 and
  the held-out run from 19,725 to 24,000. That is 56 % / 63 % of DIA-NN and 84 % of OpenSWATH SVM on the design run.
  Each fix passed a frozen plan with a held-out run and is byte-identical on Astral. See "Current HeLa result" below,
  RESULT_RTCENTRE.md and RESULT_IMGRID.md.
- **Peak-group scores (PLAN_PGSCORES, 18 Nov): failed, removed.** Ten OpenSWATH-style scores raised HeLa by about
  17 % (D 18,792, held-out H 28,020) and Astral N@FDP_c 1 % by 2-12 % (92 / 83 / 90 % of DIA-NN). They failed two
  frozen criteria: the DIA-NN agreement on D (0.9104 against a limit of 0.9249), and q calibration on Astral raw 02
  (FDP_c 0.0103 median at q ≤ 0.01). The evidence is real, but some of it separates decoys better than entrapments.
  A follow-up check (PRECHECK_PGQ) found no single score to blame: on every pg score, entrapments and decoys look the same
  (AUC 0.45-0.55). The excess is on raw 02 only. Pooled over raws 01-03 it is 0.0072, so the scorer's fit is the next
  suspect. See RESULT_PGSCORES.md and PRECHECK_PGQ.md.
  `-pg_scores` is kept as an opt-in flag for the investigation; the default is unchanged. On D, 2,359 of the 3,053
  extra IDs are DIA-NN-confirmed and real. Of the other 694, about 270-500 are false (bounded by DIA-NN's PEP and our decoys),
  so the gain is mostly real but q is somewhat optimistic in the tail. No pg score rewards peak presence over identity
  (Q6_PEAK_IDENTITY.md). The raw 02 "failure" is within one run's counting noise (Q7_NEXT.md). Pooled over raws
  01-03, NEW is at 0.0073 against CUR's 0.0088.
- **HeLa entrapment library and PLAN_PGENT (18 Nov): true FDR on diaPASEF.** Built on codon: a predicted library of
  58,998 peptide pairs. At equal true FDR, `-pg_scores` gives D 9,701 -> 12,033 (+24 %) and held-out H 23,136 -> 26,924
  (+16 %). It fails the absolute bound (FDP_c at q <= 0.01 must be <= 1.2 %; NEW is at 1.66 / 1.86 %), but CUR is just as
  far over (1.76 / 1.67 %). diaPASEF q calibration is the open problem, and it is shared by today's default
  (RESULT_PGENT.md).
- **Is the diaPASEF excess real? (PRECHECK_FOREIGN, 18 Nov):** with Arabidopsis entrapments, the accepted ones are mostly
  near-homologs of human peptides (3-4x enriched for Hamming <= 2). A post-hoc homolog filter suggested near-nominal q,
  but PLAN_PGENT2 removed homologs at library construction and still measured FDP_c 1.43 % (default) and 1.59 %
  (`-pg_scores`) at q <= 0.01. **Correction (BOUNDS.md):** those are combined estimates, which are upper bounds. The lower
  bounds are 0.87-0.99 % in every run, so diaPASEF q is neither shown to fail nor shown to hold (true FDP about
  0.9-1.9 %). The "~1.7 % true FDR" above is an upper bound.
- **Settled with r = 7 (BOUNDS_R5, 18 Nov):** with a homolog-free foreign entrapment at r = 7 the bounds are tight. diaPASEF
  q <= 0.01 gives a true FDP of 1.6-1.8 % (lower bound 1.58 %), the same with and without `-pg_scores`. The calibration
  defect is real and pipeline-wide. `-pg_scores` gives +51 % IDs at 1 % true FDP.
- **Cause found (PAIRWIN, PRECHECK_PDEC, 18 Nov):** pseudo-reverse decoys carry their target's *predicted* RT, IM and
  intensities, which do not fit their own sequence. So among strong candidates a false target beats its decoy 63-64 %
  of the time instead of 50 %. With decoys predicted like targets (DIA-NN style) it is 50 %, and paired-only q <= 0.01
  gives a true FDP of 0.4-1.0 %. Fix: one predicted decoy per target, either in the library generator or in OpenDIA.
  User to choose; then a frozen plan on the r ~7 entrapment.
- **Tool's own q with predicted decoys for 96.6 % of targets (PRECHECK_PDEC "Full coverage", 18 Nov):** the first test
  library left 32 % of targets without a decoy (decoys were reversed from FASTA entries, not digested peptides). Fixed,
  the tool's own q <= 0.01 gives true FDP 0.69-1.13 % (CUR) and 0.74-1.21 % (NEW), D seed 0, r = 1.56. IDs at 1 %
  combined FDP barely move (-3 %). The r ~7 frozen plan is the tight test.
- **Held-out H (19 Nov, RESULT_H.md, reported only):** H reproduces D for both decoy routes. At q <= 0.01 the true FDP is
  1.06-1.27 % (DIALibGen 1.09-1.27 %, OpenDIA 1.06-1.26 %), and IDs at true 1 % FDP are about 20.8 k (CUR) and 24.5 k (NEW).
  With OpenDIA's decoys on H the pairs are fair (win 0.500-0.505), but the excess is unchanged. So the remaining ~1.1x
  is not decoy asymmetry and points at q estimation at r ~7. Proposed next: PRECHECK_TDC on the existing outputs.
- **PLAN_ODP D verdict (19 Nov, 480991c28a): FAIL on criteria 1/2, narrowly.** OpenDIA's own `-decoys predict` at r ~7
  (5 seeds) gives q <= 0.01 a true FDP of 1.06-1.21 %, the same as the DIALibGen route (1.09-1.24 %). Coverage is 99.85 %
  and IDs at true 1 % FDP are 8,233 (CUR) and 10,100 (NEW, `-pg_scores` +22.7 %). OpenDIA alone is as good as the two-tool
  route, and both remove about 85 % of the excess. The remaining ~0.1-0.2 % sits in equal-fragment pairs (pair win 0.54)
  and is not located yet. H runs were submitted after the verdict, with `-trace_tensor_gb 96` (amendment 1: H has
  13.6 M windows, above the default 64 GiB cap).
- **PLAN_PDEC D verdict (19 Nov, 1f393c2184): FAIL on criteria 1/2.** At r ~7 (5 seeds), DIALibGen-predicted decoys give
  q <= 0.01 a true FDP of 1.09-1.24 % (CUR) and 1.11-1.27 % (NEW). With pseudo-reverse decoys it was 1.58-1.80 %. About
  85 % of the excess is gone, but the 90 % intervals exclude 1 %. IDs at true 1 % FDP went from 6,563 to 8,318 (+27 %).
  Pair fairness is 0.54 at r ~7 (0.50 at r 1.56), so a small asymmetry remains. H runs were submitted after the verdict.
- **PLAN_PDEC (frozen 18 Nov, 9feed18a6d):** predicted decoys from DIALibGen for now (user's choice), with OpenDIA
  generating its own decoys as the downstream goal. The test is the r ~7 foreign entrapment. Library: 99.93 % of 1.39 M
  targets paired. D seeds 0-4 are running on codon. Side fix: the PQP reader no longer drops peptides without a gene
  mapping (2982bd647d).
- **OpenDIA predicts its own decoys (DESIGN_OPENDIA_DECOYS, 18 Nov):** `OpenDIAWorkflow -decoys predict` is opt-in
  (01d5ded04a). It replaces a library's decoys with PeptDeep-predicted pseudo-reversed ones, with 99.9 % of targets paired.
  End to end on the r = 1.56 library (D seed 0), q <= 0.01 gives a true FDP of 0.62-1.02 % (CUR) and 0.53-0.86 % (NEW),
  against 1.45 / 1.55 % (upper bound) with the library's own decoys. IDs at true 1 % FDP are within 3-4 %. One seed only;
  the frozen r ~7 plan comes after PLAN_PDEC.
- **Speed and memory, same laptop, Astral raw 01:** 7:18-7:26 against DIA-NN's 7:23. Peak memory is 14.4 GB against
  6.8 GB.

## Astral snapshot, 16 November 2026

### Snapshot against DIA-NN (16 Nov, `gates/single-binary/STATUS_vs_DIANN.txt`)
- **Setup:**
  - Current `OpenDIAWorkflow`: self-calibrating, arm-E search library, `library12_full.pqp` scoring library, seeds 0-4,
    codon on 8 threads. Code `ab46f97ad4`.
  - DIA-NN 2.2.0 on the same 302,480 targets with its own decoys (`--gen-decoys`).
  - Both cut at their own run q ≤ 0.01. Entrapment FDP_c is computed only for evaluation.

| | raw 01 | raw 02 | raw 03 (held out) |
|---|---|---|---|
| OpenDIA N at FDP_c 1 % (median, range) | 1804 (1728-1827) | 1909 (1867-1935) | 1907 (1831-1957) |
| DIA-NN N at FDP_c 1 % (mzML input, the snapshot reference) | 2067 | 2283 | not run |
| OpenDIA / DIA-NN | 87 % | 84 % | |
| DIA-NN N at FDP_c 1 %, **.raw read directly** (2026-11-17) | 2204 | 2356 | 2259 |
| OpenDIA / DIA-NN (.raw) | **82 %** | **81 %** | **84 %** |
| OpenDIA IDs at q ≤ 0.01 (median), FDP_c | 1709, 0.0047-0.0095 | 1811, 0.0012-0.0088 | 1849, 0.0058-0.0096 |
| DIA-NN IDs at q ≤ 0.01, FDP_c | 2361, **0.0144** | 2322, **0.0112** | |
| DIA-NN's q ≤ 0.01 IDs that OpenDIA also reports | 67 % (1588) | 72 % (1676) | |
| Same, human (non-entrapment) IDs only | 68 % | 73 % | |
| OpenDIA's q ≤ 0.01 IDs that DIA-NN also reports | 93 % | 93 % | |
| OpenDIA only / DIA-NN only (q ≤ 0.01) | 121 / 773 | 128 / 646 | |
| Overlap at each tool's own FDP_c 1 % cut | 75 % of DIA-NN's | 76 % of DIA-NN's | |

- **At equal true error (FDP_c 1 %)**, OpenDIA reaches 84-87 % of DIA-NN run on mzML, and **81-84 %** of DIA-NN
  reading the .raw directly. DIA-NN finds 137 (raw 01) and 55-73 (raw 02) more IDs from the .raw than from mzML (its own reader, centroiding as in the
  vendor file). The .raw run is the fairer reference, since OpenDIA reads the .raw too.
- **At each tool's own q ≤ 0.01**, the gap looks larger: 72-78 %. One reason is that DIA-NN's q ≤ 0.01 is
  anti-conservative on this benchmark (FDP_c 1.1-1.4 %), while OpenDIA's is conservative (≤ 1 %).
- **93 % of what OpenDIA reports is also in DIA-NN's list.** The disagreement is mostly DIA-NN IDs that OpenDIA does
  not reach.
- **Where the missed DIA-NN IDs are lost** (current funnel, `gates/extraction-parity/gap/FUNNEL_workflow_library12full.txt`,
  medians, raw 01 / 02):

  | stage | raw 01 | raw 02 |
  |---|---|---|
  | not scored | 8 | 10 |
  | wrong peak | 116 | 106 |
  | right peak, below the cut | 570 | 463 |

  So about 80 % is ranking.

**Run time and resources, one Astral 15-min raw (3.1 GB .raw).** OpenDIA medians come from the 15 codon runs (one per
seed, 3 raws).

| OpenDIA stage | seconds |
|---|---|
| calibration (2 passes) | 128-153 |
| search | 15-19 |
| one evidence pass | 148-240 |
| features | 7-9 |
| scoring | 17-25 |
| **total wall** | **6:08-7:35** |

**Fair benchmark, same laptop (20 cores, 62 GB), raw 01, 8 threads each, run one after the other, both reading the .raw
directly** (`gates/single-binary/bench/`, 2026-11-17). DIA-NN 2.2 in docker gets the .NET 8 SDK through a read-only
mount. Its settings are the same as diann-pairlib-20261113.

| | OpenDIA (a132e9104e) | OpenDIA (5bd8553a74, same output) | DIA-NN 2.2 |
|---|---|---|---|
| Wall time | 10:14 | **7:18** | 7:23 (raw 02 7:44, raw 03 7:06) |
| CPU time (user + sys) | 1060 s (172 %) | 954 s (217 %) | 3391 s (765 %) |
| Peak memory | 15.6 GB | 15.7 GB (14.4 GB since eaa2c95cd0, 7:26 wall) | 6.8 GB (raw 02/03 6.3-6.4 GB) |
| IDs at own q ≤ 0.01 (FDP_c) | 1772 (0.0079) | 1772 (0.0079) | 2278 (0.0123) |
| N at FDP_c 1 % | 1815 | 1815 | 2204 |

- **Speed is now on par with DIA-NN on the same machine (7:18-7:26 against 7:23); memory is about 2.1× higher.**
  - The first fair run (10:14) showed most of the 8 threads idle. A profile found two hot spots, both fixed with
    byte-identical output:
    - The evidence pass tested every RT-overlapping window against all isolation windows per spectrum. Windows are now
      bucketed by precursor m/z (d7feda3dff): 271 s → 111 s.
    - With an MS-level filter, the Thermo reader stopped prefetching peaks, so the main thread decoded every scan
      itself. The read-ahead workers now read peaks for the selected levels too (5bd8553a74): calibration 185 s → 159 s.
  - Stage times now: library read and output writing about 40 s together (the library read is SQLite-bound, single thread), calibration 159 s, search 28 s, evidence pass
    111 s, features 30 s, scoring 69 s (4 fold-parallel threads in the network objective).
  - Remaining items: the calibration pass is bound by the vendor reader (three read-ahead workers, the vendor limit),
    and peak memory (14.4 GB since eaa2c95cd0, which reserves each calibration trace once instead of growing it) is
    set by the 600k-member windowed calibration.
- The earlier "OpenDIA is faster" reading was wrong. In those DIA-NN runs, 4:25 of the 13:27 went to loading a
  5.8 GB mzML, and they ran with other thread counts beside other jobs.
- The same raw 01 seed 0 run took 7:06 on codon (8 threads) against 10:14 here, with the same output. Codon wall
  times therefore look shorter than laptop wall times, and DIA-NN was not timed on codon. Compare only the same-machine rows.

- Peak memory was mis-attributed earlier. The 22 GB peak came from the second (windowed) calibration pass, which collects traces for 600k panel members. The 9.6 GB runs had reused a saved calibration, so the 12-fragment scoring library was not the cause.
  - Fix a132e9104e: the calibration traces no longer store an all-zero 1/K0 array when the library has no mobility, and finished members are trimmed. Local raw01 peak goes from 21.6 GB to **15.6 GB**, and precursors.tsv and calibration.json are byte-identical.
  - The peak is still the calibration pass (15.6 GB). After calibration the run stays at about 9.7 GB.

**Goal.** OpenDIA must identify as many precursors as DIA-NN on the same library, with honest error control and no teacher. The earlier goal was the OpenSWATH + Percolator level, and it is passed.

**Design data.** ProteoBench Astral entrapment raws 01 and 02 are the design runs. Raw 03 is the held-out run: it is evaluated only after a plan's verdict on raws 01/02 is committed (approved by the user on 16 November 2026).

**Metric.** N is the number of target precursors at an entrapment FDP_c of 1 %, with FDP_c = N_e (1 + 1/r) / N. Every count below is a median over five scorer seeds. The seed range is 50-200 IDs.

**One binary: `OpenDIAWorkflow`** (OpenMS `794ee1565e`, README `e612e34d9c`). It runs on Thermo .raw directly, and there is no Python in the analysis. The steps are:
- calibration;
- candidate search;
- one evidence pass, with 24 features (search, boundary, b/y ions, 12-fragment evidence, fragment deficit and shared-signal interference);
- a cross-fitted C++ scorer (a small network with a decoy-whitened readout and self-training);
- pooled paired TDC q-values per run, and a global q-value over runs.

Cost: 7:18 per Astral raw on 8 threads on the laptop (5bd8553a74), level with DIA-NN's 7:23 there; peak memory 14.4 GB (was 22 GB before a132e9104e; set by the windowed calibration pass). Output is byte-identical across repeats and thread counts.

| N at entrapment FDP_c 1 % (median of 5 seeds) | raw 01 | raw 02 | FDP_c at q ≤ 0.01 |
|---|---|---|---|
| OpenDIAWorkflow, reused calibration | 1819 | 1870 | 0.0024-0.0082 |
| OpenDIAWorkflow, self-calibrating, `library12_pop` scoring library | 1718 | 1941 | raw 01 0.0048-0.0104, raw 02 0.0023-0.0045 |
| OpenDIAWorkflow, self-calibrating, full 12-fragment scoring library `library12_full.pqp` (current setup; held-out raw 03: 1907) | 1804 | 1909 | 0.0012-0.0096 on all 15 runs (raws 01-03) |
| Python chain it replaces (D11FDI, C++ interference features) | 1848 | 1925 | ≤ 0.0071 |
| DIA-NN 2.2, same targets, its own decoys (reference) | 2067 | 2301 | |
| Former native base (CNN + LDA) | 1029 | 1298 | |

- **Decoys.** Since PLAN_DECOYGEN, the decoys come from OpenSwathDecoyGenerator `-predict_decoys`, which gives each decoy its own predicted spectrum and RT. This removed the entrapment excess that failed every earlier arm.
- **Interference features** (PLAN_GAP arm I). These gave +118 to +185 on raw 01 and +199 to +351 on raw 02 at 1 % over seeds. About a third of the raw 02 gain comes from entrapment pairs that only exist in this benchmark.
- **Self-calibrated RT window.**
  - On raw 01 it is wide (151 s).
  - A refined spread fit gets the spread right, but over five seeds it scores +5.8 % on raw 01 and -3.1 % on raw 02, so it fails (PLAN_CALREFINE2).
  - Two label-free width rules are undefined or point the wrong way.
  - A cross-fitted per-run width chosen on IDs (PLAN_RTWIN2) scores −2.2 % on raw 01 and −2.3 % on held-out raw 03, so it fails. That post-hoc filter was removed.
  - Recalibrating from the run's own IDs, with a window that follows the local RT error (PLAN_RTPROFILE), recovers the gradient ends: +23 to +24 IDs early, +7 to +17 late. But the narrower middle loses (raw 02 −52), the strict-FDR ranking is less stable, and it fails: raw 02 −1.05 %, raw 01 seed 4 FDP_c 0.0104, held-out raw 03 −0.6 %. Removed (`gates/single-binary/rtwin/`).
  - With the full 12-fragment scoring library, the wide window keeps FDP_c ≤ 0.0096 on all 15 runs, so the width is no longer an FDR problem.
- **Gap to DIA-NN.** It is 13-16 % at FDP_c 1 %, almost entirely in score formation on the right peak (`gates/extraction-parity/gap/`).
  - A joint (presence-weighted) interference fit has no room. Block I's explaining-away comes from claimants that are themselves absent (PRECHECK_JOINT).
  - A labelled ceiling on today's 24 inputs recovers only about +100 per raw at today's decoy rate.
  - Scorer-capacity attempts on those inputs all failed:
    - signal-level interactions and self-trained trees, offline;
    - a 3-seed ensemble (PLAN_ENSEMBLE): +0.2 % / +0.6 %, raw 03 +0.3 %. It does cut the raw 01 seed range from 99 to 17.
  - Closing the gap needs new evidence for the weak right-peak IDs.
  - Wrong- and missed-peak losses (about 100 per raw) have no cheap fix (PRECHECK_PEAK).
  - Run-specific fragment choice has no room (PRECHECK_FRAG). The searched 6 already show signal in 5.4-5.5 of 6 for the rejected right peaks, and re-selecting fragments from the run's own IDs leaves the separation from decoys at AUC 0.71-0.72. The weak IDs are low in abundance across all fragment ranks, not missing particular fragments.
  - MS1 isotope evidence (PLAN_MS1) is the first new block since block I that raises N at equal true FDR: +5.7 % / +2.1 %, every raw 01 seed above every CUR seed.
    - It fails its frozen plan. Its q ≤ 0.01 count stays flat on raw 01, and entrapments in the rejected band rank slightly higher (0.531 against the 0.527 limit).
    - Held-out raw 03 does not confirm the gain (−0.5 % at FDP_c 1 %, −2.1 % at q ≤ 0.01). It is removed per the plan (code kept in git history, 51fb7cc7ec).

**Direction.** The single voxel model was replaced in practice by the heuristic-feature network (D11 line, R-WP3). That network is now in C++ as `OpenDIA::CrossFitScorer`. OpenDIACNN and LibTorch are no longer on the default path.

## Why the direction changed

The legacy plan kept OpenSWATH's chromatogram and feature-scoring path and tried to make it cheaper with caches, tiles and admission control.
Its own experiments showed that this did not change the cost structure. The frame cache, tiled extraction and scorer micro-optimisations did not give a faster engine.
See the archived [redesign decision](archive/legacy-wp-plan-20261107/redesign-decision-20261103.md).

The rewrite replaces the whole dataflow:

```
raw run (.d / mzML / Thermo .raw)
   │
   ▼
OpenDIACalibrate      identification-free RT, fragment m/z and 1/K0 calibration (no iRT kit, no IDs)
   │  calibration.json, anchors.tsv
   ▼
OpenDIA -calibration  streaming sparse-grid search of the full library inside calibrated windows
   │  candidate peak groups (multi-apex), search features
   ▼
evidence pass         6-fragment traces, apex b/y ions, 12-fragment grids, shared-signal interference graph (24 features)
CrossFitScorer        decoy-whitened network, family cross-fitting, self-training, pooled paired TDC → q-values
                      (all inside OpenDIAWorkflow; OpenDIACNN is the earlier, optional CNN + LDA scorer)
```

## Current HeLa result (design run, not independent evidence)

**Current (2026-11-17 evening; RESULT_RTCENTRE.md, RESULT_IMGRID.md).** Three diaPASEF fixes, each frozen-plan tested on
HeLa D (2737.d) with a held-out HeLa run H (22474.d), and each byte-identical on Astral:
1. calibrated 1/K0 in the evidence pass (0e9e4fc3a5);
2. one-pass calibration stage 3 centred on the stage-1 map (6dd893b954);
3. mobility-resolved 12-fragment grid and interference (963543c7f0).

| HeLa, own q ≤ 0.01 (median of 5 seeds) | D 2737.d | H 22474.d (held out) |
|---|---|---|
| this morning | 11,653 (seed 0) | 19,725 (with fix 1) |
| now | **16,116** | **24,000** |
| DIA-NN 2.2 | 28,950 | 38,074 |
| OpenDIA / DIA-NN | 56 % | 63 % |
| OpenSWATH + PyProphet SVM | 19,132 (OpenDIA at 84 %) | not run |

About 93-95 % of OpenDIA's IDs are in DIA-NN's list on both runs. The remaining gap is mostly scoring; the b/y-ion
floor and a 12-fragment HeLa library are the next candidates (DIAG_HELA.md).

**Update 2026-11-17 (`gates/single-binary/bench/RESULT_HELA.md`).** The current OpenDIAWorkflow gets **11,653** precursors at
q ≤ 0.01 on this run. DIA-NN 2.2 gets 28,950 on the same library, OpenSWATH + PyProphet SVM 19,132, and the 2026-11-07
rewrite below 13,253. So the branch reaches about 80 % of DIA-NN only on Astral. On diaPASEF it is at 40 % of DIA-NN
and below its own earlier HeLa result: the Astral-driven changes regress this run. Wall time is 3:08 against 2:46, and
peak memory 4.5 GB against 2.3 GB.

This section records the state of 7–8 November 2026. For the current Astral state, see "Current state" above and [BOARD.md](BOARD.md).

| HeLa 200 ng diaPASEF, paired TDC, experimental q ≤ 0.01 | Targets | SVM-baseline overlap |
|---|---|---|
| OpenSWATH + PyProphet SVM (reference) | 19132 | n/a |
| OpenSWATH + PyProphet XGBoost (reference) | 19725 | n/a |
| Rewrite, no RT prior | 6531 | 5846 |
| Rewrite, CNN-anchor RT window | 8554 | 7502 |
| Rewrite, ID-free calibration filters candidates after extraction | 10679 | 9146 |
| **Rewrite, calibrate then extract (current)** | **11830** | **10158** |
| **Rewrite, teacher-free: panel-40000 calibration anchors + `-ss_iterations 1`** | **13611** | **11533** |
| **Rewrite, teacher-free end to end from one panel-40000 calibration (frozen design)** | **13253** | **11275** |

The current pipeline extracts 198,915 candidates in 40 s and runs CNN crossfit in 1.2 min with a 2.1 GB peak.

**Limits:**
- These are decoy-count estimates on the run used to design the pipeline. They are not a calibrated FDR.
- Every row except the last two uses CNN seeds from the same-run OpenSWATH teacher. The last two use no OpenSWATH input and exceed the teacher-seeded count. Neither has been checked by entrapment.
- On the independent 500 pg run, the CNN alone gave 0 IDs at 1%, paired TDC gave 230, and the calibration failed closed with 0 anchors.
- On the independent ProteoBench Astral entrapment run (raw 01, frozen teacher-free design applied once), the calibration also failed closed with 0 anchors, so no IDs and no entrapment FDP. The HeLa design does not yet transfer to weak or no-1/K0 runs.
- Option B then used raw 01 as a design run, adding a windowed second calibration pass and pooled paired TDC (both opt-in). The frozen result was evaluated once on raw 02 against a hashed gate, and it **failed**.
  - It gave 604 anchors and 360 targets at q ≤ 0.01, with 5 entrapment hits.
  - Combined entrapment FDP estimate 0.028 against a 0.015 limit (95% interval 0.009 to 0.065).
  - On this weak Astral run the experimental q is not conservative enough. See R-WP2 and R-WP4.
- An opt-in second-stage LDA (`OpenDIACNN -stage2`) combines the CNN score with co-elution, library correlation and MS1 evidence.
  - Raw 01 design run: 333 targets at q ≤ 0.01, against 321 for CNN-only.
  - Spent raw 02 re-score (a diagnostic, not a gate): 518 targets, FDP_c 0.012.
  - It has not been evaluated independently. Raw 03 is held back for that.
- The rewrite is about 5,500 to 6,100 precursors short of the reference. Closing that gap is R-WP5.

## Work packages

| WP | Scope | State, 16 November 2026 |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming candidate search | Delivered: native-scan XIC peak picking and Thermo `.raw`. Recall at calibrated windows is open. |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration | Delivered, with the windowed pass. **Open: the self-calibrated RT window width** (raw 01 too wide). |
| [R-WP3](work-packages/R3-trace-cnn.md) | Scorer | **C++ `CrossFitScorer` in `OpenDIAWorkflow`**: 98.6 / 99.0 % of the Python D11FDC median. CNN + LDA stopped. |
| [R-WP4](work-packages/R4-error-control.md) | Error control | q ≤ 0.01 holds with reused calibration (FDP_c ≤ 0.0082). Self-calibrated raw 01 reaches 0.0104 on 2 of 5 seeds. |
| [R-WP5](work-packages/R5-id-gap.md) | ID gap | **Active.** Past OSW + Percolator. 88 / 81 % of DIA-NN (raw 01 / 02, reused calibration). |
| [R-WP6](work-packages/R6-acceptance.md) | Acceptance, packaging, merge | One-command tool delivered. Raw 03 held back. PR split open. |

Dependency order: R-WP5 uses R-WP1 to R-WP4. R-WP6 needs R-WP4 and R-WP5.

**Step board:** [BOARD.md](BOARD.md) lists each step with its state (done, in progress, next, blocked) and the link to its record.

## Rules carried over

- Defaults stay byte-identical. New behaviour is opt-in until it is accepted. Check HeLa default equivalence in a fresh directory before each benchmark run.
- Write the pass rule in a PLAN file (with sha256 and time) before each run or training. Do not change a rule after a result.
- Design on raws 01 and 02 only. Do not tune on an evaluation run. Raw 03 needs approval.
- Entrapment labels never go into training. DIA-NN and OSW are references only. Do not reuse OSW scoring code.
- Never report a calibrated FDR from decoy counts alone. Call q-derived counts decoy-count estimates. Agreement with OpenSWATH is not truth.
- Native Bruker and Thermo input, with no mandatory conversion.
- Inputs and baseline OSW files are read-only.

## Issues

The epic is [#10661](https://github.com/OpenMS/OpenMS/issues/10661). The legacy WP issues #10663 to #10671 and #10675 are closed, each pointing to this plan. #10662 (WP-00) was closed earlier.
See [issue-posting.md](issue-posting.md) for the issue of each R-WP.
