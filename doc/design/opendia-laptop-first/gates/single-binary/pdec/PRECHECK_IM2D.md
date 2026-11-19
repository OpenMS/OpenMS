# PRECHECK_IM2D: RT × ion-mobility candidate picking (19 Nov 2026, diagnostic, seed 0 only)

**Code.** OpenMS 07131138b8, `OpenDIAWorkflow -im_peaks <bins>`, opt-in.
- Each (scan, fragment) keeps its most intense point per 1/K0 bin. Only the top 4 bins are kept, and the pooled
  maximum is unchanged, so the 1-D candidates are identical.
- Local maxima of the smoothed scan × bin consensus seed a pick on the 3-bin mobility band around them.
- The band peak becomes an extra candidate when no existing candidate apex is within one scan. The limits are
  `-im_extra_peaks` per precursor and `-im_min_support` fragments.
- The unit test builds a target masked in the pooled trace by a precursor at another mobility, and checks that it is
  recovered.
- The default is byte-identical to before, and an `-im_peaks` run is reproducible byte for byte. 25/25 ctests pass.

**Data.**
- r ≈ 1.56 HeLa D, seed 0, local, 8 threads, current default scorer.
- DIA-NN IDs (25,881 non-entrapment, q ≤ 0.01) and entrapment labels are used for evaluation only.

| | default | `-im_peaks 8` (≤ 2 extra, support ≥ 1) | `-im_peaks 8 -im_extra_peaks 1 -im_min_support 3` |
|---|---|---|---|
| candidate windows | 2.39 M | 3.52 M (+47 %) | 2.84 M (+19 %) |
| DIA-NN IDs with a right-RT candidate | 84.1 % | **93.5 %** | 91.6 % |
| ... lowest-abundance quartile | 77.7 % | 90.7 % | 87.3 % |
| DIA-NN IDs gaining a right candidate (none lose one) | | 2,439 | 1,936 |
| ... of these accepted at q ≤ 0.01 | | 13 % | 16 % |
| N at own q ≤ 0.01 | 13,034 | 13,220 | 13,368 |
| **N at FDP_c 1 %** | **12,357** | **12,788 (+3.5 %)** | 12,591 (+1.9 %) |
| lower / combined bound at q ≤ 0.01 | 0.79 % / 1.49 % | 0.70 % / 1.32 % | 0.73 % / 1.37 % |
| pair win | 0.524 | 0.521 | 0.505 |
| search consume | 14.3 s | 22.8 s | 35.4 s* |
| interference resolver | 40.6 s | 75.2 s | 62.6 s |
| wall / peak RSS | 5:55 / 10.9 GB | 7:15 / 14.3 GB | 7:47 / 12.1 GB* |

\* The lean run overlapped another local job, so its times are not clean.

## Reading
- **The picker does what it was built for.** Of the 15.9 % of DIA-NN IDs with no right candidate (PRECHECK_IMPICK),
  about 60 % now have one, most of them at low abundance. None loses a right candidate.
  - Calibration does not get worse: the true FDP at q ≤ 0.01 is slightly lower.
  - Pair fairness is unchanged.
- **IDs rise little (+3.5 % at true 1 %).** Only 13-16 % of the recovered right peaks are accepted. Finding the peak was
  not enough: these are faint and still score low, the same bottleneck as the 30 % "right peak, low score" group.
- **Cost:** +47 % candidate windows, +1:20 wall (+23 %), +3.4 GB.
  - Most of the extra time is downstream of the search (interference +35 s, extraction +30 s), not the 2-D pick.
  - The search itself costs +8.5 s.
  - About half of the added windows are decoys. That is expected, since targets and decoys are picked alike.
- **Single seed and only r ≈ 1.56.** The +3.5 % is within seed noise of earlier plans (±1-2 % per seed). The r ≈ 7
  library, where crowding is 4.5× denser, is where it should matter more.

## Next
- A frozen 5-seed plan on r ≈ 7 D with `-im_peaks 8`, against the current default, is the test.
- The candidate inflation is worth trimming first, because it drives the downstream cost: about half the added
  windows have support ≤ 3 and 13 % of them are right.
