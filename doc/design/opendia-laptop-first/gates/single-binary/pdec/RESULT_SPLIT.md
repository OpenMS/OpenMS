# RESULT_SPLIT: peak splitting rarely cuts a real peak; cheaper split arms (21 Nov 2026, diagnostic)

The gate is in `PRECHECK_SPLIT.md` (`a12cd8b8d1`). The part A script is `precheck_split_a.py` (sha256 `90f24fad…`). The full output is in `EVAL_SPLIT_A.txt`.

## Part A: does V8M8 split real peaks?

The common set is HeLa IDs accepted by def, V8M8 and DIA-NN, with def's apex within 6 s of DIA-NN's.

| | D r0 | D r7 | H r0 | H r7 |
|---|---:|---:|---:|---:|
| Common set | 17,130 | 14,392 | 31,251 | 28,500 |
| **Split share** (same apex ±3 s, < 0.7× width, nested) | **0.42 %** | **0.43 %** | **0.85 %** | **0.87 %** |
| Split IDs: DIA-NN extent covered, def → V8M8 | 1.00 → 0.76 | 1.00 → 0.86 | 1.00 → 0.82 | 1.00 → 0.79 |
| Split IDs: area ratio V8M8 / def (median) | 0.48 | 0.54 | 0.56 | 0.57 |
| All common: DIA-NN extent covered, def → V8M8 | 0.718 → 0.717 | 0.831 → 0.829 | 0.917 → 0.909 | 0.917 → 0.909 |
| Area-quantity Spearman, def → V8M8 | 0.702 → **0.708** | 0.707 → **0.717** | 0.784 → **0.796** | 0.785 → **0.798** |
| Lost IDs: share that are splits | 1.8 % (of 329) | 3.3 % (of 554) | 1.9 % (of 861) | 2.2 % (of 536) |
| Gained IDs: true recoveries | 61.2 % (of 1,000) | 67.2 % (of 828) | 77.3 % (of 1,737) | 66.4 % (of 1,597) |
| Median width, def / V8M8 / DIA-NN (s) | 5.0 / 5.0 / 6.9 | 5.9 / 5.9 / 6.9 | 8.6 / 8.6 / 8.6 | 9.5 / 8.6 / 8.6 |
| **Rule A** | pass | pass | pass | pass |

**How we know a real peak is not being split, in plain terms**

- **Splits are rare.** Fewer than 1 % of accepted IDs end up with a reported window that is a nested half of the default one. When a split does happen it is visible: the window covers only about 80 % of DIA-NN's peak, and the area halves. So the damage is bounded and measurable.
- **Quantification gets better, not worse.** Area tracks DIA-NN's quantity more closely with V8M8 on every cell (Spearman +0.006 to +0.013). Splitting at a real dip removes co-eluting interference from the area. That effect outweighs the rare accidental cut.
- **The IDs lost are not lost to splitting.** Only 2-3 % of the IDs V8M8 loses are split windows. The rest are competition effects, since there are more candidates per precursor.
- **Most gained IDs are true recoveries.** 61-77 % of new IDs have V8M8's apex on DIA-NN's peak where def's apex was elsewhere, so the right peak was hidden in a merge.
- **Widths line up.** On H r7 the median V8M8 width (8.6 s) moves to DIA-NN's (8.6 s), down from def's 9.5 s. Def was over-merging.

So no guard is needed for V8M8, and the same measures are applied to the part-B arms.

## Part B: cheaper split arms

The evaluator is `precheck_split_b.py`, which is `precheck_candrec_b.py` with the arm list extended and nothing else changed. Split harm comes from `precheck_split_a.py`. All 8 new runs exited 0. The full output is in `EVAL_SPLIT_B.txt`.

| Arm | D r0 HeLa | H r0 HeLa | D r7 N@FDP_c 1 % | H r7 N@FDP_c 1 % | r7 lower Δ (D / H) | Wall (D r0 / D r7 / H r0 / H r7) | Split share (max over cells) |
|---|---:|---:|---:|---:|---|---|---:|
| def | 19,082 | 34,702 | 15,465 | 30,543 | | 1 / 1 / 1 / 1 | |
| V8M8 | +4.2 % | +2.8 % | +3.2 % | +3.7 % | −0.14 / −0.05 pp | 1.45 / 1.78 / 1.24 / 1.15 | 0.87 % |
| V7M6 | +4.5 % | **+1.3 %** | **+3.5 %** | **+3.2 %** | +0.04 / −0.05 pp | 1.19 / 1.39 / 1.10 / 0.73 | **0.15 %** |
| V8M6 | +6.0 % | **+1.9 %** | **+4.2 %** | **+3.1 %** | −0.15 / +0.04 pp | **2.11** / **1.54** / 0.65 / 0.71 | 0.85 % |

- **Rule A criteria 1-3 hold for both new arms on every cell.**
  - The split share is at most 0.85 %. V7M6 is 0.08-0.15 %, five times less than V8M8.
  - The area-quantity Spearman rises by +0.003 to +0.013.
  - Lost IDs that are splits are 0.2-2.2 %.
  - Criterion 4 (gained true share ≥ 60 %), which rule B does not require for the arm, is 54-78 %.
- **Rule B fails for both arms on H r0:**
  - V7M6 gives +1.3 % and V8M6 +1.9 %, against the +3 % bar.
  - V8M6 also fails wall time (2.11× on D r0, 1.54× on D r7).
  - V7M6 passes every other criterion: D r0 +4.5 %, r7 +3.5 / +3.2 %, lower bound within +0.04 pp, wall time at most 1.39×.
- **Verdict: no arm qualifies under rule B as written.**
- **Wall time is noisy at n = 1.**
  - H r0 V8M6 ran in 0.65× the time with 0.66× the CPU of def, which is not possible for a superset search, so nodes differ.
  - CPU-time ratios (user seconds) on D, the stable run, are steadier. V7M6 is 1.21× (r0) and 1.58× (r7). V8M8 is 1.47× and 2.11×.
  - Any plan needs a fixed node type, or a repeat, before wall time can be judged.

## What the four split arms agree on

- **At r7, the frozen accuracy metric, every valley-splitting arm gains +3.1 to +4.2 % on both D and H**, with calibration unchanged or better.
- **The cap-only arm (M10) gains nothing.** The gain comes from splitting merges, not from more candidates.
- **The H r0 shortfall is consistent.** H r0 is the run where cap crowding dominates, and a split only helps once the right peak fits in the slots. V8M8, with 8 slots, does best there (+2.8 %). With 6 slots the split peaks compete for fewer places.
- **The safest arm is V7M6.** It has the fewest accidental splits (≤ 0.15 %), the lowest cost, and the most even r7 gain. Its only failure is H r0, which is a HeLa count at its own q with no entrapment control. That count is not the metric that decides a default.

A frozen 5-seed r7 plan for V7M6 is the natural next gate. That decision belongs to the user: this precheck's own rule says no plan, and changing the bar after seeing the result is not allowed. A plan would test V7M6 on what is actually decided (r7 N@FDP_c 1 %, D first), with the wall-time bar measured on a fixed node type.
