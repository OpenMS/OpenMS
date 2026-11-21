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
