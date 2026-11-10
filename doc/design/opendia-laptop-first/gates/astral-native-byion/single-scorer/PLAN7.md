# Round 7: transition refinement in peak picking and quantification

Written before any round 7 run, after `RESULT6.md`, `metric_noise.txt` and `interference_impact.txt`. Code: `dfbe38b470`. Raw 03 stays untouched.

## Why

Round 6 showed that refinement separates real from false windows well (mean r 0.89–0.91 for accepted IDs against 0.29 for decoy winners). It failed when used as extra LDA features, and when used as a mask before the CNN. EncyclopeDIA uses refinement for two things:
- choosing which fragments, and so which peak boundaries, describe the peptide;
- quantifying from interference-free fragments.

This round tests those two uses.

## Arm P1 (raws 01 and 02)

P1 is the native base command (`-stage2 -stage2_feature_set search+boundary+byion`, CNN grid2d), plus:
- `-stage2_refined_picking_min_r 0.75`. The boundary group re-picks each window's peak group on the fragments with r ≥ 0.75. Features are still measured over all fragments. This is the only change to scoring.
- `-quant_out`. Output only, with no effect on scores.

The CNN and the search, byion and extraction inputs are unchanged.

**Order of runs:**
1. HeLa default equivalence. `crossfit.tsv` must be byte-identical, or the chain stops.
2. P1 on raw 01, then P1 on raw 02.

**Internal check (part of rule 4).** These must be byte-identical to the base run:
- `crossfit.tsv` columns 1–6 and 9;
- `features.tsv` columns 1–9, which are the ids, `cnn_z` and the search features;
- `features.tsv` columns 15–16, which are byion.

## ID decision rule

This is `PLAN.md` rules 1–4, unchanged, with the score `stage2_score`:
1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes with peak RSS ≤ 30 GB and wall time ≤ 35 min, and the internal check holds.

`metric_noise.txt` showed that this rule is noise-limited. Changing the rule needs approval, so the rule stays as it is. In addition, the paired bootstrap difference against base (90 % interval, N@1% and N@2%, as in `metric_noise.py`) is reported. A clear loss is when the 90 % interval at 2 % lies entirely below 0 on either raw. A clear loss is reported as such even if the rule passes.

## Quantification rule (P1 quant_out on both raws)

**Precursor set.** Fixed before the runs: the precursors that base accepts at FDP_c 1 % on both raws and that are not entrapment. Each precursor is quantified at its base winning window (`sample_id`), so the window choice is the same for both methods.

**Methods.**
- **Naive:** `naive_area`, which is all fragments over the first-pick boundaries.
- **Refined:** the sum of `area<c>` over the fragments with r ≥ 0.9 in both raws. At least 2 fragments are required; otherwise the precursor is treated as missing for the refined method.

**Metric.** For each method:
- log2(raw 01 / raw 02) per precursor;
- centre by the median;
- robust SD = 1.4826 × MAD.

Only precursors that both methods quantify are used, so both methods are compared on the same set.

**Pass.** All three must hold:
1. The robust SD of the refined method is ≤ 0.95 × the robust SD of the naive method. This is a reduction of at least 5 %.
2. The paired bootstrap 90 % interval of (SD_refined − SD_naive), over precursors, lies entirely below 0.
3. The refined method quantifies ≥ 80 % of the precursor set.

**Reported (reference only, never used in a rule).**
- The Spearman correlation of each method with DIA-NN `Precursor.Quantity` in raw 01, over the shared precursors.
- The DIA-NN 01/02 robust SD (0.335 over 72k precursors).
- The share of precursors whose refined boundaries differ from naive, and the change in median width.
- The histogram of the number of quantitative fragments.

## Mechanism (reported only)

- The entrapment-vs-decoy winner AUC of `boundary_log_width` and `boundary_coelution`, base against P1.
- For the base 1 % IDs: the share of windows where refined picking changed the boundaries.
