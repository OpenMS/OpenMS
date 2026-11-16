# PLAN_MS1: MS1 isotope evidence block in OpenDIAWorkflow

Written 16 November 2026, before the implementation and before any run. The evidence is in
`gates/extraction-parity/gap/PRECHECK_MS1.md`.

## Block (5 features appended to the 24 scorer inputs; label-free, identical for targets and decoys)
- **Traces.** For every population window, the evidence pass collects M+0..M+3 MS1 traces on the existing 32-bin grid
  over apex ± 15 s.
  - M+i sits at the precursor m/z + i × 1.0033548 / z.
  - Each trace sums all points within ±10 ppm of the target shifted by the run's calibrated MS1 offset (the MS2 offset
    when there is no MS1 fit), as `build_voxels.py` does.
  - The fragment-sum trace is the sum of the existing 6 search-fragment traces of the window's tensor, summed over
    mobility.
- **Theoretical envelope.** M+0..M+3 come from the precursor's elemental composition: residues plus water plus the
  UniMod:4 / UniMod:35 deltas, with C/H/N/O/S heavy isotopes by convolution, as `ms1_precheck.py` does. An unknown
  modification gives an all-zero envelope, and the cosine is then 0.
- **Features:**
  1. `m1_corr_xic` = Pearson(M+0, fragment sum) over bins 6..25.
  2. `m1_iso_cos` = cosine(observed central areas M+0..3 in bins 12..19, theoretical envelope).
  3. `m1_iso_corr` = mean of Pearson(M+0, M+1) and Pearson(M+0, M+2) over bins 6..25.
  4. `m1_apex_dt` = |argmax bin of M+0 − argmax bin of fragment sum| over bins 6..25.
  5. `m1_log_m0` = log1p(central M+0 area, bins 12..19).
  - Without MS1 signal in bins 6..25, features 1-3 are 0 and feature 4 is 20.
- The bins and constants are those of the precheck (bins 6..25 is the precheck's 6:26 slice). Nothing is tuned after
  the precheck.

## Arms, runs, criteria
- **Arms:** CUR = today, with 24 inputs. NEW = `-ms1_features`, with 29 inputs. Both are OpenDIAWorkflow runs:
  self-calibrating, arm-E `library.pqp` plus `library12_full.pqp`, seeds 0-4, codon, 8 threads.
- **Runs:** raws 01/02 are the design runs. Raw 03 is held out and evaluated once, after the verdict is committed.
- **Evaluation:** `eval_run.py`.
- **Criteria**, on raws 01 and 02 separately:
  1. NEW median N@FDP_c 1 % ≥ CUR median − 1 %. MS1 must not cost IDs at equal true FDR.
  2. NEW FDP_c at q ≤ 0.01 ≤ 0.010 on every seed.
  3. NEW median count at q ≤ 0.01 ≥ CUR median + 2 % on both raws. This is the main expected effect: a less conservative
     q at controlled error.
  4. Specificity: on seed 0, the AUC of entrapment-in-band winners against their matched decoys under the NEW score is
     ≤ that AUC under the CUR score + 0.02. This is computed with the same grouping as the precheck.
- **Outcome:** pass makes `-ms1_features` the default. Fail removes the code. Nothing is retuned.

**Benchmark caveat (to be stated with any gain).** Entrapments share their parent's m/z and envelope, so MS1 evidence
for an entrapment is partly its parent's. The precheck's E > ED AUC (0.50-0.57) bounds that effect, and criterion 4
guards it.
