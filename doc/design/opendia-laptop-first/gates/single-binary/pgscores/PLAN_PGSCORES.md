# PLAN_PGSCORES: OpenSWATH-style peak-group scores in OpenDIAWorkflow

Written 18 November 2026, before the implementation and before any run. Evidence: `bench/GAP_HELA_20261118.md` and
`bench/ABLATE.txt`. Those are labelled diagnostics on OpenSWATH's own peaks, so they only rank the groups.

## Block (10 features appended to the 24 scorer inputs; label-free, identical for targets and decoys)
Behind the flag `-pg_scores`. Without the flag the output must stay byte-identical.

**Inputs.**
- The 6 search-fragment traces of each window's tensor, summed over mobility, on the scan bins.
- The picked peak `[left, right]` from the existing boundary picker (`pickWindow`).
- Weights `w_k`: library intensities of the non-padding channels, normalised to sum 1.
- Observed share `a_k`: fragment area inside the peak divided by the total over fragments.

**Features 1-4: cross-correlation.** These follow OpenSWATH's `MRMScoring` over the peak bins:
- Each trace is standardised (population sd).
- For every pair i ≤ j (diagonal included) we take the normalised cross-correlation at lags −D..D, with
  D = min(10, span). The maximum, first on ties, gives a lag and a value.
- The four features:
  1. `pg_xcorr_coelution` = mean + sample sd of |lag|.
  2. `pg_xcorr_coelution_w` = Σ_i w_i² lag_ii + 2 Σ_{i<j} w_i w_j |lag_ij|.
  3. `pg_xcorr_shape` = mean of the values.
  4. `pg_xcorr_shape_w`, weighted like 2.
- Without a picked peak of ≥ 3 bins and ≥ 2 channels with signal: 1-2 = 10, 3-4 = 0.

**Features 5-6: signal.**

5. `pg_log_sn` = log1p(mean_k apex_k / (median_k + floor)).
   - apex_k is fragment k at the consensus apex bin.
   - median_k is its median over all scan bins.
   - floor is the run's median nonzero per-bin fragment intensity, from up to 20,000 evenly strided windows. It is
     derived from the file.
6. `pg_intensity_fraction` = signal of all fragments inside the peak divided by signal over all scan bins.

**Features 7-10: core-window spectrum evidence.**
- The evidence pass collects these only with the flag.
- They come from spectra within the existing core half width (4 s) of the apex, within the calibrated precursor 1/K0
  tolerance window of the tensor.
- Tolerance t = the run's calibrated `ms2_tolerance_ppm`, around the calibrated fragment m/z.
- Per fragment k we sum:
  - the monoisotope intensity m_k and its intensity-weighted signed ppm;
  - isotopes +1..+4 at +i × 1.0033548 / z_frag, where z_frag is the library fragment charge (1 if unset);
  - preceding peaks at −1.0033548 / c for c = 1..4.

7. `pg_massdev` = mean_k min(1, |ppm_k| / t). A fragment without signal counts as 1.
8. `pg_massdev_w` = Σ_k w_k min(1, |ppm_k| / t).
9. `pg_isotope_corr` = Σ_k a_k Pearson(observed M+0..M+4, averagine(mz_k × z_frag) scaled to max 1).
   - The averagine comes from `CoarseIsotopePatternGenerator`, tabulated per 1 Da.
   - A fragment with no signal scores 0.
10. `pg_isotope_overlap` = Σ_k a_k × #{c : preceding peak at charge c > m_k}.

These are OpenSWATH's definitions (dia_nr_isotopes 4, dia_nr_charges 4) with our extraction geometry. Nothing is tuned
after this plan.

## Arms and runs
- **CUR** = tip 963543c7f0, which has no flag. **NEW** = the same code with `-pg_scores`.
- **HeLa D (design):** NEW seeds 0-4, 4 threads, two at a time. The CUR seeds are reused from PLAN_IMGRID (`D_imgrid_s*`,
  median 16,116, agreement 0.9299).
- **Astral raws 01/02 (design):** NEW seeds 0-4 on codon, 8 threads, node-local HOME, in a container of the NEW commit.
  The CUR seeds are reused from PLAN_MS1 CUR. Tip output on Astral is byte-identical to those runs (A_imgrid s0
  = 1772 at q ≤ 0.01).
- **Held out**, evaluated only after the design verdict is committed:
  - HeLa H, NEW seeds 0-4, against CUR `H_imgrid_s*` (median 24,000);
  - Astral raw 03, NEW and CUR seeds 0-4.

## Criteria
**Gate 0 (precheck).** Seed 0 on D with `-features_out`, before any other NEW run.
- The labelled GBM ceiling on the common set (the `ablate.py` method) uses our 24 features plus the 10, against the
  24 alone.
- It must rise by ≥ +300.
- Otherwise stop: no further runs, and the code is removed.

**Design.** All must hold:
1. D: NEW median N ≥ 1.03 × CUR median (≥ 16,600).
2. D: NEW median agreement with DIA-NN ≥ CUR median − 0.005 (≥ 0.9249).
3. Astral raws 01 and 02, each:
   - NEW median N@FDP_c 1 % ≥ CUR median − 1 %;
   - NEW median FDP_c at q ≤ 0.01 ≤ 0.010;
   - no NEW seed has FDP_c above 0.0125.

**Outcome.**
- Pass: `-pg_scores` becomes the default and the flag is removed.
- Fail: the code is removed. Nothing is retuned.
- Held-out runs are reported in either case: H median against 24,000, agreement, and raw 03 N@FDP_c 1 % against CUR.

**Caveats to be stated with any result.**
- diaPASEF has no entrapment, so agreement with DIA-NN is only a proxy for specificity.
- The OpenSWATH ablation ranks evidence measured on OpenSWATH's own peaks; our features use our picked peak and binned
  traces.
