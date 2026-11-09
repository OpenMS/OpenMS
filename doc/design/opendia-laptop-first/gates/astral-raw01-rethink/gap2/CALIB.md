# Why decoy-estimated q is anti-conservative (raws 01 and 02, OSW-like scorer)

This is a design diagnostic. Entrapment is used here only to diagnose. Any choice made from it must be confirmed on an untouched run (raw 03).

The scripts are `calib.py`, `calib_curve.py`, `calib_tail.py`, `calib_borrow.py`, `calib_hits.py`, `calib_il.py`, `calib_shift.py`, `calib_struct.py`, `calib_entr.py`, `calib_entr2.py`, `calib_lowmass.py`, `calib_y1.py`, `calib_scan.py`, `calib_libreal.py`, `calib_mitig.py` and `matched_fdp02.py`. Each writes a matching `.txt` output.

## 1. Size and shape (`calib_curve.txt`)

Decoy-estimated FDR against entrapment FDP_c along the accepted list. The 90% intervals are exact Poisson.

| run, scorer | q ≤ 0.01: N, FDP_c [90% CI] | q ≤ 0.05: N, FDP_c [90% CI] |
|---|---|---|
| 01 LDA | 993, 0.010 [0.004, 0.021] | 1429, 0.058 [0.044, 0.075] |
| 01 MLP | 1222, 0.012 [0.005, 0.022] | 1691, 0.078 [0.062, 0.095] |
| 02 LDA | 1115, 0.024 [0.014, 0.037] | 1463, 0.079 [0.062, 0.098] |
| 02 MLP | 1231, 0.020 [0.011, 0.032] | 1573, 0.072 [0.057, 0.090] |

**The excess is systematic, not noise.** At q ≤ 0.05 every interval lies above 0.05 on both runs, for both learners. At q ≤ 0.01 it is significant on raw 02 only. OSW+Percolator's own m_score ≤ 0.01 export has FDP_c 0.032 (raw 01) and 0.034 (raw 02) on the same library and decoys (`matched_fdp02.txt`).

**Bound caveat.** FDP_c is the *combined* entrapment estimate, which is an estimated upper bound (Wen et al., Nat Methods 2025). With r ≈ 1 it is about twice the *lower bound* N_e/N. The paired bound cannot be used here: entrapment peptides are not shuffles of this library's human peptides, and only 4.2% share a human precursor's composition. Lower bounds at q ≤ 0.01:

| run | LDA | MLP | OSW export |
|---|---|---|---|
| raw 01 | 5/993 = 0.005 | 7/1222 = 0.006 | 23/1462 = 0.016 |
| raw 02 | 13/1115 = 0.012 | 12/1231 = 0.010 | 25/1491 = 0.017 |

At 1% our scorer is therefore *inconclusive* in Wen et al.'s sense: the upper bound is above the line and the lower bound is at or below it. OSW's export is invalid even by the lower bound. DIA-NN 2.3.1, with its own different library, reports FDP_c 0.007 at 1% on all three runs (lower bound 0.0036). So valid control is achievable on these raw files. The direct evidence of a real excess comes from section 2, not from FDP_c > 1% alone.

## 2. Where it comes from (`calib_tail.txt`)

Before competition, let S_e be the share of entrapment target winners above a score cut, and S_d the same share for decoy winners. Valid target-decoy competition needs S_e/S_d ≤ about 1. The 90% Poisson band is in brackets.

| run, scorer | cut at q ≤ 0.01 | cut at q ≤ 0.02 | cut at q ≤ 0.05 |
|---|---|---|---|
| 01 LDA | 1.22 [0.32, 2.21] | 1.12 | 1.12 [0.83, 1.40] |
| 01 MLP | 1.14 | 1.67 [1.10, 2.25] | 1.52 [1.21, 1.83] |
| 02 LDA | 2.55 [1.38, 3.75] | 1.60 | 1.55 [1.21, 1.89] |
| 02 MLP | 2.14 [1.12, 3.19] | 1.65 | 1.41 [1.10, 1.72] |

**In the score tail, false targets beat decoys by about 1.1–1.6×, and by up to 2.5× near the 1% cut on raw 02.** The LDA shows it too, so this is mainly not MLP overfitting. It is a mismatch between the decoy model and false targets. The MLP adds to it at the 2–5% cuts on raw 01.

## 3. What it is not

- **Not a single feature.** Over all winners, every subscore has an entrapment-vs-decoy AUC within 0.48–0.53 (`calib_auc_raw0*.tsv`). Among the winners above the 5% cut, no subscore differs at p < 0.02 (`calib_hits.txt`, n = 41–85). The closest is MS1 shape on raw 02 (p 0.02–0.08), and MS1 adds nothing to N. So dropping a feature is not a supported fix.
- **Not true peptides in disguise.** No library entrapment sequence equals a human one under I = L or Q/K, N/D swaps (`calib_il.txt`). Only 1 of the tail hits on raw 02 was reported by DIA-NN at 1%.
- **Borrowed signal: not shown on raw 01, weakly suggested on raw 02.** Here neighbour overlap means the number of top-6 fragments shared with a human precursor within ±1 Th and ±30 s (`calib_borrow.txt`). On raw 02 the MLP's 1% entrapment hits share ≥ 2 fragments with a co-eluting human precursor in 2 of 12 cases, and the 5% hits in 6 of 56. The decoys above the same cuts show 0 of 11 and 2 of 78. Across random precursors the rates are 6.0% (entrapment) and 4.9% (decoys). On raw 01 there is no enrichment: 0 of 5 and 0 of 41 for the LDA. With these counts this is a hypothesis, not a finding.

## 4. Mitigation test (predeclared, `calib_mitig.txt`)

**Pass rule:** q ≤ 0.01 FDP_c ≤ 0.015 on **both** runs, with N at FDP_c 1% within 3% of the single-seed MLP or better.

| variant | raw 01: q ≤ 0.01 N (FDP_c), N at FDP_c 1% | raw 02: q ≤ 0.01 N (FDP_c), N at FDP_c 1% |
|---|---|---|
| A LDA | 993 (0.010), 1011 | 1115 (0.023), 885 |
| B MLP seed 0 | 1222 (0.012), 1221 | 1231 (0.020), 1040 |
| C MLP 5-seed ensemble | 1171 (0.010), **1234** | 1251 (0.024), **1136** |
| D rank average of A and C | 1138 (0.012), 1099 | 1208 (0.027), 1009 |

**The test fails.** The ensemble is the better ranker: N at matched true error goes up on both runs (1234 against 1221, and 1136 against 1040), and it calibrates raw 01 better at 2–5%. It does not fix raw 02 at 1% (FDP_c 0.024). On raw 02 the ensemble still reaches only 88% of OSW+Percolator's 1287.

## 5. Structural checks (`calib_struct.txt`, `calib_entr*.txt`, `calib_lowmass.txt`, `calib_y1.txt`, `calib_scan.txt`)

- **Groups match on precursor m/z, charge, length, library RT, oxidation, C content, fragment type and ordinal, XIC points and RT span.**
- **An out-of-scan-range artefact in the library.** No signal exists below 150.06 Th in these XICs, which is the acquisition's lower limit. Yet the library's top-6 fragments include ions below 150 Th:
  - in 9.7% of decoys, against 5.7% of entrapment and 5.2% of human precursors;
  - the decoy set even includes transitions at 55–100 Th, below the 100 Th assay limit;
  - the y1 of a K-terminal peptide is 147.11 Th. Decoy generation swaps the terminal K and R (`switchKR`, the default). An R-terminal target's y1 at 175.12 Th is intense and often in the top 6 (21,905 slots), and its decoy inherits a dead 147.11 Th y1 in the same slot. The reverse case (11,396 K-target slots at 147.11 Th) is rarer, so the dead slots fall unevenly on the decoys.

  A decoy with a dead slot never reaches the score tail (0 of 23,446 above the 5% cut), which explains why decoys show more empty traces. **It is a real defect, but it does not explain the excess.** Dropping every pair with a slot below 150 Th changes S_e/S_d at the 1%/5% cuts only slightly: LDA raw 01 1.22/1.12 → 1.22/1.07, LDA raw 02 2.55/1.55 → 2.82/1.56, ensemble raw 02 2.67/1.22 → 2.66/1.20. FDP_c at q ≤ 0.01 is unchanged.
  - **Fix for the native pipeline:** pick the top-N fragments only within the acquired scan range, and apply the same rule to targets and decoys. This also returns a dead slot to 5.2% of human targets, a sensitivity gain.
- **Entrapment traces are less often empty than both human and decoy traces at the same fragment m/z** (for example, 300–500 Th: entrapment 0.069, human 0.099, decoy 0.091). Re-weighting entrapment winners to the decoys' strata (empty-trace count × candidate count) lowers S_e/S_d by only about 5% (`calib_entr2.txt`). Entrapment is therefore slightly more signal-prone than human false targets, which biases FDP_c a little upward, but not enough to explain the excess.
- **Decoy library intensities are not realistic** (`calib_libreal.txt`). Among peptides with an internal proline, the top fragment is a y ion starting at P for 47% of human and entrapment precursors, but for only 9.6% of decoys. Decoys inherit the target's intensities ordinal by ordinal on the permuted sequence, so they lose the proline effect. This is a real decoy-realism defect.

  It is **not shown to carry the excess**. Removing the library-intensity features, and then also the library-weighted xcorr, still leaves S_e/S_d at the 5% cut at 1.15–1.22 (raw 01) and 1.22–1.23 (raw 02), against 1.12 and 1.55 with all features. The q ≤ 0.05 FDP_c stays at 0.060–0.062. On raw 02 the 5% excess shrinks, from 1.55 to 1.23. On raw 01 it does not, and the counts at 1–2% are too small to tell.
  - **Fix for the native pipeline:** predict decoy fragment intensities with the same model used for targets, rather than inheriting them by ordinal.

## 6. Decoy-family regeneration test (in progress; scripts in `../decoyregen/`)

**Setup.** The library's targets are kept, and three decoy families are generated from them with `OpenSwathDecoyGenerator`. All three are merged into one PQP, so a single extraction serves every family:

- **O**: pseudo-reverse with `switchKR`. This is the control: it reproduces the library's decoys exactly (100% sequence overlap).
- **P**: pseudo-reverse without `switchKR`, so no terminal K/R swap and no precursor mass shift.
- **S**: shuffle without `switchKR`.

OSW is used for XIC extraction only, locally, with `onDisc` reading and a 26 GB memory cap, at 34 ppm MS2 and 55 ppm MS1 with a ±100 s RT window. Our own scorer is then run with each family in turn as the decoy set (`score2.py`).

**Pass rule, set before running:** a family fixes calibration if, on both runs, its q ≤ 0.01 FDP_c is ≤ 0.015 and its S_e/S_d at the 5% cut lies within the 90% band of 1. At the same time, family O on the same extraction must reproduce the excess; otherwise the test is inconclusive.

**First attempt: invalid, not used** (`score02_offcentre.txt`). OSW's own auto-iRT gave a linear fit with rsq 0.78 (the default gate is 0.95). The RT windows came out off-centre: the median |DIA-NN RT − window centre| was 50 s, against 19 s for the reference XICs, and only 77% of confirmed peaks fell inside, against 95%. As a result, even the control family O found only 243 IDs at FDP_c 1%, against 885 on the reference XICs. Calibration comparisons from that extraction are not meaningful. For the record, family P gave 737 (LDA) and 870 (MLP) against the control's 243 and 330, which is an unexplained difference in itself.

**Second attempt.** Each run gets its own library copy, whose library RT is replaced by the run-calibrated RT used by the reference analysis (`calib_lib.py`; this RT is a monotone function of library RT). OSW calibration is switched off. One remaining difference from the reference: no m/z correction is applied.

### Raw 02 results (raw 01 pending)

**Control reproduces.** With family O, the MLP gives 1055 IDs at FDP_c 1%, against 1040 to 1119 on the reference XICs. The LDA gives 842, against 885. At q ≤ 0.01 the decoy-count estimate is again anti-conservative: FDP_c is 0.019 for the MLP and 0.029 for the LDA, and the S_e/S_d band at the 5% cut lies above 1 (1.46 [1.15, 1.77] for the MLP).

**The predeclared scorer is undefined for family S.** Under `score2.py` the shuffle family gives 0 IDs. Its iteration-0 ranking puts no target winner at q ≤ 0.15, so the semi-supervised loop has no positives. A post-hoc fallback seeds with the top 200 target winners when fewer pass the cut. It is applied to every family alike (`score2_seed.py`), and O and P barely change under it.

| Family, raw 02 (seed fallback) | q ≤ 0.01 N (FDP_c) | q ≤ 0.05 FDP_c | S_e/S_d at the 5% cut | N at FDP_c 1% |
|---|---|---|---|---|
| O LDA | 1031 (0.029) | 0.104 | 2.11 [1.71, 2.52] | 842 |
| O MLP | 1137 (0.019) | 0.073 | 1.46 [1.15, 1.77] | 1055 |
| P LDA | 933 (0.017) | 0.099 | 2.01 [1.62, 2.41] | 856 |
| P MLP | 1008 (0.010) | 0.062 | 1.23 [0.93, 1.53] | 1035 |
| S LDA | 912 (0.018) | 0.077 | 1.49 [1.15, 1.84] | 865 |
| S MLP | 797 (0.003) | 0.038 | 0.81 [0.57, 1.05] | 1130 |

**Family S calibrates partly by borrowing its own target's signal.** Shuffled decoys keep the target's precursor m/z (97.8% within 0.01 Th). Many also keep fragments: 9.5% share at least 3 of their top-6 fragments (within 34 ppm) with their own target, against 1.3% for P and 0.1% for O (`family_overlap.txt`). In the decoy tail this is enriched. Among S decoy winners above the q ≤ 0.05 cut, 36% share at least 3 fragments, against 9.7% of all S decoy winners (`s_tail02.txt`). A decoy that reads its target's own XICs is not a model of a false target. It also lowers q for the wrong reason: it raises the decoy score where the target is real.

**Clean-pair set (`score3.py`, reading declared before running).** All pairs whose S decoy shares at least 2 of its top-6 fragments with its target are dropped, for every family. That removes 20% of candidates and leaves the same targets in all three. Each family is then scored with the MLP at 3 seeds:

| Family, clean, raw 02 | MLP q ≤ 0.01 FDP_c (seeds 0/1/2) | q ≤ 0.05 FDP_c | S_e/S_d at the 5% cut |
|---|---|---|---|
| O | 0.008 / 0.011 / 0.010 | 0.063 / 0.066 / 0.065 | 1.29 / 1.34 / 1.30, every band reaching just below 1 (0.92 to 0.98) |
| P | 0.005 / 0.009 / 0.013 | 0.046 / 0.049 / 0.045 | 0.92 / 0.99 / 0.90 |
| S | 0.005 / 0.003 / 0.010 | 0.040 / 0.044 / 0.041 | 0.82 / 0.88 / 0.83 |

On the clean set, the P and S families meet the per-run rule (q ≤ 0.01 FDP_c ≤ 0.015, and the 5% band includes 1). The control condition fails on this subset, though: O's 5% bands all reach just below 1. The O–P/S contrast at q ≤ 0.05 is consistent across seeds (FDP_c 0.063–0.066 against 0.040–0.049). By the declared reading, this run is inconclusive.

**Where O's excess sits (`score4.py` (a), declared before running).** On the full set, at the q ≤ 0.05 cut, O's S_e/S_d is split by stratum. In the 20% of pairs that drop out of the clean set it is 2.44, 1.87 and 2.41 (bands lie above 1). In the remaining 80% it is 1.19, 1.51 and 1.10. The excess is therefore concentrated in that stratum (2 of 3 seeds clearly). Those targets are short: median length 10 against 18, 34% at 8 residues or fewer against 5%, and 55% at charge 2 against 18%. Entrapment share is equal (0.50; `strat_char.txt`). Short peptides have few distinct b/y ions, so any rearrangement collides with them. For O this means its pseudo-reverse decoys (no fragment overlap) are easier to tell from these targets than false targets are.

**MS1 borrowing is not the explanation (`score4.py` (b)).** P and S decoys sit at the target's precursor m/z, so their MS1 subscores read the target's precursor trace. Without the MS1 group on the clean set, P and S still sit below O at q ≤ 0.05: median FDP_c 0.043 (P) and 0.041 (S), against 0.060 (O). That is not within 0.01 of O, so P/S calibration is not attributed to MS1 borrowing on this run.

**Interim reading (raw 02 only; superseded by the raw 01 section below, where P, the short-peptide concentration and the MS1 result do not replicate):**
- Decoys without the precursor mass shift (P, and S on clean pairs) track false targets better than the library's `switchKR` decoys.
- The library decoys' excess is concentrated in short, low-diversity peptides.
- Shuffling borrows target signal on short peptides and needs an overlap filter. OpenSwathDecoyGenerator has a similarity check against the target, and the native generator should enforce one.
- This is not yet a confirmed fix. It needs raw 01, and the control fails on the clean subset.
- **Confound for P:** without `switchKR`, the pseudo-reverse keeps the C-terminal residue, so a P decoy shares its target's y1 and sits at the same precursor m/z. P therefore also borrows some target signal (mean 0.30 shared top-6 fragments, against 0.15 for O). The clean-pair filter is defined by S overlap, not P's. P calibrating better may thus partly be the same borrowing that S shows, not better modelling of false targets. This tension with section 4 (the mass shift sign/size showed no dose response) is unresolved.

### Raw 01 replication and verdict

The raw 01 outputs are `score01_seed.txt`, `s_tail01.txt`, `score3_01.txt` and `score4_01.txt`.

**Control:** O MLP gives 1113 IDs at FDP_c 1%, against 1222 at q ≤ 0.01 on the reference XICs. Its q ≤ 0.01 FDP_c is 0.011, and its q ≤ 0.05 FDP_c is 0.081, with the 5% band at 1.66 [1.32, 2.01]. The excess shows at 5% and not at 1%, as on the reference XICs.

| Clean pairs, MLP, seeds 0/1/2 | raw 01 q ≤ 0.01 FDP_c | raw 01 q ≤ 0.05 FDP_c | raw 01 5% S_e/S_d | raw 02 q ≤ 0.05 FDP_c |
|---|---|---|---|---|
| O | 0.009 / 0.000 / 0.012 | 0.078 / 0.074 / 0.062 | 1.64 / 1.51 / 1.26 (2 of 3 bands above 1) | 0.063 / 0.066 / 0.065 |
| P | 0.013 / 0.012 / 0.015 | 0.086 / 0.086 / 0.077 | 1.73 / 1.76 / 1.59 (all above 1) | 0.046 / 0.049 / 0.045 |
| S | 0.010 / 0.011 / 0.014 | 0.056 / 0.051 / 0.064 | 1.15 / 1.02 / 1.28 (all include 1) | 0.040 / 0.044 / 0.041 |
| S without MS1 | 0.023 / 0.022 / 0.013 | 0.063 / 0.063 / 0.067 | – | 0.039 / 0.041 / 0.050 |
| O without MS1 | 0.002 / 0.007 / 0.010 | 0.073 / 0.070 / 0.078 | – | 0.064 / 0.060 / 0.055 |

**What does not replicate:**
- **P (pseudo-reverse without `switchKR`)** is worse than the control on raw 01. Removing the precursor mass shift is not a fix, which agrees with section 4.
- **The short-peptide stratum.** On raw 01, O's 5% excess is as large in the clean 80% of pairs (1.73 / 1.41 / 1.41) as in the short-peptide 20% (1.48 / 1.84 / 1.39). The raw 02 concentration was not robust.
- **"MS1 borrowing is not the explanation."** On raw 01, removing the MS1 group moves S's median q ≤ 0.05 FDP_c to 0.063, within 0.01 of O's 0.073. By the declared reading, S's advantage on raw 01 is attributed to MS1 borrowing: S decoys sit at their target's precursor m/z and read its MS1 trace. Without MS1, S also loses q ≤ 0.01 validity (0.022–0.023).

**Verdict against the predeclared rule.**
- **S on clean pairs** meets the per-run calibration rule on both runs. But the control condition holds only on raw 01 (on raw 02 O's bands reach just below 1), so the formal result is inconclusive.
- **S's improvement is confounded** by signal borrowing: fragment overlap on the full set, and the MS1 trace on raw 01.
- **No decoy family tested is shown to fix the calibration.**
- **Without any borrowing channel** (O on clean pairs, and S without MS1), every family sits at about the same excess: q ≤ 0.05 FDP_c 0.055–0.078 over both runs. This points away from the decoy sequence rule as the cause. Decoys that share no signal with any real peptide under-represent false targets in the tail, whichever rearrangement makes them. The likely candidate is that false targets (entrapment) pick up chimeric signal from co-eluting real peptides more often than decoys do. That is untested.

**For the native pipeline:**
- A decoy at the target's precursor m/z must not read the target's MS1 trace. Use a shifted decoy MS1 m/z, or keep MS1 features out of the decoy comparison.
- Shuffled decoys need a fragment-overlap filter against their own target.
- Neither change is expected to remove the excess on its own.

**Limits:**
- Decoy intensities are inherited.
- There is no scan-range fragment filter and no m/z correction.
- Extraction is OSW's.
- The seed fallback and the clean-pair, stratum and MS1 analyses are post-hoc (each reading was declared before its own run).
- S_e/S_d counts in the tail are small (tens of entrapment winners).

## 7. Co-eluting real signal (chimera) test: rejected (`../decoyregen/chimera.py`)

**Hypothesis.** False targets (entrapment) pick up fragment signal from co-eluting real peptides more often than decoys do.

**Measure.** For every winner (family O, full set, MLP seeds 0/1/2), k is the number of its top-6 library fragments within ±17 ppm of a theoretical b/y ion of a co-eluting real precursor. The real precursors are the DIA-NN Q ≤ 0.01 IDs of that run: about 82k, used as a diagnostic map only. Co-eluting means the same 2 Th isolation window, with the winner's apex inside that precursor's RT.Start to RT.Stop ± 5 s. The decoy winners are then reweighted to the entrapment k distribution, and S_e/S_d is recomputed at the q ≤ 0.05 cut. The reading was declared before running: the hypothesis is supported if the reweighting removes at least half the excess, with its band including 1, for ≥ 2 of 3 seeds on both runs.

| | k ≥ 1 share, entrapment / decoy | 5% S_e/S_d pooled → reweighted (seeds 0/1/2) | excess removed |
|---|---|---|---|
| raw 01 | 0.356–0.369 / 0.319–0.331 | 1.66→1.60, 1.50→1.45, 1.41→1.36 | 0.10–0.11 |
| raw 02 | 0.364–0.370 / 0.327–0.334 | 1.46→1.40, 1.60→1.54, 1.41→1.35 | 0.10–0.15 |

**Result.** A small mechanism exists: entrapment winners overlap co-eluting real peptides 3.7 points more often than decoys, consistently. But it explains only 10–15% of the excess, and every reweighted band stays above 1. On raw 02 the excess is largest in the k = 0 stratum (1.75–1.93), the winners with no overlap at all. **Rejected on both runs.**

**Limits.**
- Only DIA-NN-identified peptides are mapped. Unidentified co-eluting signal and non-b/y ions are not.
- The window grid is approximate.
- Fragments come from theory, not observed spectra.

### Where this leaves the diagnosis

Rejected or minor (each at most about 15% of the excess):
- single features;
- I/L;
- the precursor mass shift and decoy sequence rule (section 6);
- borrowed signal;
- dead low-m/z slots;
- entrapment stratum mismatch;
- library-intensity features;
- chimeric co-elution with identified peptides.

What remains:
- The excess (decoy-count estimate about 1.4–1.7× too optimistic in the 5% tail) is robust over the decoy family, the ranker (LDA and MLP) and both runs.
- DIA-NN, on the same raws and with its own library and decoys, reaches FDP_c 0.007 at its 1% cut.

Untested candidates:
1. **The library itself.** Entrapment and human targets share the same predicted RT/intensity model. Decoys inherit their target's RT and ordinal intensities, though, so their norm_rt and library-match features follow the target's prediction quality, not that of a random sequence. The library-feature ablation (section 5) moved the excess only on raw 02.
2. **Entrapment that is not false.** S_e/S_d is computed on precursor winners before any target-decoy competition, so the pairing in the TDC cannot explain it. If a few entrapment peptides are actually present in the sample, or match real peptides, they would inflate the small entrapment tail counts (tens) directly. Examples are sequences identical to human peptides under a modification or isobaric substitution, and contaminants. Test: do the tail entrapment winners coincide with DIA-NN entrapment IDs at the same RT more often than the bulk do?

Candidate 2 is tested in section 8.

## 8. Entrapment truth, feature-space separability and fragment commonness (`../decoyregen/`)

**Entrapment that is not false: rejected (`ent_real.py`).** DIA-NN reports 297 entrapment IDs per run at Q ≤ 0.01, a rate of 0.0036. Of our 55–64 tail entrapment winners, DIA-NN also reports at most 1 (0–1.8%). Removing those moves the 5% S_e/S_d by 0.00–0.03, and every band stays above 1 on both runs.

**Classifier two-sample test (`c2st.py`, diagnostic, uses entrapment labels).** A gradient-boosted classifier tells entrapment winners from decoy winners on the 20 subscores with out-of-fold AUC **0.556 on both runs**. The features carrying it are log_width, lib_corr, log_sn, intensity and xcorr_coel_w. Entrapment winners are wider (median log_width 2.43 against 2.29) and have slightly higher library correlation. Reweighting decoy winners by the out-of-fold density ratio brings the 5% S_e/S_d from 1.46 to 1.10 [0.87, 1.34] on raw 02, and from 1.66 to 1.28 [1.02, 1.54] on raw 01. Within winners above the q ≤ 0.20 cut, though, the AUC is 0.50–0.52. So most of the difference is a weak, bulk shift in what false targets look like, not something visible only in the tail. Entrapment winners look like broader, more library-like background than decoy winners do.

**Fragment commonness (`commonness.py`, reading declared before running).** Commonness is the mean log count of human-library top-6 fragments (same 2 Th window, ±17 ppm) matching a precursor's top-6 fragments. Entrapment top-6 fragments are more common than decoy fragments: 0.513 against 0.438, AUC 0.562 on both runs. Ion-type mix and the share below 300 Th are identical (0.118). Reweighting decoy winners to the entrapment commonness distribution removes 15–29% (raw 02) and 26–37% (raw 01) of the 5% excess. Only one band out of six reaches 1. The declared rule needed half, so it is **not supported as the main cause, but it is the largest single label-free factor found so far.** The reason for the difference is not resolved here. Real tryptic peptides share sub-sequences, and so b/y ions, with other real peptides more than reversed sequences do; this fits the result but is not tested.

**Label-free mitigation, commonness-weighted TDC (`com_tdc.py`, rule declared before running): FAIL.** Each decoy is counted with weight P_target(bin)/P_decoy(bin) over the library's commonness distribution. This uses no run data and no entrapment labels: library mean 0.792 for targets against 0.668 for decoys, with weights 0.59–1.54. Results:
- q ≤ 0.05 FDP_c: 0.068 / 0.075 / 0.067 (raw 02) and 0.066 / 0.070 / 0.060 (raw 01), against 0.070–0.081 unweighted.
- q ≤ 0.01 FDP_c: 0.019 / 0.006 / 0.010 and 0.007 / 0.008 / 0.013.

The rule needed q ≤ 0.05 FDP_c ≤ 0.06 for ≥ 2 of 3 seeds on both runs. It moves in the right direction, with a cost of 0–8% of IDs at q ≤ 0.01.

**Where this leaves it.** The anti-conservativeness is a weak, broad difference between false targets and decoys. False targets look more like real peptides: more common fragments, wider and more library-like background peaks. No single library or run property tested explains more than about a third of it, and the decoy sequence rule does not explain it. A label-free correction based on fragment commonness recovers about a quarter. For the native pipeline, the honest options are:
1. Decoys that match target fragment commonness: generate candidate decoys and select or weight them by commonness against the target library. This is the direct form of the weighting above.
2. A calibration margin in the reported q, validated by entrapment on design runs.
3. Accept the decoy-count estimate as anti-conservative at 5% and report entrapment FDP_c beside it.

None of these is implemented.

## 9. Label-free bulk-null density-ratio TDC: PASS on both design runs (`../decoyregen/bulknull*.py`)

**Idea.** Section 8 found that false targets differ from decoys by a weak, broad shift, not one confined to the tail. Low-scoring target winners are almost all false (DIA-NN-identified precursors are 0.3% of them), so that shift can be learnt without entrapment labels.

**Method.**
1. A gradient-boosted classifier is trained on target winners against decoy winners. It uses only the bulk: winners scoring below the median decoy-winner score. It is cross-fitted by pair fold on the 20 subscores.
2. Each decoy winner gets the weight w = p/(1−p) · n_d/n_t. This is clipped to [0.2, 5], applied out-of-fold, and normalised to mean 1 over the bulk decoys.
3. The paired TDC counts decoys by weight: q = (Σ w over kept decoys + 1) / kept targets.

The ranking is unchanged, so the N at each entrapment FDP_c is unchanged. Only q moves. Bulk AUC is 0.56–0.57. The mean decoy weight in the tail is 1.2–1.5, so decoys that look like false targets count more.

**Exploratory (`bulknull.py`, rule declared before running).** With single MLP seeds 0/1/2, the 5% S_e/S_d goes from 1.41–1.66 to 1.02–1.26, and every band includes 1. The rule needed q ≤ 0.01 FDP_c ≤ 0.015 and q ≤ 0.05 FDP_c ≤ 0.06. That held for 2 of 3 seeds on each run (seed 2 on raw 02 and seed 1 on raw 01 fail at 5%: 0.067 and 0.069).

**Confirmation (`bulknull_confirm.py`, fixed configuration and rule declared before running).** Score: the 5-seed MLP ensemble. Primary: median threshold, boosting, clip [0.2, 5].

| | plain TDC q ≤ 0.01: N (FDP_c, lower) | q ≤ 0.05 FDP_c | **bulk-null q ≤ 0.01: N (FDP_c, lower)** | **q ≤ 0.05 FDP_c** |
|---|---|---|---|---|
| raw 01 | 1143 (0.016, 0.008) | 0.067 | **1038 (0.006, 0.003)** | **0.050** |
| raw 02 | 1112 (0.013, 0.006) | 0.082 | **996 (0.004, 0.002)** | **0.055** |

The rule passes on both runs. Entrapment-matched N at FDP_c 1% is 1101 (raw 01) and 1083 (raw 02); the ranking is untouched. At q ≤ 0.01 the weighted estimate is now conservative (FDP_c 0.004–0.006), and it reports 9–10% fewer IDs than plain TDC. Most of that is IDs that plain TDC reported above its nominal rate.

**Sensitivity (no rule).** Bulk threshold at the 30th or 70th percentile, or clip [0.1, 10]: q ≤ 0.05 FDP_c 0.049–0.058, all pass. **Logistic-regression weights do not correct** (0.064 and 0.083), so the shift is nonlinear in the subscores.

**Contamination (`bulknull_contam.py`, reading declared before running).** True targets in the bulk could teach the classifier "looks like a real peptide" and over-correct. Dropping every DIA-NN-identified target from the bulk training set changes q ≤ 0.05 FDP_c by +0.006 (raw 02, to 0.061) and +0.002 (raw 01). The decoy-weight Spearman correlation is 0.94. So contamination is not driving it. On raw 02, the cleaner version lands just above the 0.06 bar.

**Status and limits.**
- This is a design-run result on two runs (raw 03 is held). It is the first label-free mitigation to pass a predeclared rule, but it is not a guarantee.
- The weights estimate a density ratio from about 180k bulk winners per class, at AUC 0.57. They correct a weak shift, and if the shift changes with score the correction is extrapolated into the tail. The p30/p70 sensitivity is reassuring, but that remains the assumption.
- The 1% statements rest on small entrapment counts. The lower bound at q ≤ 0.01 is 0.002–0.003 on both runs.
- It uses OSW extraction and the OSW-like subscores. Before the native scorer adopts it, it must be re-validated on native features.

**For the native pipeline.** After the cross-fitted scorer, fit the bulk-null classifier on the same features (a boosted-tree or small-MLP head) and report the weighted q. Keep plain TDC q as a secondary column. This adds one model fit per run and needs no new data.

## Reading

1. The anti-conservativeness is a property of **decoy winners under-representing false targets in the score tail**, by about 1.1–1.6×. The pipeline's ranker only adds to it. OSW+Percolator, with the same decoys, is worse: invalid even by the lower bound. No ranker change alone has fixed it. At 1% our scorer sits in the inconclusive zone between the entrapment bounds (section 1). The cause is still unidentified. These have been checked and are not it: single features, I/L, the precursor mass shift, borrowed signal (weak), dead low-m/z slots, entrapment stratum mismatch (each about 5% at most), and library-intensity features (raw 02 partly, raw 01 not).
2. **The decoy precursor mass shift is not supported as the cause.** The library's decoys are pseudo-reversed with a precursor mass shift of ±12, ±28 or ±44 Da (99.9% of pairs; checked in the PQP). The entrapment peptides carry no shift; they are scored at their own m/z. If the shift caused the excess, tail survival should vary with the sign or size of the shift (`calib_shift.py`). It does not vary consistently:
   - for negative against positive shifts, S_e/S_d at the 5% cut is 1.13 against 1.10 (raw 01 LDA) and 1.35 against 1.83 (raw 02 LDA);
   - the |12| and |44| groups hold only 0–7 decoys above the cut.

   The hypothesis is not excluded, but there is no in-data support. A clean test would need regenerated decoys at the target's own m/z.
3. **Levers, in order of how directly they hit the cause:**
   - (a) A **decoy-generation test**: regenerate decoys without the mass shift and with model-predicted fragment intensities, and separately as shuffles instead of pseudo-reverse. In both, restrict fragments to the scan range. Then extract and rescore, and repeat the tail test. This needs a new library and an extraction run on raws 01 and 02: OSW used 1 h 28 min and 37.6 GB of RAM per 3-run job on the cluster, and OpenDIATraces is capped at 16,384 samples. Native OpenDIA extraction is the natural vehicle, since OSW's XICs are fixed to these decoys. None of the in-data tests found a single cause, so the expected effect is uncertain.
   - (b) **More decoys per target** (for example two decoy families, with q computed from their mean). This lowers the variance of the tail estimate but does not remove a bias.
   - (c) Report **a decoy-count estimate corrected by a measured tail factor** only as a diagnostic, never in production. Entrapment-based calibration would not be available on real data.
   - (d) Keep the **5-seed ensemble** as the ranker. It is the best N at matched true error on both runs.
4. The pass rules stay as they are: q ≤ 0.01 must hold at FDP_c ≤ 0.015, and on raw 02 nothing so far meets that.
