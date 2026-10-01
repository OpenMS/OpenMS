# Prototype fixes for #10378 and #10379: 20-file benchmark (2026-10-01)

Supplement to OpenMS issue [#10364](https://github.com/OpenMS/OpenMS/issues/10364). Follow-up to the
[untangled-PR benchmark](../pr-untangled-20261001/README.md). Each prototype adds switches to the head of one PR,
so every variant can be run as its own arm on that PR's code. Each prototype's defaults reproduce its PR head, which
the identity checks below confirm.

The changes are prototypes, not pushed to the PRs. Several of them follow what ANDES (`bigbio/andes` `b7eaece`) does
in the same places:
- ANDES models fragment charges 1–3 on ion-trap data and only charge 1 on deisotoped data.
- ANDES narrows only the precursor window, to median ± max(2, 3σ + 0.5) ppm.
- ANDES fits no fragment kernel.

## Prototypes

| Build | Base | Commit (patch in this folder) | New parameters (default reproduces the PR head) |
| --- | --- | --- | --- |
| `proto10379` | #10379 head `16ff0344` | `377d3c7a` (`proto-10379-calibration-switches.patch`) | `calibration:apply` = `both` / `precursor` / `fragment` / `none`; `calibration:precursor_window` = `quantile` / `robust` (median ± max(2 ppm, 3·1.4826·MAD + 0.5 ppm)); `calibration:subset` = `top_tic` / `stride`; `scoring:mass_error_kernel_fit` = `full` / `shift` / `none` |
| `proto10378` | #10378 head `241f9a70` | `129e4d4c` (`proto-10378-ion-prior-context.patch`) | `annotate:ion_prior_fragment_charges` = `1` / `auto` (`auto`: charge 1 if deisotoped, otherwise up to min(z − 1, 3)); `annotate:ion_prior_residue_context` = `false` / `true` (adds a level for proline C-terminal / D,E N-terminal of the cleavage below the position context) |

Both builds pass `HyperScore_test`, `FragmentIndex_test` and `ProSEAlgorithm_test`, and for `proto10378` also
`FragmentIonLikelihoodModel_test`, which gains sections for the residue context. `sources.json` gives the arms with
their parameter deltas and files, and the build hashes.

**Identity checks** (`checks.tsv`, 17 of 17 pass):
- `proto10379` with default parameters gives byte-identical PIN and native TSV to the #10379 head on all 14
  high-resolution files.
- `proto10378` with default parameters does the same against the #10378 head on two files.
- On deisotoped HF-X data, `fragment_charges=auto` gives a PIN identical to charge 1.

## Protocol

Unchanged from the untangled-PR benchmark:
- Frozen `native_dedup:<dataset>:default` parameters plus the listed delta.
- Percolator 3.09.0 `-Y -U`, seeds 1, 42 and 137; target PSMs at q ≤ 0.01.
- Mean over seeds within a file, then over files within a group.
- Native counts use rank-one TDC, (D+1)/T ≤ 0.01.
- The develop and PR-head arms are the runs from the untangled-PR benchmark.
- #10379 arms run on the 14 ppm (high-resolution) files.
- Fragment-charge arms run on the 6 low-resolution files, plus `hfx_A2` as the deisotoped control.

A container restart interrupted 2 of the 14 `m_shift_frag` searches. Their incomplete folders were deleted and both
searches were rerun from the start (`run_proto2.log`). 161 searches, 0 failures.

Seed-to-seed relative SD within a file is 0.1–2.1% (Lumos LFQ 0.1–0.2%, Velos/Astral/timsTOF up to 2%).

## Results

### A. #10379 calibration windows, HyperScore
Cells are % vs develop / native TDC delta vs develop.

| Arm | HF-X HCD | Astral HCD | Lumos HCD LFQ | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- |
| develop (PSMs) | 4380.3 | 2189.6 | 5496.1 | 2244.9 | 1015.0 |
| #10379 head: `auto` = precursor + fragment windows | -3.34% / -4 | +5.09% / +331 | -4.16% / -142 | +0.61% / +174 | +0.79% / -4 |
| precursor window only | -2.83% / -9 | -3.46% / -5 | -3.71% / -148 | -0.87% / -2 | +0.79% / -4 |
| **fragment window only** | **-0.91% / +27** | **+5.48% / +349** | **-0.43% / +15** | **+1.94% / +206** | **+0.00% / +0** |
| robust precursor window only (ANDES-style) | -4.73% / -68 | -5.15% / +159 | -7.24% / -302 | -5.97% / -67 | -0.57% / -6 |
| robust precursor window, stride sample | -6.40% / -119 | -4.44% / +167 | -7.09% / -287 | -7.34% / -69 | -3.12% / -26 |

**Fragment window only:**
- Per file, it beats the #10379 head on 11 of 14 files (exceptions: `astral_A2` −0.6%, timsTOF −0.2/−1.4%), by up to +7.5% on `lumos_lfq_5192`.
- It keeps the whole Astral gain and adds +1.9% on TMTpro.

**Precursor narrowing loses in every form:**
- Most clearly on the robust window, which is much narrower (Lumos LFQ 3.4–4.0 ppm, Astral 4.4–4.7 ppm, HF-X and TMTpro 7.7–9.2 ppm).
- On Astral the PR's quantile window barely moves (`astral_A2`: [−19.8, +18.2] ppm). The native TDC count stays
  unchanged, yet Percolator loses 3.2–3.7% on each file.
- A likely reason (not tested): Percolator learns from the wrong candidates at the window edge.

### B. Mass-accuracy scorer
Cells are % vs develop (in parentheses: vs the scorer with the fixed 7 ppm kernel and no calibration) / native TDC
delta vs develop.

| Arm | HF-X HCD | Astral HCD | Lumos HCD LFQ | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- |
| scorer, fixed 7 ppm kernel, no calibration | +0.34% / -13 | +4.13% / +417 | -0.14% / +37 | -0.07% / +231 | +0.62% / +10 |
| #10379 head: fitted kernel + both windows | -3.04% (-3.37) / -601 | +3.10% (-0.98) / -186 | -3.39% (-3.26) / -269 | +1.22% (+1.29) / +72 | +1.25% (+0.62) / +48 |
| fitted kernel (center + width), windows unchanged | +0.25% (-0.09) / -714 | +0.42% (-3.56) / -204 | -0.17% (-0.03) / -194 | -0.07% (-0.00) / +52 | -0.72% (-1.34) / +52 |
| fitted center only, windows unchanged | +0.09% (-0.25) / +10 | +4.95% (+0.79) / +383 | -0.08% (+0.06) / +35 | -0.45% (-0.38) / +227 | -0.57% (-1.19) / +61 |
| **fitted center only + fragment window** | **-0.85% (-1.19) / +12** | **+6.52% (+2.30) / +404** | **-0.44% (-0.30) / +18** | **+2.22% (+2.29) / +232** | **-0.57% (-1.19) / +61** |
| fitted center only + robust precursor window, stride | -5.39% (-5.71) / -138 | -0.97% (-4.90) / +412 | -7.15% (-7.02) / -312 | -8.10% (-8.03) / +86 | -3.96% (-4.55) / +22 |

**Fitted width:**
- Against the fixed 7 ppm kernel, it lowers the native TDC count by 701 (HF-X), 621 (Astral) and 231 (Lumos LFQ)
  per file on average. The fitted SDs are 1.8–3.5 ppm on the Orbitrap and Astral files.
- It costs 3.6% Percolator yield on Astral (`astral_B3` −9.4%, fitted SD 1.76 ppm).

**Fitted center only:** about neutral against the fixed kernel (−1.2 to +0.8% per group). Native counts are within −34 to +51 of it.

**Best arm in this table:** the fixed-width kernel with the fitted center, plus the fragment window. It reaches +6.5%
on Astral and +2.2% on TMTpro, and −0.4% to −0.9% on HF-X, Lumos LFQ and timsTOF. Compared with HyperScore and the
fragment window (table A), the scorer adds about +1.0% on Astral and +0.3% on TMTpro, and −0.6% on timsTOF. The two
arms' fragment windows differ slightly, because each calibration pass scores with its arm's scorer.

### C. #10378 ion priors
Cells are % vs develop (in parentheses: vs the #10378 head). Native scores are unchanged by design.

| Arm | Velos CID | HF-X HCD | Astral HCD | Lumos HCD LFQ | Lumos CID TMT | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- | --- | --- |
| #10378 head (fragment charge 1) | -0.37% | +0.30% | +4.60% | -0.11% | +0.92% | +1.76% | +2.35% |
| **fragment charges `auto`** | **+6.20% (+6.59)** | identical (1 file, deisotoped) | – | – | **+3.14% (+2.19)** | – | – |
| cleavage-residue context | +0.54% (+0.91) | +0.63% (+0.33) | +4.82% (+0.20) | -0.04% (+0.07) | +1.92% (+0.99) | +1.44% (-0.32) | +1.72% (-0.61) |
| both | +6.16% (+6.55) | – | – | – | +4.03% (+3.08) | – | – |

**Fragment charges `auto`:**
- Gains on every low-resolution file: Velos +5.4/+6.5/+8.0%, Lumos CID TMT +1.6/+2.2/+2.8% against the PR head.
- It turns Velos CID, the PR's weakest group, into its largest gain.

**Residue context:** within seed noise overall (−0.6 to +1.0% per group). Adding it on top of `auto` changes −0.7
to +3.0% per file (mean +0.4%).

## Files

| File | Contents |
| --- | --- |
| `proto-10379-calibration-switches.patch`, `proto-10378-ion-prior-context.patch` | The prototype commits (`git am` onto the PR heads) |
| `sources.json`, `checks.tsv` | Builds, arms with files and parameter deltas, identity checks |
| `proto_per_group.tsv`, `proto_per_file.tsv`, `tables_proto.md` | Results against develop and against the PR heads, per group and per file, with resolved windows and kernels |
| `pr-prototypes-reproduction.zip` (+ `.sha256`) | Scripts (`run.py` arms, `evaluate_proto.py`, `tables_proto.py`, pipelines, harness), build records, logs and per-job `summary.json` / parameters / search and Percolator logs for all 161 searches |

The inputs, nightly and Percolator are those of the 2026-09-30 reproduction bundle referenced in the untangled-PR
supplement.

## Final PR commits, per-file robustness and entrapment (added 2026-10-01)

The prototypes were applied to the PRs: #10378 `b02b538` (fragment charges from the deisotoping decision, no new
parameter) and #10379 `eb74991` (`auto` applies the fragment window only; `scoring:mass_error_kernel_fit=shift`
by default). Each final build reproduces its prototype arm byte for byte (PIN and native TSV):
- `b02b538` against `p_zauto` on all six low-resolution files, and against the #10378 head on `hfx_A2` and
  `astral_A2`.
- `eb74991` against `h_frag` (default) and `m_shift_frag` (mass-accuracy) on all 14 high-resolution files, and
  against develop on two low-resolution files.

**Per-file robustness** (`robustness.py`, `robustness.tsv`):
- Each file's mean over seeds is compared with develop.
- z is the difference over the seed-to-seed standard error. A file counts as up at z ≥ 2 and down at z ≤ −2.
- Seed variance does not capture data-sampling variance; on Lumos LFQ it is tiny, so small changes reach large z.

| Candidate | Files up / flat / down | Where it gains | Where it loses |
| --- | --- | --- | --- |
| #10377 (cutoff) | 0 / 20 / 0 | intensity-scaled input only (see the untangled-PR supplement) | – |
| #10378 `b02b538`, priors on | 10 / 10 / 0 | Velos +4.1 to +7.8%, Lumos CID TMT +2.8 to +3.5%, Astral +2.1 to +6.7% | – |
| #10379 `eb74991` default | 2 / 15 / 3 | Astral +2.4 to +8.8%, TMTpro +1.1 to +3.8% | HF-X −0.7 to −1.3%, Lumos LFQ −0.2 to −0.6% (every file negative) |
| #10379 `eb74991` + mass-accuracy | 5 / 6 / 3 (of 14) | Astral +3.1 to +8.7%, TMTpro +1.6 to +2.8% | HF-X −0.6 to −1.2%, Lumos LFQ −0.7% on 2 files |
| #10335 merged | 4 / 16 / 0 | Velos +8.9 to +16.3%, Lumos CID TMT +1.5 to +5.4% | – |

**Entrapment on the Velos UPS1/yeast files** (`entrapment_velos.py`, `entrapment_velos.tsv`):
- **Design:** the samples (PXD001819) contain yeast and the 48 human UPS1 proteins, but the database holds all
  20.5k reviewed human proteins. A PSM whose peptide (I = L) occurs only in non-UPS1 human proteins is false.
- **Estimator:** the combined FDP estimate is N_E·(1 + 1/r)/N, with r = 3.38 entrapment-only per sample peptide of
  the tryptic space.
- **Scope:** sums over 3 files × 3 seeds at Percolator q ≤ 0.01.

| Arm | PSMs | Entrapment PSMs | FDP (combined) | Extra PSMs vs develop | FDP of the extra PSMs |
| --- | --- | --- | --- | --- | --- |
| develop | 23403 | 172 | 0.95% | – | – |
| #10378 `b02b538`, priors on | 24854 | 185 | 0.97% | +1451 | 1.2% |
| #10335 merged | 26136 | 231 | 1.15% | +2733 | 2.8% |

**What the entrapment check shows:**
- **#10378:** the Velos gain keeps the nominal 1% FDR.
- **#10335:** most of its Velos gain is real, but its FDP rises to 1.15% because the extra PSMs are less clean.
- **Coverage:** the other groups have no entrapment design (HYE samples contain all three species; TMTpro and
  plasma are human only).

## #10378 cross-fitted (added 2026-10-01, later)

Review of #10378 found that the ion-prior model, trained on a run's confident PSMs, also scored those PSMs. The PR now
cross-fits the model (`523214b`):
- Spectra fall into 3 folds by scan index.
- Each fold is scored by a model whose training PSMs are selected, gated and fitted on the other two folds only.

`ef0dcee` was searched on all 20 files (`f78x_priors` in `robustness.tsv`). `523214b` gives an identical PIN and
native TSV on `velos_125_R1` and `astral_A2`; every fold model has 503 or more training PSMs, against a gate of 100.

**Per-group results** (% vs. develop):

| | Velos | HF-X | Astral | Lumos LFQ | Lumos CID TMT | TMTpro | timsTOF |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `b02b538` (not cross-fitted) | +6.20% | +0.30% | +4.60% | -0.11% | +3.14% | +1.76% | +2.35% |
| `523214b` (cross-fitted) | +5.51% | +0.05% | +3.54% | -0.19% | +3.05% | +0.76% | +0.66% |

- **Per file:** 8 up, 11 flat, 1 down (`lumos_lfq_5192`, −0.34%).
- **Velos entrapment:** combined FDP 0.91% (develop 0.95%, `b02b538` 0.97%). The 1290 extra PSMs bring one extra
  entrapment hit.

## Why the gains differ between instruments (added 2026-10-01, later)

Folder `pattern/`. The scripts run from the harness directory of the untangled-PR benchmark with `PYTHONPATH=.`
(they import `evaluate.py`); the ones that read spectra or idXML also need the pyOpenMS nightly of the reproduction
bundle. Effects are % Percolator PSMs at q ≤ 0.01 against develop, mean over seeds 1/42/137, then over files.

| File | Contents |
| --- | --- |
| `spectra.py` → `spectra.tsv` | Raw MS2 properties of the 8,000 benchmark spectra per file: peak count, peaks per 100 Da, share of peaks and intensity a top-20-per-100-Da filter keeps, intensity below m/z 200 |
| `filter_loss.py` → `filter_loss.tsv` | Peaks and intensity at m/z ≥ 200 removed by a top-20 / 40 / 100 per 100 Da filter, and the share of spectra with a 100 Da bin above 20 peaks |
| `fragment_tails.py` → `fragment_tails.tsv` | Signed ppm errors of the matched ions in develop's search, for accepted targets and decoys: robust SD, SD of weak (< 5% of base peak) and strong (≥ 20%) peaks, and the share beyond the calibrated fragment window of #10379 (high-resolution files) |
| `covariates.py` → `covariates.tsv` | Per file: ID rate, near-threshold share N(q ≤ 0.05) / N(q ≤ 0.01) − 1, share of precursors with charge ≥ 3, and the effects of #10378 (`523214b`), #10335, the #10379 fragment window and the mass-accuracy scorer; prints Spearman correlations |
| `classes.py` → `classes.tsv` | The above by analyzer class and by group, with ANDES (auto model) against develop |
| `peak_filter.py` → `peak_filter.tsv` | New arms `wt40` and `wt100`: the develop build with `peaks:window_top` 40 or 100 instead of 20, on the 14 high-resolution files; seeds, z and native TDC change |
| `peak-filter-results.zip` (+ `.sha256`) | `run.py` with the arms, the run log and, for the 28 searches, `summary.json`, parameters, and search and Percolator logs |

**Findings:**
- **Near-threshold PSMs.** Across files, every PR's gain rises with the near-threshold share (Spearman ρ: #10378
  +0.62, #10379 fragment window +0.83, mass-accuracy scorer +0.52, #10335 +0.41). High-ID-rate Lumos LFQ and HF-X
  gain nothing from any of them.
- **Fragment charges (low resolution).** Ion-trap spectra are not deisotoped and 39–52% of their precursors have
  charge ≥ 3. #10335 (+7.4%) and #10378 after its fragment-charge fix (+4.3%) gain there; #10378 was flat on Velos
  before the fix.
- **Fragment window.** The calibrated window helps where weak peaks stay accurate: Astral 2.8 ppm, Exploris TMTpro
  3.5 ppm (+5.5% / +2.0%), against HF-X 5.8 ppm and Lumos LFQ 6.0 ppm (−0.9% / −0.4%). On timsTOF (5.6 ppm) the window
  never narrows below the 20 ppm cap.
- **Peak filter.** The top-20-per-100-Da filter removes 41% of the fragment intensity on Astral, 23% on the ion
  traps, 13% on timsTOF and 4–8% on the Orbitraps. Keeping 100 peaks per 100 Da gains Astral +8.8% (3 of 3 files up)
  and timsTOF +4.4% (2 of 2); the Orbitraps are flat (HF-X +1.2%, Lumos LFQ +0.2%, TMTpro −0.6%). Native TDC
  counts fall on most files, so the gain comes through Percolator's features. Ion-trap files were not run.
- **ANDES.** It selects one of its trained models by instrument class (LowRes, Q Exactive, Orbitrap Astral, TOF,
  timsTOF). Its high-resolution models deconvolute, its low-resolution models model fragment charges 1–3, and it
  applies the top-20 filter only to TMT/iTRAQ ("regresses high-res Astral ~14%"). No TOF model is bundled, so timsTOF
  falls back to `cid_lowres_tryp`; that is the only group where ANDES trails develop (−5.1%; +1.3% with
  `hcd_qexactive_tryp` forced). Its largest lead, +41% on Astral, comes with a dedicated `hcd_astral_tryp` model.

## Dense-spectrum peak quota (added 2026-10-02)

Folder `peak-filter/`. Prototype `20b6e75` on branch `claude/quirky-galileo-tkpwe4-prose-peak-filter` (off develop
`ea3f2c1`; `peakfilter-dense-quota.patch`). Where the local filter keeps full quotas (high resolution), a spectrum in
which `peaks:window_top` (20) peaks per 100 Da window would remove more than `peaks:dense_intensity_loss` (0.2) of its
deisotoped intensity keeps `peaks:dense_window_top` (100) peaks per window. It was benchmarked as build
`peakfilter_bench`, the same patch on `f5ea2d04`: develop's Boost change (#10223) adds library symbols the pinned nightly
lacks. All three class tests pass, and the TOPP replays differ from develop only by the two recorded parameters.

**20-file suite** (`pf_eval.py` → `pf_per_file.tsv`, `pf_eval.log`; z ≥ 2 counts as up):

| Arm | Files up / flat / down | Astral | timsTOF | HF-X | Lumos LFQ | TMTpro | Ion trap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `peaks:window_top` 40, all spectra | 7 / 10 / 3 | +7.0% | +4.4% | +0.6% | +0.2% | +0.9% | −1.1 to −2.8% |
| `peaks:window_top` 100, all spectra | 7 / 11 / 2 | +8.8% | +4.4% | +1.2% | +0.2% | −0.6% | −0.5 to −2.3% |
| Dense quota, loss 0.1 | 4 / 10 / 0 | +8.2% | +2.2% | +0.6% | +0.3% | 0.0% | not run (unchanged) |
| **Dense quota, loss 0.2** | 3 / 17 / 0 | +9.0% | +2.1% | +0.2% | −0.15% | −0.7% | identical |
| Dense quota, loss 0.3 | 2 / 12 / 0 | +7.0% | +1.3% | identical | −0.06% | +0.4% | not run (unchanged) |

- **Spectra counted dense** at loss 0.2: Astral 96%, timsTOF 26%, Orbitrap 1–5% (`filter_loss_processed.tsv` has the
  per-spectrum distribution).
- **Identity:** with `peaks:dense_window_top` 0 (4 files), and on all ion-trap files, PIN and native TSV are
  byte-identical to develop.
- **Native TDC:** falls 1–2% on Astral in every variant.

**Doubled search spaces** (`entrapment_E.tsv`, `entrapment_F.tsv`; Percolator q ≤ 0.01, 3 seeds; combined FDP 2·N_E/N):
- `_E`: targets plus one shuffled entrapment protein per target (`make_entrapment.py`). K, R and P stay fixed and the
  other residues of each segment are shuffled, so every target peptide has a same-mass, same-composition twin (r = 1).
- `_F`: the same size, but with fully shuffled proteins of unrelated peptide masses, like a foreign proteome
  (`make_foreign.py`).
- Database hashes are in `entrapment-databases.sha256`.

| Arm | Astral (A2 / B1 / B3) | timsTOF (30 / 50 min) | Astral combined FDP |
| --- | --- | --- | --- |
| Dense quota vs develop, benchmark database | +18.3 / +3.3 / +5.4% | −0.7 / +4.8% | — |
| Dense quota vs develop, `_E` | −5.2 / +2.9 / +2.2% | −1.2 / +2.9% | 1.88% → 1.87% |
| Dense quota vs develop, `_F` | +7.3 / −2.4 / −1.3% | +2.9 / +4.4% | 1.40% → 1.20% |
| #10378 (`523214b`, priors on) vs develop, `_E` | +2.5 / +4.2 / +0.9% | −1.3 / +1.4% | 1.88% → 1.67% |

Neither change inflates the error rate. The dense quota's Astral gain, however, does not survive a change of the search
space, while #10378's does. The dense quota is therefore not proposed.

`peak-filter-results.zip` (+ `.sha256`) holds `run.py` with all arms, the build records, run logs and, for every
search, `summary.json`, parameters, and search and Percolator logs.
