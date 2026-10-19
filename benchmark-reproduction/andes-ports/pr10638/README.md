# PR #10638: ProSE `chimeric:enabled` (co-isolated precursors)

Question: does the gain bigbio/andes reports for its chimeric cascade (`--chimeric`) also show up in OpenMS ProSE, and is
the FDR of the extra (`chimeric_rank` 2) PSMs controlled?

The summary is in the PR description. The tables here are the source.

## Builds

| arm | commit | binary |
|---|---|---|
| base | `77cccda7` (merge-base of the PR with develop) | `/opt/openms-base` |
| pr_off | PR head `2c4aef5d`, `chimeric:enabled` false (default) | `/opt/openms-pr` |
| pr_on | PR head `2c4aef5d`, `-Search:chimeric:enabled` (max_coisolated 4, max_kl 0.3, mz_tolerance 10 ppm: defaults) | `/opt/openms-pr` |

Ubuntu 24.04 with system packages, Release, `-DWITH_GUI=OFF -DWITH_THERMO_RAW=OFF`; 4 cores, 15 GB RAM.
`ProSEAlgorithm_test` and `CoIsolationDetector_test` pass at the PR head.

## Data

| file | source | MS2 searched |
|---|---|---|
| HF-X | ProteoBench M2 `LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML` (archive.openms.de), whole run | 109,507 |
| Astral | ProteoBench M8 `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML` (archive.openms.de; the file andes benchmarks), RT slice 623–1048 s (`FileFilter -rt 623:1048`, MS1 and MS2 kept, the central 50% of the MS2 scans) | 61,750 |
| FAIMS | `testfiles/FAIMS/ExplorisFAIMS-DDA-LFQ.mzML` (Exploris 480, FAIMS CV −45/−55/−65, E. coli sample), whole run | 78,206 |

Why a slice for Astral: on the whole run (123,539 MS2) `pr_off` peaks at 11.55 GB RSS and `pr_on` was killed by the
OOM killer at 13.95 GB (15 GB container). See *Memory* below.

Database: the #10364 HYE target-decoy FASTA (`make_td.py`, sha256 `efbf5ac2…`, source FASTA `d9ac434d…`) with one
shuffled mass-matched entrapment protein per target (`frozen_make_entrapment.py`, the frozen
`sage-comparison-20261002/make_entrapment.py`; sha256 `fb77e3df…`, the same file as in pr10629). 63,778 targets
(31,889 + 31,889 `ENTRAP_`), 63,778 reversed `DECOY_` decoys. Input hashes: `results/input_sha256.txt`.

## Search and rescoring

`scripts/run_search.sh`: the #10364 parameters — Trypsin/P, 2 missed cleavages, length 7–40, charge 2–5, fixed
Carbamidomethyl (C), at most one Oxidation (M), precursor ±20 ppm, fragment 20 ppm, `scoring:max_candidates_per_spectrum`
50, `report:top_hits` 10, `peptide:deduplicate` true, calibration off, `-threads 4`; every other parameter at its
default. ProSE writes `-out_pin` with `-Search:rescore false`.

Rescoring (`scripts/perc.sh`): external Percolator 3.06 (OpenMS/THIRDPARTY Linux binary) `-Y -U --seed {1,42,137}`, the
same in all arms. In the PIN a `chimeric_rank` 2 PSM has the SpecId and ScanNr of its spectrum and its own ExpMass (the
co-isolated m/z), so Percolator's target-decoy competition (grouped by scan and ExpMass) treats it as its own precursor.
For the evaluation only, `perc.sh` appends `#<row>` to every SpecId (all arms alike), so that rank 1 and rank 2 can be
told apart in Percolator's output.

The default route (ProSE's in-process Percolator, `-Search:FDR:PSM 0.01`; `scripts/run_inproc.sh`,
`scripts/count_idxml.py`) was run on HF-X as a check.

## Evaluation (`scripts/evaluate.py`)

Per seed: target PSMs at q ≤ 0.01 (entrapment PSMs count as targets), split by `chimeric_rank` from the PIN; unique
peptides (modified sequence); entrapment PSMs (all proteins `ENTRAP_`); combined FDP 2·N_E/N (r = 1); PSMs and peptides
at an entrapment FDP of 1% (largest q-ordered prefix with 2·N_E/N ≤ 0.01, the measure andes reports). Means over the
three seeds. On `pr_off` the counts equal those of the frozen `evaluate_entrapment.py` (pr10629 `frozen_eval.py`; it keys
by scan, so it cannot count two PSMs of one spectrum and does not apply to `pr_on`).

## Results (means over seeds 1, 42, 137; `results/percolator_external.jsonl`)

| file | arm | PSMs q ≤ 1% (rank 1 / rank 2) | peptides | FDP all | FDP rank 1 | FDP rank 2 | PSMs at FDP 1% | peptides at FDP 1% |
|---|---|---|---|---|---|---|---|---|
| HF-X | off | 60,631 | 41,463 | 1.93% | | | 56,593 | 38,751 |
| HF-X | on | 101,052 (59,947 / 41,105) | 47,214 | 1.88% | 1.68% | 2.17% | 93,399 | 43,861 |
| Astral slice | off | 24,618 | 16,690 | 1.96% | | | 22,183 | 15,171 |
| Astral slice | on | 56,917 (23,302 / 33,615) | 20,410 | 1.85% | 1.25% | 2.28% | 49,968 | 18,088 |
| FAIMS | off | 25,533 | 10,088 | 1.09% | | | 25,457 | 10,050 |
| FAIMS | on | 30,366 (25,392 / 4,974) | 10,502 | 1.00% | 0.91% | 1.46% | 30,371 | 10,504 |

Seed SD of the PSM count ≤ 142 in every arm, so every file is up by the #10364 rule (|z| ≥ 2).

- Peptides found only by rank-2 PSMs: 6,196 (HF-X), 4,610 (Astral), 490 (FAIMS).
- Rank-2 PSMs at q ≤ 1% with the same peptide as an accepted rank-1 PSM of the same spectrum (what andes drops before
  Percolator): 33 (HF-X), 82 (Astral), 4 (FAIMS).
- In-process Percolator (HF-X, one run, `results/percolator_inprocess_hfx.txt`): 59,969 → 100,030 PSMs (rank 2: 40,262),
  peptides 40,930 → 46,838, FDP 1.79% → 1.80% (rank 1 1.75%, rank 2 1.88%).

## Default off and unchanged primary search (`results/default_off_sha256.txt`)

- `base` and `pr_off` write byte-identical PINs (HF-X, Astral) and the same idXML apart from the date (HF-X).
- `pr_on` vs `pr_off` (HF-X): the 528,346 rank-1 PIN rows are the same rows; all columns are equal except `ion_prior_llr`
  and `ion_prior_explained`, because the self-trained ion priors are trained on all PSMs including rank 2. After
  Percolator, rank-1 PSMs at 1% change by −1.1% (HF-X), −5.3% (Astral), −0.6% (FAIMS), with a lower rank-1 FDP.

## Runtime and memory (`results/runtime.tsv`)

Median of 3 interleaved runs, 4 threads, search only (no Percolator):

| file | off wall | on wall | off peak RSS | on peak RSS |
|---|---|---|---|---|
| HF-X | 67 s | 118 s | 6.04 GB | 7.47 GB |
| Astral slice | 108 s | 186 s | 8.57 GB | 11.33 GB |
| Astral whole run | 193 s (one run) | killed (OOM) | 11.55 GB | > 13.95 GB |

The index restriction that the option turns off saves nothing here (the spectra reach almost the whole index). The
memory comes from `buildCoIsolatedSpectra_()`: every co-isolated precursor gets a full copy of its spectrum's residual
peaks, plus a copy of the parent's evidence spectrum (ion priors, local fragment evidence); 219,251 co-isolated
precursors on the whole Astral run.

## andes reference

bigbio/andes `fa0195a` README and `docs/benchmarks/README.md`, Astral `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1` (the
same file, whole run), Percolator 3 seeds, 1:1 shuffled entrapment, at an equal entrapment FDP of 1%:
43,076–43,145 → 60,403–60,435 PSMs (+40%, #118) and 66,377–66,424 (+54%, top-N #136); peptides 26,723–26,746 →
29,384–29,400 (+10%; flat before #136). At nominal q ≤ 1% the FDP rose from 1.05% to 1.18%. UPS1 (ion trap): +8.0% PSMs.

andes itself was not run here (building it was not permitted in this environment); its numbers are for the whole run
with its own engine and 10 ppm / N-terminal acetylation settings, so only the relative change is compared.

## Reproduce

```bash
python3 -I scripts/make_td.py ProteoBenchFASTA_MixedSpecies_HYE.fasta hye_td.fasta
python3 -I scripts/frozen_make_entrapment.py hye_td.fasta hye_entrap_td.fasta
FileFilter -in LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML -out astral_A1_rt623-1048.mzML -rt 623:1048
SKIP_PERC=1 scripts/run_search.sh <prefix> pr_off <mzML> hye_entrap_td.fasta <dir>/pr_off
SKIP_PERC=1 scripts/run_search.sh <prefix> pr_on  <mzML> hye_entrap_td.fasta <dir>/pr_on -Search:chimeric:enabled
PERCOLATOR=<THIRDPARTY>/Linux/x86_64/Percolator/percolator scripts/perc.sh <dir>/pr_on   # and pr_off
python3 -I scripts/evaluate.py <dir>/pr_off <dir>/pr_on
scripts/timing.sh                                                    # interleaved timing (paths inside)
```
