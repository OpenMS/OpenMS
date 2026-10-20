# PR #10638: ProSE `chimeric:enabled` (co-isolated precursors)

Question: does the gain bigbio/andes reports for its chimeric cascade (`--chimeric`) also show up in OpenMS ProSE, and is
the FDR of the extra (`chimeric_rank` 2) PSMs controlled?

The summary is in the PR description. The tables here are the source. **Round 1** (below the round-2 section) measured
the PR as opened (`2c4aef5d`); **round 2** measures the fixes pushed after it (memory, andes port gaps) and the
ablations that chose the final variant (`f2800c61`).

## Round 2: memory fix, andes port gaps, ablations (`results/round2/`, `scripts/round2/`)

Builds: merge-base `842adef9` (`base2`); `f28b682d` (`pr2`: memory fix + andes semantics); `0e0a9640` (`pr3`: three
andes semantics reverted); `f2800c61` (`pr4`, final: `chimeric:top_hits` 3). The local ablation build is `f28b682d` +
`scripts/round2/ablation.patch` (switches in `PROSE_CHIMERIC_ABLATION`; not part of the PR). Same protocol as round 1
(`scripts/run_search.sh`, `scripts/perc.sh`, `scripts/evaluate.py`), plus the whole Astral run.

What changed in the code (PR branch):
- Memory: co-isolated precursors are searched in rounds (k-th co-isolated precursor of every spectrum), their residual
  spectra scored in batches of 8,192 and kept only as peak masks while the fragment index exists; the residual spectra
  are built after the index is released, and use the parent's peak list for the local fragment evidence (no copies).
- andes semantics kept: sequential peak claiming (the best hit claims its peaks, then each co-isolated winner before
  the next round), no hit with the best hit's peptide, a co-isolated precursor whose winner repeats a PSM of the
  spectrum (peptide and charge) is dropped.
- andes semantics reverted after the ablation: isotope error 0 only, one hit per co-isolated precursor, MS1 tolerance
  max(precursor tolerance, 0.01 Th), claiming only the most intense peak per singly charged ion.
- New `chimeric:top_hits` (default 3): hits per co-isolated precursor (memory).

HF-X ablation (the variant was chosen here; 3 seeds; "@1%" = at an entrapment FDP of 1%):

| arm | PSMs q ≤ 1% (r1 / r2) | peptides | FDP all / r1 / r2 | PSMs @1% | peptides @1% |
|---|---|---|---|---|---|
| off | 60,631 | 41,463 | 1.93% | 56,593 | 38,751 |
| round 1 (`2c4aef5d`) | 101,052 (59,947 / 41,105) | 47,214 | 1.88 / 1.68 / 2.17% | 93,399 | 43,861 |
| andes-faithful (`f28b682d`) | 93,194 (59,959 / 33,235) | 45,886 | 2.02 / 1.64 / 2.70% | 85,553 | 42,218 |
| + top hits for rank 2 | 93,314 | 45,903 | 1.88 / 1.78 / 2.04% | 87,170 | 42,999 |
| + search's isotope errors | 97,838 | 46,627 | 1.92 / 1.56 / 2.48% | 90,356 | 43,183 |
| + all peaks in tolerance | 93,095 | 45,860 | 1.98 / 1.60 / 2.66% | 85,951 | 42,410 |
| + 10 ppm MS1, no floor | 97,068 | 46,403 | 1.99 / 1.62 / 2.58% | 89,206 | 42,797 |
| all four (`0e0a9640`, top hits 10) | 100,971 (60,082 / 40,889) | 47,181 | 1.83 / 1.72 / 1.99% | 94,071 | 44,152 |
| final, `chimeric:top_hits` 1 | 100,169 | 46,974 | 1.88 / 1.53 / 2.39% | 92,641 | 43,579 |
| final, `chimeric:top_hits` 2 | 100,742 | 47,111 | 1.83 / 1.67 / 2.07% | 94,304 | 44,238 |
| **final, `chimeric:top_hits` 3 (`f2800c61`)** | 100,908 (60,111 / 40,797) | 47,173 | 1.83 / 1.72 / 1.99% | 94,217 | 44,226 |

Validation on the other files (not used for the choice):

| file | arm | PSMs q ≤ 1% (r1 / r2) | peptides | FDP all / r1 / r2 | PSMs @1% | peptides @1% | peak RSS |
|---|---|---|---|---|---|---|---|
| Astral whole run (123,539 MS2) | off | 42,469 | 26,349 | 1.89% | 37,808 | 23,810 | 11.55 GB |
| | round 1 | OOM | | | | | > 13.95 GB |
| | andes-faithful | 80,492 (40,870 / 39,622) | 30,045 | 2.42 / 1.47 / 3.40% | 65,098 | 24,605 | 11.77 GB |
| | top hits 10 (`0e0a9640`) | OOM | | | | | > 13.95 GB |
| | **final** | 91,875 (40,821 / 51,054) | 31,703 | 1.90 / 1.48 / 2.23% | 81,296 | 28,233 | 12.95 GB |
| Astral slice (61,750 MS2) | off (round 1) | 24,618 | 16,690 | 1.96% | 22,183 | 15,171 | 8.57 GB |
| | round 1 | 56,917 | 20,410 | 1.85 / 1.25 / 2.28% | 49,968 | 18,088 | 11.33 GB |
| | andes-faithful | 48,272 | 19,327 | 2.61 / 1.29 / 3.86% | 37,052 | 15,058 | |
| | **final** | 55,419 (23,484 / 31,936) | 20,462 | 1.91 / 1.32 / 2.34% | 48,432 | 17,986 | |
| FAIMS (78,206 MS2) | off | 25,533 | 10,088 | 1.09% | 25,457 | 10,050 | 6.07 GB |
| | round 1 | 30,366 | 10,502 | 1.00 / 0.91 / 1.46% | 30,371 | 10,504 | |
| | andes-faithful | 29,504 | 10,433 | 1.04 / 0.98 / 1.47% | 29,468 | 10,413 | 6.89 GB |
| | top hits 10 | 30,338 | 10,502 | 1.00 / 0.92 / 1.42% | 30,341 | 10,505 | 6.90 GB |
| | **final** | 30,331 (25,401 / 4,930) | 10,499 | 0.98 / 0.90 / 1.41% | 30,344 | 10,504 | |

- The ordering holds on the held-out files: final ≥ andes-faithful everywhere (Astral whole run +25% PSMs and +15%
  peptides at FDP 1%); final ≈ round 1 (Astral slice −3% PSMs / −0.6% peptides, FAIMS equal), but round 1 does not fit
  the whole Astral run.
- In-process Percolator (default rescoring, HF-X, final): 100,456 PSMs (rank 2: 40,362), FDP 1.75% (rank 1 1.75%,
  rank 2 1.74%); without the option 59,969, 1.79% (`percolator_inprocess_hfx_final.json`).
- Default off: `base2`, `pr2_off` and `pr4_off` write byte-identical PINs (HF-X, Astral, FAIMS) and the same idXML
  apart from the date (HF-X, FAIMS); `default_off_sha256.txt`. `pr3_on` on HF-X is byte-identical to the ablation arm
  with all four switches.

Runtime and memory of the final build (median of 3 interleaved runs, search only, `runtime_final.tsv`):

| file | wall off → on | peak RSS off → on |
|---|---|---|
| HF-X | 55 s → 90 s (round 1: 67 → 118 s) | 6.04 → 6.75 GB (round 1: 7.47) |
| Astral whole run | 171 s → 252 s | 11.55 → 12.95 GB (round 1: OOM) |

Where the remaining overhead sits (`rss_profile_astral_*.tsv`, RSS sampled once per second): during the scoring, with
the fragment index, the final build is 0.25 GB above the run without the option (11.44 vs 11.19 GB). The peak is now in
the post-processing: 12.6–12.9 GB against 10.7 GB, from 140k more PSM records (1.48 M hits instead of 1.12 M) and their
residual spectra. With 10 hits per co-isolated precursor the post-processing went past 13.95 GB
(`rss_profile_astral_pr3_on_oom.tsv`).

andes and the FDR of rank 2: andes gives secondaries no feature and no competition of their own (its PIN columns
`PrecursorIsotopeKL`/`PrecursorSNR` are always 0); they enter Percolator like any PSM, with the scan and their own
ExpMass. Nothing to port. The rank-2 FDP stays above rank 1 (Astral 2.23% vs 1.48%); class-wise FDR (#10636,
FalseDiscoveryRate `-PSM_class` with `chimeric_rank`) would control it per class.


# Round 1: the PR as opened (`2c4aef5d`)

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
