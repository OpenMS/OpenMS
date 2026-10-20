# Benchmark of PR #10633: PIP-ECHO `donor_max_q` (port of bigbio/andes #172)

PR #10633 adds `PipEcho:donor_max_q` (ProteomicsLFQ `-PipEcho:donor_max_q`, default 1.0 = off). With a value below 1,
PIP-ECHO transfers a precursor (sequence and charge) to other runs only if its best identification over all runs has a
q-value at or below that value. andes #172 reported that this donor gate (best donor PSM q <= 0.001) cut the transfers
of entrapment precursors on the ProteoBench HF-X data from 728 to 125 at transfer q 0.01.

Question: does the same gate show the same effect in OpenMS (ProteomicsLFQ, MBR with PIP-ECHO)?

## Setup

| | |
|---|---|
| Data | ProteoBench DDA quant ion module, Q Exactive HF-X, 3 + 3 runs (PXD028735), the ProteoBench mzML used by #10585 (`proteobench.cubimed.rub.de/raws/DDA/mzml/`, sha256 in `results/inputs.sha256`; the archive.openms.de copies are a different conversion and not byte-identical) |
| Builds | PR head `786fc7874475b13c93e71207707459ffb52e4e02`; baseline = merge-base `74048ef59f2bccf480c272428a2a43b19d2a25f2`. Both Release builds with system packages (Ubuntu 24.04), configured with `-DGIT_TRACKING=OFF -DOPENMS_GIT_SHORT_SHA1=bench-pr10633`: feature checkpoints are tied to the revision string, and the identical label lets both builds combine the same checkpoints (the PR changes PIP-ECHO linking only, not detection). The identification tools come from the PR-head build (the PR does not touch them). EntrapmentEvaluator (PSM-level check only) and DecoyDatabase come from `claude/andes-entrapment` 2f5b38d1 (#10629); that DecoyDatabase writes the same `HYE_td.fasta` as the PR head's. |
| Identification | #10585 recipe: CometAdapter (Comet 2025.01 rev. 1) 10 ppm, isotope error 0/1, 0.02 Da, Trypsin, 2 missed cleavages, length 7-40, charge 2-5, Carbamidomethyl (C) fixed, Oxidation (M) + Acetyl (Protein N-term) variable (max 3) → PSMFeatureExtractor → PercolatorAdapter (Percolator 3.06, `-score_type pep -post_processing_tdc -subset_max_train 300000 -seed 42`) → Percolator q <= 0.10 → PEP as main score (`scripts/id_one.sh`). Per run 50,906 / 52,121 / 55,935 / 51,758 / 51,605 / 56,258 PSMs at q < 0.01, the numbers of #10585. |
| Databases | ProteoBench arms: ProteoBench HYE FASTA (sha256 d9ac434d…) + `DecoyDatabase -method reverse -enzyme Trypsin` (63410dbe…), as #10585. Entrapment arms: frozen #10364 protocol: HYE target-decoy FASTA (efbf5ac2…) → `make_entrapment.py` (one shuffled entrapment protein per target, K/R/P in place, all 766,991 searchable peptides mass-matched, full-reversal `DECOY_` decoys of targets and entrapment; fb77e3df…), the same file as the #10637 benchmark. |
| Quantification | `ProteomicsLFQ -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation -targeted_only false -pip_echo true` (`scripts/lfq_arm.sh`), PIP-ECHO defaults (`PipEcho:fdr 0.05`). Feature detection once per database into `-feat_dir`; every arm combines the same checkpoints and IDs, only the option under test changes. Then `ConsensusMapNormalizer -algorithm_type rt_local` → `MSstatsConverter`. |
| Metrics | Transfers: PIP-ECHO's accepted transfers (`mbr_transfer_map_indices` / `mbr_transfer_qvalues` of the consensusXML) on consensus features that keep an identification after ProteomicsLFQ's FDR filters, i.e. the reported ones (`scripts/transfer_stats.py`). A transfer is an entrapment transfer when every protein of its peptide is an `ENTRAP_` protein (I = L: a peptide also in a target protein is a target). Transfer FDP = 2·N_E / (N_T + N_E) (r = 1). ProteoBench: official `proteobench` 0.18.6, module `quant_lfq_DDA_ion_QExactive`, format `quantms`, min 3 runs (`scripts/score_proteobench.py`; it reproduces the public quantms submission d01e87b…: 33,040 / 0.2088 / 0.3866 / 0.2236). |

Transfers at q <= 0.01 do not depend on `PipEcho:fdr`. The arms with `PipEcho:fdr 0.01` give the same 11,099 (off) and 7,917 (0.001) transfers as the q <= 0.01 subset of the `PipEcho:fdr 0.05` arms, so both thresholds are read from one run.

## Results

### Option off = baseline

`base_off` (merge-base) and `pr_off` (PR, option at its default) give the same output on both databases. The MSstats CSVs, unnormalized and rt_local, are byte-identical. consensusXML and mzTab are identical up to random unique ids, dates and the parameter list, where the PR adds `PipEcho:donor_max_q = 1.0` (`results/identity_*.tsv`, `scripts/canon_hash.sh`).

### Entrapment (HYE + 1:1 shuffled entrapment)

| arm (`PipEcho:fdr` 0.05) | precursors not transferable | reported transfers q <= 0.01 | entrapment | transfer FDP | reported transfers q <= 0.05 | entrapment | transfer FDP | transfers on unreported precursors (q <= 0.05) |
|---|---|---|---|---|---|---|---|---|
| off (= baseline) | – | 11,152 | 21 | 0.38 % | 27,001 | 48 | 0.36 % | 16,240 |
| `donor_max_q 0.01` | 29,406 of 96,195 (31 %) | 10,396 | 20 | 0.38 % | 26,810 | 45 | 0.34 % | 1,078 |
| `donor_max_q 0.005` | 33,907 (35 %) | 9,533 | 17 | 0.36 % | 24,675 | 43 | 0.35 % | 538 |
| `donor_max_q 0.001` | 42,592 (44 %) | 8,194 | 14 | 0.34 % | 20,524 | 39 | 0.38 % | 118 |
| off, `-proteinFDR 1` | – | 11,263 | 121 | 2.1 % | 27,597 | 476 | 3.4 % | 15,216 |
| `donor_max_q 0.001`, `-proteinFDR 1` | 42,592 (44 %) | 8,207 | 22 | 0.54 % | 20,615 | 66 | 0.64 % | 0 |

Identified precursors (the same in every arm): 151 of 64,034 are entrapment-only with `-proteinFDR 0.01`
(FDP 0.47 %), 1,267 of 66,727 without the protein filter (3.8 %). Entrapment-only precursors quantified in at least
3 runs without the protein filter: 577 of 51,076 (1.13 %) → 490 of 49,325 (0.99 %) with the gate at 0.001. PSM level (EntrapmentEvaluator, all six runs):
combined FDP 0.31 % at Percolator q <= 0.001, 1.6 % at 0.01, 11.6 % at 0.10 (the input filter, i.e. the donor
candidates; `results/entrapment_psm_comet.tsv`).

### ProteoBench (plain HYE database, #10585 setup)

| arm | precursors (min 3) | median \|eps\| | E. coli | CV | rt_local: median \|eps\| | E. coli | CV |
|---|---|---|---|---|---|---|---|
| #10585 PIP-ECHO (reference, nightly before #10585) | 46,375 | 0.2214 | 0.4090 | 0.2356 | 0.1570 | 0.3378 | 0.0950 |
| off (= baseline), `PipEcho:fdr 0.05` | 46,415 | 0.2214 | 0.4088 | 0.2358 | 0.1572 | 0.3353 | 0.0951 |
| `donor_max_q 0.01` | 46,364 | 0.2217 | 0.4102 | 0.2355 | 0.1564 | 0.3405 | 0.0949 |
| `donor_max_q 0.001` | 45,071 | 0.2204 | 0.3996 | 0.2349 | 0.1540 | 0.3271 | 0.0942 |
| off, `PipEcho:fdr 0.01` | 43,573 | 0.2206 | 0.4052 | 0.2327 | 0.1512 | 0.3112 | 0.0904 |
| `donor_max_q 0.001`, `PipEcho:fdr 0.01` | 43,038 | 0.2213 | 0.4061 | 0.2326 | 0.1516 | 0.3119 | 0.0902 |

Gate at 0.001: 41 % of the 95,484 precursors are not transferable. Reported transfers fall from 11,099 to 7,917 at
q <= 0.01 and from 27,602 to 21,081 at q <= 0.05. Transfers on precursors that ProteomicsLFQ does not report fall from
14,055 to 75. ProteoBench metrics on the entrapment database (entrapment rows counted with their source species) are in
`results/summary_entrap.tsv`.

### andes #172 (same data, its own pipeline; `results/andes_pr172_reference.tsv`)

main → donor PSM q <= 0.001: entrapment transfers 728 → 125 (q 0.01) and 1,199 → 229 (q 0.05); precursors (min 3) /
median |eps| / E. coli / CV 67,586 / 0.164 / 0.302 / 0.120 → 64,150 / 0.164 / 0.297 / 0.121 (q 0.01) and
73,203 / 0.177 / 0.471 / 0.134 → 68,197 / 0.174 / 0.439 / 0.133 (q 0.05); about 40 % of the entries not transferred.

| | andes #172 | OpenMS (this benchmark) | verdict |
|---|---|---|---|
| share of precursors gated at 0.001 | ~40 % | 41 % (HYE), 44 % (entrapment DB) | reproduced |
| entrapment transfers, recipe (`-proteinFDR 0.01`) | 728 → 125 / 1,199 → 229 | 21 → 14 / 48 → 39, FDP 0.38 → 0.34 % / 0.36 → 0.38 % | not reproduced: the output has almost no false-donor transfers to remove |
| entrapment transfers without protein FDR (`-proteinFDR 1`) | | 121 → 22 (−82 %) / 476 → 66 (−86 %), FDP 2.1 → 0.54 % / 3.4 → 0.64 % | reproduced (−83 % / −81 % in andes) |
| precursors (min 3) | −5.1 % (q 0.01), −6.8 % (q 0.05) | −1.2 % (`PipEcho:fdr 0.01`), −2.9 % (0.05) | same direction, smaller |
| accuracy | flat at q 0.01, E. coli 0.471 → 0.439 at q 0.05 | flat at 0.01; at 0.05 median \|eps\| 0.157 → 0.154 and E. coli 0.335 → 0.327 (rt_local) | partly (same direction at 0.05, smaller) |

Why the entrapment effect depends on the protein filter: ProteomicsLFQ applies its PSM FDR (pooled, after linking) and the picked protein FDR to the linked map. With `-proteinFDR 0.01`, entrapment-only precursors are mostly removed together with their entrapment proteins, and their transfers go with them. In the off arm, 16,240 accepted transfers sit on precursors that are not reported. The reported transfers then contain only 21 / 48 entrapment transfers, and the gate cannot lower that further measurably (counts this small differ by a few by chance). andes reports every library entry that passes its precursor-level gate, without a protein-level filter. That is what `-proteinFDR 1` approximates, and there the OpenMS gate removes 82–86 % of the entrapment transfers, as in andes. The pipelines also differ in rescoring (Percolator per run here, one model over all runs in andes), in the MBR transfer model (PIP-ECHO SVM vs andes' GBDT with identity decoys) and in the donor candidates (PSMs up to Percolator q 0.10 here, andes' 1 % entry gate), so the absolute numbers are not comparable, only the directions and relative sizes.

Runtime: with the gate the combine step is faster, because PIP-ECHO scores fewer candidates: 26,619 instead of 49,545 on HYE. The four entrapment arms of the first batch ran at the same time on the same machine and took 1,793 s (merge-base, off), 1,795 s (PR, off), 917 s (0.01) and 793 s (0.001) (`results/arm_runtime.txt`).

## Second dataset (Agilent 6550, PXD069532)

Too shallow for PIP-ECHO: 505–719 PSMs per run at q < 0.01, 618 precursors, 6 transfer candidates and no decoy
transfers, so PIP-ECHO keeps no transfer in any arm (its `min_decoys` fallback). The gate leaves 369 (0.01) or 280
(0.001) of the 618 precursors transferable, with nothing to compare. Details in `results/agilent.md`. The SCIEX set of
the same folder was not run: FragPipe finds only 312 MBR values there.

## Reproduce

`scripts/run_all.sh` lists the steps in order (inputs, databases, identification, detection, arms, evaluation).
Detection needs 10–13 GB per run and must run alone on a 15 GB machine. An arm combining checkpoints needs about 1.5 GB.
Detection took about 20 min per run, an arm 13–40 min.
