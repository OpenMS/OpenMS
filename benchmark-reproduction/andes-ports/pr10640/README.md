# PR #10640 benchmark: knockoff q-values for variant PSMs on PXD022787 (donor genotypes)

Question: does the calibration andes (bigbio/andes #173/#175/#180) reports for knockoff variant q-values show up
in the OpenMS port (`VariantDatabase`, `FalseDiscoveryRate -PSM_variant_database`)?

## Data and ground truth

- PXD022787 (EN-TEx, four donors, one TMT10 set, Orbitrap Fusion, ion-trap CID MS2), the dataset andes used.
  51 of the 69 runs (`runs_used.txt`; 1.83 M MS2 spectra). The 18 runs `Fusion_Entex10plex_P16` … `RT40-0X`
  were not downloaded: the downloader hit the session's 2 h background limit.
- Variants: UniProt `humsavar.txt` (current release, downloaded 2026-10-19), single-residue substitutions of reviewed human proteins whose
  reference residue matches the Swiss-Prot sequence, I↔L swaps left out: 83,473 variants (andes: 83,472 single-residue
  variants of UniProt). `scripts/prepare_variants.py`.
- Truth: the personal protein FASTA deposited with the dataset (`EntexPersonalProteins_Contaminates_Isobaric.fasta`)
  tags every sequence with the donor haplotypes that carry it (`1H1` … `4H2`, `A` = all, `R` alone = reference
  carried by nobody). A variant PSM is *carried* if its peptide, extended to the substitution site when the site lies
  next to it (new or lost cleavage site), occurs (I = L) in a carried sequence of the same gene; *absent* if instead
  the same stretch with the reference residue occurs there; otherwise *unlabelled*. True FDP = absent / (absent +
  carried), on PSMs and on distinct variants. `scripts/evaluate.py`. No knockoff PSM was labelled carried (the null
  is clean). Tier arm: as in andes, variants a donor carries get tier `sample`, all others `curated`.

## Builds

- PR head e7cc3807 (as benchmarked) and fc78916f (fix, see below); baseline = merge-base 2ef86937.
- Class-FDR competitor: `origin/claude/andes-class-fdr` 349328c1 (#10636), `FalseDiscoveryRate -PSM_class variant_class`
  with a meta value set by `scripts/annotate_variant_class.py` (variant = every protein a variant entry).
- Ubuntu 24.04 system packages, `-DWITH_THERMO_RAW=ON` (raw → mzML with FileConverter, MS2 only).

## Pipeline

1. `VariantDatabase -in human_sprot.fasta -variants variants_tiers.tsv -out var_ko4_tiers.fasta -knockoffs 4 -enzyme Trypsin/P`
   (83,473 variant + 65,044 knockoff entries; 53,814 variants without a knockoff position) and `-knockoffs 0` for the
   no-knockoff database (ko0). `var_ko4_onetier.fasta` = the same file without the `tier=` tokens.
2. Search DB: the variant/knockoff entries are full-length protein copies (247 M residues, 22× Swiss-Prot; TTN
   alone 47.6 M), which ProSE could not search in 15 GB. `scripts/trim_variant_entries.py` cuts each entry to a
   window around its substitution that begins and ends at Trypsin/P sites, ≥ 4 sites and ≥ 45 residues from the site
   (accessions unchanged; every ≤ 2-missed-cleavage peptide holding or bordering the site is kept): 27.8 M residues.
   Reversed `DECOY_` decoys (DecoyDatabase).
3. ProSE per run (`scripts/search_params.sh`, `scripts/process.sh`): frozen #10364 Lumos CID TMT parameters
   (Trypsin/P, 2 MC, 7–40, z 2–5, ±20 ppm, isotope −2..1, 0.5 Da, CAM fixed, TMT6plex N-term fixed, Ox(M)),
   except that K carries TMT6plex **or** Carbamyl as variable modifications (≤ 2 per peptide): the samples are
   heavily carbamylated (the deposited PD search has more Carbamyl-K than TMT-K), and with fixed TMT-K almost nothing
   passes 1 %. Top 1 hit, `-Search:rescore false`, `-Search:database:chunk_size 100000`.
4. `scripts/rescore_fdr.sh`: PercolatorAdapter (in-process) over all 51 runs pooled, `-score_type svm
   -post_processing_tdc -subset_max_train 500000`, seeds 1, 42, 137; PeptideIndexer `-IL_equivalent`; PSMs with
   Percolator q > 0.15 dropped (`scripts/filter_percolator_q.py`; a score threshold, so TDC and knockoff q-values up to
   0.10 are unchanged — verified identical on seed 1); then per cutoff 0.01/0.05/0.10:
   - `off`: `FalseDiscoveryRate -PSM true -protein false -FDR:PSM <q>` (global target-decoy q),
   - `ko_onetier` / `ko_tiers`: + `-PSM_variant_database var_ko4_onetier.fasta` / `var_ko4_tiers.fasta` (PR e7cc3807),
   - `*_fix`: the same with fc78916f,
   - `class` (ko0 database): #10636 `-PSM_class variant_class`,
   - baseline FalseDiscoveryRate on the same input as `off` at 0.01 (md5 in `results/baseline_md5.txt`).

## Results (`results/summary.tsv`, mean of 3 seeds; per seed in the `_per_seed` columns)

| DB | Arm | q | Variant PSMs | carried / absent | True FDP (PSMs) | Variants carried / absent (FDP) | Non-variant PSMs |
|---|---|---|---|---|---|---|---|
| ko0 | global TDC | 0.01 | 4,167 | 1,586 / 1,997 | 55.7% | 460 / 601 (56.6%) | 274,256 |
| ko0 | class FDR (#10636) | 0.01 | 2,834 | 1,328 / 1,129 | 45.9% | 391 / 306 (43.9%) | 271,930 |
| ko0 | global TDC | 0.05 | 6,413 | 1,889 / 3,687 | 66.1% | 628 / 1,613 (72.0%) | 310,396 |
| ko0 | class FDR (#10636) | 0.05 | 3,781 | 1,521 / 1,729 | 53.2% | 442 / 484 (52.3%) | 307,192 |
| ko4 | global TDC (option off) | 0.01 | 4,033 | 1,573 / 1,893 | 54.6% | 457 / 565 (55.3%) | 272,351 |
| ko4 | knockoffs, one tier (fix) | 0.01 | 80 | 75 / 1 | 1.3% | 26 / 1 (3.7%) | 272,333 |
| ko4 | knockoffs, tiers (fix) | 0.01 | 583 | 561 / 0 | 0.0% | 184 / 0 (0%) | 272,333 |
| ko4 | global TDC (option off) | 0.05 | 6,100 | 1,839 / 3,445 | 65.2% | 601 / 1,483 (71.2%) | 306,724 |
| ko4 | knockoffs, one tier | 0.05 | 182 | 170 / 7 | 3.9% | 61 / 6 (9.0%) | 306,697 |
| ko4 | knockoffs, one tier (fix) | 0.05 | 178 | 169 / 4 | 2.3% | 60 / 3 (4.8%) | 306,697 |
| ko4 | knockoffs, tiers | 0.05 | 683 | 650 / 4 | 0.6% | 215 / 4 (1.8%) | 306,697 |
| ko4 | knockoffs, tiers (fix) | 0.05 | 679 | 649 / 1 | 0.15% | 214 / 1 (0.5%) | 306,697 |
| ko4 | global TDC (option off) | 0.10 | 8,219 | 2,138 / 5,086 | 70.4% | 803 / 2,618 (76.5%) | 329,536 |
| ko4 | knockoffs, one tier | 0.10 | 400 | 359 / 30 | 7.8% | 134 / 25 (15.8%) | 329,509 |
| ko4 | knockoffs, one tier (fix) | 0.10 | 390 | 358 / 22 | 5.7% | 132 / 17 (11.6%) | 329,509 |
| ko4 | knockoffs, tiers | 0.10 | 778 | 731 / 10 | 1.3% | 260 / 9 (3.2%) | 329,509 |
| ko4 | knockoffs, tiers (fix) | 0.10 | 768 | 730 / 1 | 0.14% | 258 / 1 (0.4%) | 329,509 |

At 1 % the fix changes nothing (rows identical to e7cc3807). Per-run counts: `results/per_file_seed1_q0.05.tsv`.
Knockoff log lines per stratum: `results/knockoff_strata_logs.txt`.

andes (README, same dataset, 69 runs, its own engine): no knockoffs, class q ≤ 0.01 / 0.05: 65 % / 71 % FDP;
knockoffs at 0.05: 201 PSMs, 4 / 160 absent, **2.5 %** (variants 3.4 %); at 0.10: 307 PSMs, 14 / 233, **6.0 %**
(8.2 %); prior tiers at 0.05: one tier 5.6 % (218 / 13), tiers split 2.4 % (329 / 8).

Other checks: option off is byte-identical to the baseline build for all 3 seeds (`results/baseline_md5.txt`).
Knockoff entries in the searched database cost 0.7 % of the non-variant PSMs at 1 % (274,256 → 272,351; andes:
−0.30 % TMT, −0.43 % Astral) and 29 % search time (ProSE 49.8 → 64.0 min for 51 runs); FalseDiscoveryRate with the
option is within run-to-run noise of without it (269–296 s vs 276–300 s on 1.7 M PSMs, seed 1).

## Reproduce

`scripts/fetch_one.sh` (PRIDE download), `scripts/convert.sh` + `scripts/process.sh` (convert, search),
`scripts/rescore_fdr.sh` (rescoring, FDR arms, evaluation), `scripts/summarize.py results/eval results/summary.tsv`.
Database checksums: `db_sha256.txt`. Paths in the scripts are those of the benchmark machine.
