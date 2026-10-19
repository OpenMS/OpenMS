# PR #10630: `scoring:candidate_tie_break intensity` on the frozen #10364 ProSE suite

Benchmark of OpenMS PR [#10630](https://github.com/OpenMS/OpenMS/pull/10630) (port of bigbio/andes #99/#116):
ProSE ranks fragment-index candidates with equal vote counts by the summed intensity of their matched peaks
before the `scoring:max_candidates_per_spectrum` cut. Measured 2026-10-19, one machine (4 cores, Ubuntu 24.04, GCC 13.3).

## Summary

- **Default unchanged:** with the option at `none`, PIN and idXML are identical to the merge base on all 20 files.
- **`intensity` does not raise identifications in OpenMS:** 65,870.4 → 65,730.0 target PSMs over 20 files (−0.21%;
  1 file up, 18 flat, 1 down). Astral −1.2% (astral_B3 z = −3.2), HF-X +0.3%, other groups within ±0.4%. On the
  entrapment databases −0.00% over 17 files. Entrapment FDP moves both ways (Astral 1.40 → 1.53%, TMTpro 1.36 → 1.06%,
  Velos 1.18 → 1.12%), within the noise of 40–110 entrapment PSMs per group.
- **Why:** ProSE has one candidate per (peptide, charge, isotope error), and the tie keys after the vote count were
  |isotope error|, sign, charge, index. The PR puts intensity before the isotope keys. In the spectra whose rank-1
  peptide changes, the new rank-1 has isotope error 0 in 11% of cases instead of 46% (Astral), and about half of
  them are decoys either way. andes has no such multiplicity and its replaced key was the form id (mass order).
- **Exploratory variant** (intensity after the isotope keys; local patch, tuned on this suite): +0.05% over 20 files
  (0/20/0), +0.04% on the entrapment databases. It removes the Astral loss but gains nothing either.
- **andes on the same inputs** (695acf7, flag off → on, entrapment databases): Astral +2.0% (flat per file) at FDP
  1.77 → 1.72%, HF-X −1.0% (2 files down), Lumos LFQ +0.2% (FDP 1.17 → 1.36%). andes reported +1.7% on the full
  Astral file PXD070049. The Astral direction holds in andes, but not in OpenMS.
- **Cost:** +1.3% wall time, +1.0% CPU time, unchanged peak RSS.
- **Second dataset, MHC class I** (id_bruker, unspecific, `results_idb/`): +0.8% (3 files flat, noisy seeds) on the
  target-decoy database, −0.4% on the entrapment database. The FDP of 0.52 → 0.68% rests on 6–7 entrapment PSMs. No
  gain either, although 30–40% of the (per-chunk) queries end inside a vote tie.

## Builds

| key | source | role |
| --- | --- | --- |
| `base` | `74048ef59f` (merge base of the PR with develop) | default-identity check |
| `pr` | PR head `eb878e0440` | main A/B: option `none` (default) vs `intensity` |
| `diag` | PR head + `scripts/tiebreak_diag.patch` (local, not pushed) | counts vote ties at the candidate cap; its output equals `pr` + `intensity` byte for byte |
| `iso` | PR head + `scripts/variant_intensity_after_isotope.patch` (local, not pushed) | exploratory: intensity key **after** the isotope-error keys |
| andes | bigbio/andes `695acf7` (parent of #116, which made the tie-break the default) | like-for-like reference, `--fragment-index-intensity-tiebreak` off/on |

OpenMS: system packages, `-DCMAKE_BUILD_TYPE=Release -DWITH_GUI=OFF`, target `ProSE`; FragmentIndex_test and
ProSEAlgorithm_test pass on the PR head. andes: `cargo build --release --locked --bin andes` (Rust 1.87.0).

## Protocol

- Inputs: the frozen 20 mzML of #10364 (8,000 MS2 each) regenerated with the bundle's `prepare.py`; every selected
  mzML matches `inputs.json` (`selected_sha256`). FASTAs from the HYE ProteoBench FASTA (sha256 `d9ac434d…`) with
  `scripts/make_dbs.py`, a pyOpenMS-free copy of `prepare_databases.py`; the hy/hye/human hashes match `databases.json`.
- Entrapment databases: `scripts/make_entrapment.py` (from `sage-comparison-20261002`): one shuffled entrapment
  protein per target with K/R/P fixed (paired masses, r = 1), reversed decoys of both. hye/human hashes equal the
  earlier supplements (`fb77e3df…`, `76d25071…`). Velos uses natural entrapment (UPS1 in yeast against all human
  proteins of `hy`, r = 3.376).
- Search: TOPP `ProSE` with the frozen `native_dedup:<dataset>:default` parameters (jobs.json) mapped to the current
  parameter names (`precursor:isotopes` and `snes_enabled` no longer exist and are dropped); everything else at the
  ProSE default. `-threads 4`, `-out_pin`, no `-rescore`. `scripts/run_one.py` writes the full command line.
- Rescoring: external Percolator 3.06.0 (OpenMS/THIRDPARTY `fd8c89f`, `Linux/x86_64/Percolator`), `-Y -U --seed s`
  for s = 1, 42, 137, on ProSE's PIN (no PIN normalisation). Target PSMs at q ≤ 0.01; mean over seeds per file, mean
  over files per group. z = Δ / sqrt(sd_a²/3 + sd_b²/3); up/down when |z| ≥ 2.
- Entrapment FDP: combined estimator, shuffled 2·ΣN_E/ΣN (seed means summed over a group's files), natural
  N_E(1+1/r)/N over 3 files × 3 seeds.
- Timing: 3 replicates per file and arm (base, none, int) on an otherwise idle machine, arm order rotated per
  replicate; median wall time, CPU time and peak RSS (`/usr/bin/time -v`).
- andes: frozen andes search space of #10364 (20 ppm, z 2–5, top-n 10, 2 MC, length 7–40, ≤1 Ox(M), isotope −1..2,
  `--precursor-cal off`, `--protocol standard`, the target-decoy FASTA with `--decoy-strategy none`), retrieval left at
  the commit's `auto` (fragment-ion index for high-resolution data), entrapment databases, Percolator as above on
  andes' native PIN. Astral, HF-X and Lumos LFQ only (the files on which andes uses the fragment index and that have
  an unlabelled entrapment database).

## Results (`results/`)

| file | content |
| --- | --- |
| `per_group_int_vs_none.tsv`, `per_file_int_vs_none.tsv` | main A/B, frozen databases |
| `per_group_int_vs_none_entrapment_db.tsv` (+ per file) | same on the entrapment databases |
| `entrapment_velos_natural.tsv`, `entrapment_shuffled.tsv` (+ per file) | FDP per group and arm, incl. andes and the variant |
| `per_group_none_vs_base.tsv`, `identity_base_vs_none.tsv` | default identity: PIN (main runs) and idXML (timing replicates, without `date` and `search_engine_version`) |
| `vote_ties.tsv` | per query: candidates above the cap and whether the cap falls inside a vote tie (diag build) |
| `top1_changes.tsv`, `accepted_changes.tsv`, `mechanism_per_group.tsv`, `rank1_isotope_error.tsv` | where the changes come from |
| `per_group_intiso_vs_none*.tsv` | exploratory variant |
| `per_group_andes_on_vs_off_entrapment_db.tsv` (+ per file) | andes like-for-like |
| `resources.tsv` | timing and memory |
| `run_summaries.json` | per run: input/output hashes, seeds, wall/CPU time, RSS |

## Second dataset: MHC class I immunopeptidomics (id_bruker)

- Data: archive.openms.de `testfiles/id_bruker/classI/high/` (PXD058436, MHCquant2, timsTOF Pro), techRep1 of
  DN02_Spleen, DN13_Corpus and DN12_Liver (`scripts/idb_manifest.json`, with byte counts). 8,000 MS2 per file with
  the frozen `prepare.py` sampling (copied next to the manifest; charge 2–5 only, so singly charged precursors are
  excluded). `scripts/idb_inputs.sh` resumes truncated downloads and checks the byte count.
- Database: `testfiles/human_sp.fasta` (20,431 entries, sha256 `b205ed99…`) plus reversed `DECOY_` decoys
  (`321b4b05…`), and the shuffled entrapment version from `make_entrapment.py` (`93f43ae9…`). For unspecific
  digestion, the shuffle keeps protein composition and length (hence r ≈ 1 in candidate count) but does not pair
  peptide masses.
- Search: the frozen tims parameters with `enzyme 'unspecific cleavage'`, `peptide:enzyme_specificity none`, length
  8–12, no fixed modification, ≤1 Ox(M), and `database:chunk_size 5000` (9 chunks; each chunk applies its own cap of
  50, so `vote_ties.tsv` counts 9 queries per spectrum). Same Percolator protocol. About 60 s and 3 GB per search.
- Only 3 files and about 700–1,300 PSMs per file. Seed-to-seed spread is large (corpus 501–559 vs 479–646), and the
  entrapment estimate rests on single-digit counts.

| arm | PSMs (sum of 3 files) | Δ vs `none` | entrapment DB PSMs | FDP (combined) |
| --- | ---: | ---: | ---: | ---: |
| `none` | 2616.7 | | 2171.3 | 0.52% |
| `intensity` | 2638.0 | +0.81% (0/3/0) | 2163.3 (−0.37%) | 0.68% |
| variant (intensity after isotope keys) | 2632.7 | +0.61% (0/3/0) | 2088.3 (−3.8%) | 0.77% |

The default is identical to the merge base on all three files (PIN). The rank-1 peptide changes in 6–9% of the
spectra. Among the changed rank-1 PSMs, isotope error 0 falls from 48–61% to 18–23%, and about half are decoys,
as on the tryptic suite.

## Reproduce

Paths in the scripts assume the root `/home/user/bench` (`BENCH` in each script).

```bash
# inputs (frozen bundle prose-andes-reproduction.zip unpacked as $BENCH/frozen, sha256 e24cbe74…)
bash scripts/inputs.sh                       # FASTAs + 20 mzML, hash-checked
for db in hye hy human; do python3 scripts/make_entrapment.py frozen/broad_benchmark/data/${db}_td.fasta frozen/broad_benchmark/data/${db}_E_td.fasta; done
git clone --depth 1 --filter=blob:none --sparse https://github.com/OpenMS/THIRDPARTY && git -C THIRDPARTY sparse-checkout set Linux/x86_64/Percolator
# builds: install bin/ProSE, lib/, share/OpenMS of each build under /opt/openms-<key> (base, pr, diag, iso)
bash scripts/queue.sh jobs/jobs_pr.txt; bash scripts/queue.sh jobs/jobs_bd.txt; bash scripts/queue.sh jobs/jobs_iso.txt
while read a d; do python3 scripts/andes_run.py $a $d; done < jobs/jobs_andes.txt
bash scripts/queue.sh jobs/jobs_timing.txt   # on an idle machine
python3 scripts/evaluate.py results; python3 scripts/isotope_check.py <datasets> > results/rank1_isotope_error.tsv
bash scripts/idb_inputs.sh; bash scripts/queue.sh jobs/jobs_idb.txt; python3 scripts/evaluate.py results_idb --idb
```
