# PR #10629: EntrapmentEvaluator and DecoyDatabase -entrapment_ratio

The benchmark checks the measuring tool that the other andes-port benchmarks rely on. It covers three questions:

1. Does EntrapmentEvaluator reproduce the entrapment FDP of the frozen #10364 harness on the same results?
2. Are DecoyDatabase's entrapment proteins sound (mass-matched, deterministic, ratio honoured), and do searches against them estimate the FDP like the frozen databases do?
3. How does this compare with the andes method (bigbio/andes `docs/benchmarks`, #48, #112)?

The results tables are in `results/`. The PR description has the summary.

## Builds

| arm | commit | notes |
|---|---|---|
| baseline | `77cccda7` (merge-base with develop) | DecoyDatabase only, for the default-off check |
| PR head | `2f5b38d1` | as opened |
| fix | `f0bddcb3` = `2f5b38d1` + `8371c2b1` (I = L in EntrapmentEvaluation) + `f0bddcb3` (mass-matched entrapment in DecoyDatabase) | pushed to the PR branch |

The build environment was Ubuntu 24.04 with system packages, Release, `-DWITH_GUI=OFF -DWITH_THERMO_RAW=OFF`, 4 cores and 15 GB RAM.

## Data and search (the #10364 protocol)

- **Spectra:** 12 of the 20 frozen files, 8,000 MS2 spectra each, re-sampled with the frozen `broad_benchmark/prepare.py`. Every source file and every sampled mzML matches the sha256 in the bundle's `inputs.json`. The files are Velos PXD001819 (3), HF-X ProteoBench M2 (3), Astral ProteoBench M8 (3) and Exploris TMTpro PXD054559 (3). The frozen protocol has no entrapment design for the Lumos files; timsTOF was skipped for time.
- **Databases:** the frozen `hy`/`hye`/`human` target-decoy FASTAs, regenerated without pyopenms (`prepare_databases_nopyopenms.py`; all three sha256 match `databases.json`; HYE source sha256 `d9ac434d…`).
  - Velos uses natural entrapment: every reviewed human protein other than the 48 UPS1 proteins is entrapment. `make_velos_natural.py` puts `ENTRAP_` in front of those accessions, because EntrapmentEvaluator detects entrapment by prefix.
  - The shuffled r = 1 entrapment DBs come from the frozen `sage-comparison-20261002/make_entrapment.py` (`frozen_make_entrapment.py`; sha256 `fb77e3df…` for hye and `76d25071…` for human, matching `entrapment-databases.sha256`).
- **Search:** ProSE with the frozen `native_dedup:<dataset>:default` parameters from `jobs.json` (`run.py`). The keys `snes_enabled` and `precursor:isotopes` do not exist in this build and are skipped. ProSE runs with its default in-process Percolator (`inproc.idXML`, one seed) and also writes `-out_pin`.
- **External rescoring:** Percolator 3.06 (OpenMS/THIRDPARTY Linux binary) runs on that PIN with `-Y -U --seed {1,42,137}`. The only change to the PIN is the SpecId, rewritten to `<dataset>_<scan>_<rank>` as the frozen `export.py` writes it, so that the frozen scripts can parse it. `target.tsv` then goes through `IDFileConverter -score_type qvalue` to idXML.
- **Determinism:** re-running ProSE and Percolator with the fix build gives a byte-identical PIN and `target.tsv`, and the same idXML apart from the date.

## Reproduce

Layout: put the scripts in `<root>/scripts/`. Unpack the frozen bundle `prose-andes-reproduction.zip` (sha256 `e24cbe74…`) into `<root>/frozen/z/`. Put a sparse clone of OpenMS/THIRDPARTY (`Linux/x86_64/Percolator`) in `<root>/thirdparty/`. Create `<root>/work/broad_benchmark/` with the bundle's `prepare.py` and `manifest.json`, and put `inputs.json` in `<root>/work/`.

```bash
bash scripts/fetch_prepare.sh work/broad_benchmark velos_125_R1 ... astral_B3      # download, check, sample, check
python3 -I scripts/prepare_databases_nopyopenms.py work/broad_benchmark/data      # needs data/hye_targets.fasta
python3 -I scripts/frozen_make_entrapment.py work/broad_benchmark/data/hye_td.fasta work/broad_benchmark/entrapment/hye_entrap_td.fasta  # same for human
python3 -I scripts/make_velos_natural.py work/broad_benchmark/data/hy_td.fasta work/broad_benchmark/entrapment/hy_natural_td.fasta
DecoyDatabase -in hye_targets.fasta -out hye_ddfix_r1.fasta -entrapment_ratio 1      # ratio 2, -seed 2, -Decoy:keepPeptideNTerm false as listed below
python3 -I scripts/swap_decoys.py hye_ddfix_r1.fasta hye_ddfix_r1_rev.fasta          # control: #10364 full-reversal decoys
bash scripts/batch_frozen.sh <root> <bin>                                            # Velos natural + frozen shuffled DBs
bash scripts/batch_dd.sh <root> <bin> <arm> <hye db> <human db>                      # DecoyDatabase arms
python3 -I scripts/compare.py <bin> shuffled 1 - runs/hfx_*/frozenE                  # EntrapmentEvaluator vs frozen script
python3 -I scripts/compare.py <bin> velos 0.296196 work/broad_benchmark/data/hy_td.fasta runs/velos_*/natural
python3 -I scripts/summarize.py <bin> runs work/broad_benchmark results/entrapment
python3 -I scripts/db_stats.py <db.fasta>                                            # database properties
bash scripts/andes_count.sh runs/hfx_*/frozenE/s*/target.tsv                         # andes run.sh rule
python3 -I scripts/andes_te.py <target+entrapment fasta without decoys>              # andes build_databases.sh T/E
```

## Scripts

| script | what |
|---|---|
| `frozen_eval.py` | The frozen evaluation logic, copied verbatim: `pr-prototypes-20261001/entrapment_velos.py` (Velos natural; classification by sequence with I = L; r = 3.376) and `sage-comparison-20261002/evaluate_entrapment.py` (r = 1: entrapment if all proteins start with `ENTRAP_`, FDP = 2·N_E/N). |
| `compare.py` | Runs `frozen_eval.py` and EntrapmentEvaluator on the same 3 seeds × files and compares N, N_E and FDP row by row. |
| `diff_velos.py` | Lists the accepted Velos PSMs on which the two classifications disagree, and why. |
| `check_inproc.py` | Independent count (plain XML) on ProSE's in-process Percolator idXML, against EntrapmentEvaluator. |
| `summarize.py` | Per file and group: PSMs, entrapment PSMs, FDP with T/E from the searched database. |
| `db_stats.py` | Entrapment DB properties: same composition, moved cleavage sites, entrapment peptides without a target peptide of the same composition (mass), collisions with targets, T/E. |
| `andes_count.sh`, `andes_te.py` | The andes counting rule (`run.sh` awk) and T/E (`build_databases.sh`), verbatim, from bigbio/andes `fa0195a`. |

## Results files

| file | content |
|---|---|
| `evaluator_vs_frozen.tsv` | Question 1, row by row. |
| `velos_disagreements.txt` | The 57 Velos PSMs that EntrapmentEvaluator counted as entrapment at the PR head and the frozen script did not. 52 are I/L variants of a sample peptide; 5 are non-tryptic substrings of a yeast protein. |
| `entrapment_per_group.tsv`, `entrapment_per_file.tsv` | All arms. |
| `database_properties.tsv` | Question 2, database side. |
| `andes_rule.tsv`, `andes_te.tsv` | Question 3. |
| `default_off_sha256.txt` | DecoyDatabase without `-entrapment_ratio`, for baseline, PR head and fix (`-method reverse` and `shuffle`): byte-identical. |
| `runtime.tsv` | Runtime and peak RSS. |
