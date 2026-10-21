# PR #10634: ProSE `report:collapse_positional_isomers`

Benchmark of the option that reports only the best of a spectrum's hits with the same charge, residues and
modifications (positional isomers). It answers whether what bigbio/andes#143 reported for the same idea also holds in
OpenMS, and whether the collapse changes the FDR estimate.

The results tables are in `results/`. The PR description has the summary.

## Builds

| arm | commit | notes |
|---|---|---|
| baseline | `74048ef5` (merge-base of the PR with develop) | the option does not exist |
| PR | `5639017e` (PR head) | option off (default) and on |

Ubuntu 24.04 with system packages, Release, `-DWITH_GUI=OFF -DWITH_THERMO_RAW=ON` (only to convert the phospho
`.raw`), targets `ProSE FileConverter DecoyDatabase ProSEAlgorithm_test`; `ProSEAlgorithm_test` passes.
4 cores, 15 GB RAM. Both builds were copied to `/opt/openms-{base,pr}` (bin, lib, share).

## Data and search

**Frozen #10364 suite (≤ 1 Ox(M)).** All 20 files (8,000 MS2 each), re-sampled with the bundle's `prepare.py`
(`prep_inputs.py` drives it and resumes PRIDE downloads); every sampled mzML and all three target/decoy FASTAs match
`inputs.json` / `databases.json` (HYE source sha256 `d9ac434d…`). ProSE TOPP runs with the frozen
`native_dedup:<dataset>:default` parameters from `jobs.json`, passed as `-Search:<key>` (`run_arm.py`); `snes_enabled`
and `precursor:isotopes` do not exist in this build and are skipped. ProSE writes `-out_pin` without rescoring.

**Phospho (andes's phospho benchmark).** PXD007653 `20151014_QEp6_NaKr_SA_totalliver_control2_phospho.raw` (SHA-1
matches PRIDE), converted with `FileConverter -RawToMzML:reader inprocess` (vendor centroiding; 102,819 MS2). Database
as in andes `docs/benchmarks`: UniProt mouse reviewed (UP000000589, 17,277 entries, downloaded 2026-10-19) plus a 1:1
shuffled entrapment twin per protein (andes `glyco/build_shuffled_entrap.py`, seed 1), plus full-reversal `DECOY_`
decoys of both (`mouse_entrap_td.fasta`, sha256 in `results/inputs.sha256`). Search: the frozen high-resolution
parameters (Trypsin/P, 2 MC, length 7–40, 20/20 ppm, top 10, candidate cap 50) with the andes phospho mod set:
fixed CAM(C), variable Oxidation(M), Acetyl(Protein N-term), Phospho(S/T/Y), at most 3 per peptide (andes: 4), and
`database:chunk_size 8000` (unchunked, the 211 M peptidoforms need ~45 GB). `results/phospho_params.json` has all
parameters.

**Rescoring.** External Percolator 3.09.0 (bioconda `percolator-3.9-he32743f_1`, the frozen version) with
`-Y -U --seed {1,42,137}` on each PIN; target PSMs at q ≤ 0.01; mean over seeds per file, mean over files per group;
z = Δ / sqrt(var_on/3 + var_off/3), a file is up/down at |z| ≥ 2. Peptides = distinct residue sequences of the accepted
PSMs, peptidoforms = distinct modified sequences.

**Entrapment.** (1) the frozen designs: Velos natural entrapment (frozen `entrapment_velos.py` via
`../pr10629/scripts/frozen_eval.py velos`, r = 3.376) on the normal runs, and paired-mass shuffled r = 1 databases
(`../pr10629/scripts/frozen_make_entrapment.py`, sha256 `fb77e3df…` hye, `76d25071…` human) for HF-X, Astral and
TMTpro (arms `*_frozenE`); (2) a protein-level 1:1 shuffle (andes's design, `make_entrapment_db.py`) for all 20
files (arms `*_entrap`). FDP = 2·N_E/N (r = 1); a PSM is entrapment when all its proteins are `ENTRAP_`.

## Reproduce

Layout `<root>=/home/user/bench`: the frozen bundle unpacked into `<root>/frozen/`, Percolator 3.09 wrapper at
`<root>/percolator`, builds in `/opt/openms-base` and `/opt/openms-pr`, each with `valid_params.json` (the `Search:`
parameter names of `ProSE -write_ini`).

```bash
python3 -I scripts/prep_inputs.py <root>/frozen <python with pyopenms>          # 20 inputs + 3 FASTAs, hash-checked
for db in hy hye human; do python3 -I scripts/make_entrapment_db.py frozen/broad_benchmark/data/${db}_td.fasta frozen/broad_benchmark/data/${db}_entrap_td.fasta; done
for db in hye human; do python3 -I ../pr10629/scripts/frozen_make_entrapment.py frozen/broad_benchmark/data/${db}_td.fasta frozen/broad_benchmark/data/${db}_frozenE_td.fasta; done
python3 -I scripts/suite.py                     # arms base, off, on, off_entrap, on_entrap, interleaved per file
python3 -I scripts/suite_frozenE.py             # off_frozenE, on_frozenE (HF-X, Astral, TMTpro)
python3 -I scripts/evaluate.py runs/suite results --require-complete
python3 -I ../pr10629/scripts/frozen_eval.py velos frozen/broad_benchmark/data/hy_td.fasta runs/suite/velos_*/{off,on}/s*/target.tsv
# phospho
FileConverter -in control2_phospho.raw -out phospho/control2_phospho.mzML -RawToMzML:reader inprocess -process_lowmemory
python3 -I scripts/build_shuffled_entrap.py mouse.fasta mouse_entrap.fasta 1   # then append DECOY_ reversals of all entries
scripts/phospho.sh off on base
python3 -I scripts/phospho_table.py runs/phospho
python3 -I scripts/rank1_check.py runs/phospho/off/input.pin runs/phospho/on/input.pin
python3 -I scripts/rank_switch.py runs/phospho/off/input.pin runs/phospho/off/s42/target.tsv
scripts/timing.sh > results/timing.tsv          # 3 interleaved repetitions per arm, 4 threads
```

## Results

| file | content |
|---|---|
| `per_file.tsv` | suite per file: base/off/on PSMs (seed mean), Δ, z, call, seeds, peptides, PIN rows, isomer rows, base == off hashes |
| `per_group.tsv` | suite per group and total |
| `entrapment.tsv` | combined FDP per group, both entrapment designs |
| `velos_natural_entrapment.txt` | Velos natural entrapment, off and on |
| `phospho.tsv` | phospho per arm and seed: PSMs, peptides, peptidoforms, phospho PSMs, entrapment, PIN rows, isomer sets |
| `phospho_rank1_check.txt` | best PIN row per spectrum is identical in off and on (peptide, label, all features) |
| `phospho_site_switch.txt` | accepted PSMs that are not ProSE's rank 1 of their spectrum, and how many are a site isomer of it |
| `timing.tsv` | search wall time, peak RSS, Percolator wall time |
| `inputs.sha256` | inputs |

| `mhc_site.tsv` | spike-in site accuracy per arm, file and seed (see below) |
| `mhc_wrong_sites.txt` | wrong-site PSMs per arm/file/seed: true → reported label residue, and the peptides |
| `mhc_params.json`, `mhc_inputs.sha256` | spike-in search parameters and inputs |

## Known-position test: MHCquant spike-in (PXD011628)

The phospho benchmark has no ground truth for the site. The MHCquant spike-in has: 66 synthetic HLA class I peptides,
each with exactly one heavy residue (V or P: `Label:13C(5)15N(1)`, +6.0138; L: `Label:13C(6)15N(1)`, +7.0172; table
`spiked_peptides_labels.tsv` from https://archive.openms.de/openms/testfiles/lfq/, from Bichmann et al. 2019 Table S3),
spiked into a JY immunopeptidome. Files: the three 100 fmol replicates `…TUMAPs_100fm_…msms52/53/54.raw` (Lumos,
FTMS CID MS2), converted with `FileConverter -RawToMzML:reader inprocess`.

Search (`results/mhc_params.json`): human Swiss-Prot `human_sp.fasta` from the same archive (20,431 entries) +
reversed `DECOY_` decoys; unspecific cleavage, length 8–12, charge 1–4, `low_memory true`, 10 ppm precursor, 20 ppm
fragment, no fixed mods, variable Oxidation(M), Label:13C(5)15N(1) (V), Label:13C(5)15N(1) (P), Label:13C(6)15N(1) (L),
≤ 2 per peptide, `report:top_hits 10`; other parameters as the frozen high-res set. Same Percolator protocol.
V and P labels have the same modification name, so the PR's key treats them as isomers of each other.

Evaluation (`mhc_site.py`): accepted PSMs (q ≤ 0.01) whose residue sequence is a spiked peptide; correct = one heavy
label at the true position; wrong site = one heavy label elsewhere. Also ProSE's own rank-1 row of the same spectrum
(best `score` in the PIN) and, where Percolator accepted another row than ProSE's rank 1, which of the two was right.

```bash
scripts/mhc.sh                                         # off/on for msms52/53/54
python3 -I scripts/mhc_site.py spiked_peptides_labels.tsv runs/mhc/{off,on}_5{2,3,4}
python3 -I scripts/mhc_wrong.py spiked_peptides_labels.tsv runs/mhc/{off,on}_5{2,3,4}/s*/target.tsv
```

Result, sum over the three files of the seed means: off 4,308.7 correct / 133.0 wrong site (97.0%), on 4,248.7 /
180.7 (95.9%); wrong-site PSMs +36%, up in all 9 file × seed pairs. All accepted PSMs 23,696.7 → 23,761.3 (+0.3%).
With the option off, Percolator accepted 111.7 spike-in PSMs that are not ProSE's rank-1 row; 89.0 of them had the
correct site (ProSE's rank 1: 27.7). The collapse removes those alternatives, so ProSE's rank-1 site errors stay.
Most of the change is one peptide, VLLAGFKPPL (label on P8 next to P9): 17–22 → 34–47 wrong-site PSMs per file (seed 42).

Isolation: with the option off, PIN and idXML (date/version lines removed) are byte-identical to the baseline on all 20
suite files and on the phospho file.
