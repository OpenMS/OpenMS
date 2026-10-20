# PR #10639 — ProSE `database:mass_slices` (port of bigbio/andes#77, #120–#123)

The section "Total peptide order" describes the current head. The tie-order caveat of the "Rebuild" section is resolved there; the sections after "Rebuild" describe the first version (`cd287f54`) and are kept as measured.

Question: does the memory/speed gain andes reported for searching in precursor-mass slices also show up in OpenMS, and
do the results stay identical?

## Total peptide order (PR head `a3464870`, decided by Timo after the rebuild)

`FragmentIndex::sortPeptides_()` sorted by (precursor m/z, protein) with an unstable sort, so the order of equally heavy
peptidoforms of one protein (positional isomers, repeats) depended on the algorithm and on which other peptides were
in the vector: a slice index could arrange them differently from the full index, and the
`scoring:max_candidates_per_spectrum` cut (ties broken by peptide index) then kept another one. `a3464870` sorts by
(precursor m/z, protein, start, length, mod bitmask), the order the `low_memory` path already used for its hits. This
is a behaviour change of the default path: wherever the candidate cut runs through such a tie, another candidate can
now be kept. `448ea1e0` (other session) moves four FragmentIndex_test sections into functions so that `main()` fits
the 1 MB Windows stack (the test crashed there since `4522ade7`).

Tests: no test reference changes. FragmentIndex_test, ProSEAlgorithm_test and all 57 TOPP_ProSE tests pass (64 ctest
entries), and the outputs of all 23 ProSE TOPP test calls are byte-identical between `4522ade7` and `a3464870` (the
`FuzzyDiff` tolerance therefore hides nothing); the test data has no tie at the cut.

**Slices are now byte-identical** (`results/v3_identity.txt`, `scripts/identity3.sh`; sha256 of every PIN and of every
idXML without date and database path):

| search | slices compared | PIN and idXML |
|---|---|---|
| frozen `hfx_A2` | 1 = 4 = 16 | identical |
| 3 files in one call | 1 = 4 = 16 | identical (and `hfx_A2` in it = `hfx_A2` alone) |
| whole `hfx_A2` run (111,413 MS2) | 1 = 4 = 16 | identical |
| `phos2`, `hye_E` | 4 = 16 (1 does not fit) | identical |
| `mhc` | 4 = 16 (1 does not fit) | identical |

**Default change, frozen #10364 suite** (`results/suite/`): all 20 files, re-sampled with the bundle's `prepare.py`
(every input sha256 = `inputs.json`; the timsTOF files without `normalize_tims_metadata.py`, whose output does not have
that hash), frozen `native_dedup:<file>:default` parameters, arms old = develop merge-base `74048ef5`
(`/opt/openms-base`), new = `a3464870`; Percolator 3.06 `-Y -U`, seeds 1/42/137, q ≤ 0.01. z = Δ / sqrt(var_new/3 +
var_old/3), a file is up/down at |z| ≥ 2. Entrapment: frozen paired-mass r = 1 databases for HF-X, Astral, TMTpro and
timsTOF (sha256 `fb77e3df…` / `76d25071…`, FDP = 2·N_E/N), protein-level 1:1 shuffle for Lumos LFQ/TMT (`hy`), frozen
Velos natural entrapment (r = 3.376) on the normal run; N and N_E summed over the files of a group.

| group | files | PSMs old | PSMs new | Δ | up / down | FDP old | FDP new |
|---|---:|---:|---:|---:|---|---:|---:|
| Velos | 3 | 8,932.7 | 8,946.7 | +0.16 % | 0 / 0 | 1.18 % | 1.20 % |
| HF-X | 3 | 13,792.7 | 13,801.0 | +0.06 % | 0 / 0 | 1.58 % | 1.59 % |
| Astral | 3 | 9,149.3 | 9,128.7 | −0.23 % | 0 / 0 | 1.40 % | 1.38 % |
| Lumos LFQ | 3 | 16,660.7 | 16,663.0 | +0.01 % | 0 / 0 | 0.96 % | 0.96 % |
| Lumos TMT | 3 | 7,371.3 | 7,372.3 | +0.01 % | 0 / 0 | 0.84 % | 0.85 % |
| Eclipse TMTpro | 3 | 7,866.7 | 7,859.3 | −0.09 % | 0 / 0 | 1.36 % | 1.37 % |
| timsTOF | 2 | 2,097.0 | 2,097.7 | +0.03 % | 0 / 0 | 0.87 % | 0.87 % |
| **total** | 20 | 65,870.3 | 65,868.7 | **−0.003 %** | **0 / 0** | | |

Per file Δ −0.44 … +0.37 %, max |z| 1.82 (`hfx_A2`, +8 PSMs). Every PIN changes, in 10–190 of 29,037–79,134 rows (1,218
of 1,104,777 rows in total); 138 spectra change their hits, and all are explained by the tie order
(`scripts/suite/tie_partner.py`, `results/suite/tie_partner.txt`): 103 keep the same peptides (only
`hyperscore_zscore` or the order of equal scores changes), and in the other 35 every peptide that is only in one
version has a partner of the same protein with the same composition and number of Ox(M) (another sequence, or the
oxidation on another M), or is the lowest-scoring of the 10 reported hits (pushed over the `report:top_hits` cut).

Larger searches (`results/v3_perc.jsonl`, same Percolator protocol): whole `hfx_A2` run 63,449 → 63,391 PSMs (−0.09 %),
FDP 1.71 → 1.68 %; `phos2`/`hye_E` (4 slices) 4,421 → 4,420, FDP 1.22 → 1.22 %; `mhc` 2,011 → 2,011.

Verdict: the total order makes `mass_slices` exactly equal to one index, and its effect on the default output is
neutral (no file up or down, total −0.003 %, entrapment FDP within ±0.04 points).

Reproduce (`<root>=/home/user/bench`; suite scripts in `scripts/suite/`, which expect `<root>/suite`):

```bash
python3 -I scripts/suite/prep_inputs.py frozen/z/prose-andes-reproduction <python with pyopenms 3.6>   # 20 inputs, hash-checked
for db in hye human; do python3 -I scripts/suite/frozen_make_entrapment.py frozen/.../data/${db}_td.fasta suite/db/${db}_frozenE_td.fasta; done
python3 -I scripts/suite/make_entrapment_db.py frozen/.../data/hy_td.fasta suite/db/hy_entrap_td.fasta
python3 -I scripts/suite/valid_params.py /opt/openms-base; python3 -I scripts/suite/valid_params.py /opt/openms-v3
python3 -I scripts/suite/suite.py /opt/openms-v3          # old/new x normal/entrapment database, interleaved per file
python3 -I scripts/suite/evaluate.py suite/results         # per_file.tsv, per_group.tsv
python3 -I scripts/suite/tie_partner.py <file's _td.fasta> runs/<file>/old_normal/input.pin runs/<file>/new_normal/input.pin
scripts/identity3.sh                                        # slices off / 4 / 16 with a3464870
```

## Rebuild of the feature (PR head `4522ade7`, decided by Timo after the first benchmark)

The first version (`cd287f54`, benchmarked below) digested the whole database and held all of its peptides in every
slice, and cut slices by spectrum count. `4522ade7` ports what andes does:

- `FragmentIndex::build(fasta, searched_spectra, bounded_generation = true)` generates only the peptides between the
  lowest and highest bound of the slice's precursor windows; a peptide none of whose peptidoforms can reach that range is
  skipped before its variable modifications are enumerated (andes#77, "mass-bounded form walk").
- Slices are cut at equal numbers of the database's peptides (`FragmentIndex::countPeptidesByMass()`, 1 Da bins, one
  pass that holds 1/slices of the database at a time). andes cuts fixed widths in Da (10–150 Da, from its memory
  budget); equal peptide counts give equal index sizes for any database and modification set (phospho, 16 slices: the
  largest slice holds 13 % of the indexed peptides, `results/v2_slice_balance.txt`).
- A warning when `mass_slices` is ignored (`chunk_size > 0`, open search, `low_memory` without enzyme), and the
  scoring time of sliced runs in the search report.

Same protocol (median of 3 interleaved runs, 4 threads, `/opt/openms-v2` = `4522ade7`; `c4`/`c16` =
`database:chunk_size` giving 4 / 16 protein chunks). Peak RSS GiB / wall s (`results/summary_v2.md`):

| search | base | off | 4 slices | 16 slices | chunk_size 4 | chunk_size 16 | other |
|---|---:|---:|---:|---:|---:|---:|---|
| frozen `hfx_A2` | 2.69 / 10.2 | 2.69 / 9.9 | **1.29 / 11.9** | **0.68 / 15.6** | 1.41 / 13.6 | 0.73 / 14.2 | |
| 3 files in one call | 4.64 / 25.5 | 4.66 / 25.5 | **2.03 / 25.3** | **1.18 / 28.8** | | 1.65 / 53.5 | |
| whole `hfx_A2` run | 5.98 / 83.0 | 5.98 / 81.7 | 4.88 / 89.6 | 4.89 / 93.2 | | 4.91 / 142.8 | spectra dominate |
| `phos2`, `hye_E` | | killed > 12.5 | **5.70 / 44.9** | **2.12 / 47.4** | 9.36 / 73.2 | 3.01 / 65.4 | andes 10/40/150 Da: 1.43 / 134 · 2.20 / 132 · 5.07 / 137 |
| `mhc` (unspecific) | | killed > 12.5 | 5.38 / 79.4 | 2.28 / 87.6 | | 2.35 / **57.6** | `low_memory` **2.06** / 164.2 |

First version for comparison: frozen 1.92 / 14.2 (4) and 1.09 / 30.6 (16); 3 files 2.87 / 30.2 and 1.90 / 46.7;
`phos2`/`hye_E` 11.77 / 73.3 and 7.08 / 194.3; `mhc` 9.29 / 81.9 and 6.43 / 208.5.

Identity: off = base (byte-identical idXML and PIN). With slices the outputs are no longer byte-identical: the peptides
are sorted by (precursor m/z, protein) with an unstable sort, so equally heavy peptidoforms of one protein (positional
isomers) can be arranged differently in a slice index than in the full one, and the candidate cut
(`scoring:max_candidates_per_spectrum`, ties broken by peptide index) can keep the other isomer
(`results/v2_pin_diff.txt`): frozen `hfx_A2` 10 / 30 of 47,986 PIN rows differ (4 / 16 slices), all only in
`hyperscore_zscore`; 3 files 10–42 rows per file, 0–5 spectra with another peptide in their top 10; whole run 252 / 181
of 663,269 rows, 14 spectra. A local build with a total sort order (`precursor_mz_, protein_idx, mod_bitmask_,
sequence_`) makes 1 and 16 slices byte-identical on frozen `hfx_A2`, the 3 files and the whole run, so this is the only
cause; it would also change the default output vs. the merge-base (frozen `hfx_A2`: 1 peptide and 20 z-scores; whole run
63,391 vs. 63,449 PSMs), so it is not in the PR.

PSMs at q ≤ 0.01 (mean of 3 seeds) and entrapment FDP (`results/percolator.jsonl`):

| search | off | 4 slices | 16 slices | chunk_size 4 / 16 |
|---|---:|---:|---:|---:|
| frozen `hfx_A2` | 4,492 (1.48 %) | 4,492 (1.47 %) | 4,492 (1.47 %) | 4,510 (1.55 %) / 4,491 (1.38 %) |
| whole run | 63,449 (1.71 %) | 63,435 (1.71 %) | 63,425 (1.69 %) | – / 63,508 (1.76 %) |
| `phos2`, `hye_E` | – | 4,421 (1.22 %) | 4,424 (1.25 %) | 4,419 (1.30 %) / 4,428 (1.26 %) |
| `mhc` | – | 2,011 | 2,011 | – / 1,987; `low_memory` 2,011 |

Verdict for the rebuild: slicing now beats `chunk_size` on memory at similar or lower time on every tryptic search
(phospho 16 slices 2.12 GiB / 47 s vs. 3.01 GiB / 65 s; 3 files 1.18 / 29 vs. 1.65 / 54) and is 2.8× faster than andes
on the same phospho input at a memory between andes's 10 and 40 Da slices. On the unspecific `mhc` search it is close
to `chunk_size` in memory (2.28 vs. 2.35 GiB) but 1.5× slower, and `low_memory` stays the smallest (2.06 GiB, 1.9×
slower). Results equal those without slices up to the tie order above (PSM counts within ±0.04 %).

## Arms and builds

| arm | build | what |
|---|---|---|
| `base` | `74048ef5` (merge-base with `develop`) | baseline |
| `s1` | `cd287f54` (PR head as benchmarked) | `database:mass_slices 1` (default, off) |
| `s2`, `s4`, `s8`, `s16` | `cd287f54` | `database:mass_slices N` |
| `cN` | `cd287f54` | `database:chunk_size` giving N protein chunks (existing alternative) |
| `fix_*` | `2b9ee0be` (PR head + the fix from this benchmark) | same, slice indices get only FragmentIndex parameters |

Release builds with Ubuntu 24.04 system packages (`-DWITH_THERMO_RAW=OFF`, Arrow 26), the same container (4 cores,
15 GB RAM, no swap). `ProSE` TOPP tool, whole process under `/usr/bin/time -v`, one search at a time, 4 threads (and 1
thread for one file), timing = median of 3 interleaved runs (`bench.py`: rep 1 of every arm, then rep 2, …).
A watchdog (`watchdog.sh`) kills a search above 12.5 GiB RSS so the container never swaps; such runs are `killed`.

## Data and search spaces

- Spectra: the frozen #10364 subsets (8,000 MS2 each, `reproduce.py inputs`, hashes match `inputs.json`): `hfx_A2`,
  `hfx_B1` (HF-X, ProteoBench M2), `astral_B1` (Astral, M8), `velos_125_R1` (Velos PXD001819, 0.5 Da), and the whole
  `hfx_A2` run (archive `LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02.mzML.gz`, 111,413 MS2). `results/inputs.sha256`.
- Databases: HYE ProteoBench FASTA (`d9ac434d…`) → `hye_td`/`hy_td` (reversed `DECOY_`, `prepare_databases.py`);
  `hye_E` = hye targets + one shuffled, mass-matched `ENTRAP_` protein per target (r = 1) + reversed decoys of both
  (`make_entrapment.py` from the #10364 supplements; 127,556 entries, sha256 `fb77e3df…`). Velos uses `hy_td`.
- `frozen`: the frozen `native_dedup:<file>:default` parameters (`gen_ini.py` copies them from `jobs.json` into a
  ProSE INI): Trypsin/P, 2 MC, 7–40, z 2–5, CAM fixed, ≤1 Ox(M), ±20 ppm, isotope error −2..1, fragment 20 ppm
  (Velos 0.5 Da), 50 candidates, top 10, deduplicate, calibration off. On `hye_E` (Velos: `hy_td`).
- `phos2`: frozen + variable Phospho (S, T, Y), ≤2 variable mods per peptide; on `hye_td` (`phos2_hye`) and on
  `hye_E` (`phos2_E`, the large space: 162 M peptidoforms before the restriction).
- `multi3`: `hfx_A2`, `hfx_B1`, `astral_B1` in one `ProSE -in` call (frozen, `hye_E`).
- `chunk`: frozen `hfx_A2` with `database:chunk_size 32000` (4 chunks), alone and with `mass_slices 4`.
- `calib`: frozen `hfx_A2` with `calibration:enabled true`.
- `mhc` (added on request of the coordinating session): immunopeptidomics, MHCquant2 timsTOF Pro class I
  `DN08_Brain_classI_techRep1.mzML` (PXD058436, archive `testfiles/id_bruker/classI/middle/`, 67,944 MS2, 59 % z = 1)
  against `human_sp.fasta` (archive `testfiles/`, 20,431 proteins) + reversed `DECOY_` decoys, unspecific cleavage
  (`enzyme "unspecific cleavage"`, `enzyme_specificity none`), length 8–12, z 1–4, no fixed mods, ≤1 Ox(M), max
  mass 2000, otherwise frozen (±20 ppm, 20 ppm fragments): 134 M peptides. Arms with the fixed build `2b9ee0be`; the
  existing alternatives `low_memory true` (peptide-mass mode; `mass_slices` is not used there) and `chunk_size 2554`
  (16 chunks). No entrapment (human only): PSMs only.

Search only: `-out_idxml` and `-out_pin`, no `-rescore` (the timings exclude rescoring). Identity: sha256 of the PIN
and of the idXML without run date, engine version, file paths and the `ProSE:` parameter echo. Where results differ:
Percolator 3.06 (OpenMS THIRDPARTY, `-Y -U`, seeds 1/42/137; the frozen protocol used 3.09) on the PIN, target PSMs at
q ≤ 0.01 (mean of seeds) and entrapment FDP = 2·N_E/N (`perc_eval.py`).

andes like-for-like (`run_andes.sh`): bigbio/andes `fa0195a2` (current main; `--candidate-index mmap` no longer exists
there, the out-of-core index is the only path), the frozen #10364 andes command, `phos2` mods
(`andes_mods_phos2.txt`, `--max-mods 2`), `hye_E`, `hfx_A2`, 4 threads, `--fragment-index-slice-da` 10 / 40 / 150.

## Results (4 threads unless noted; `results/summary.md` has every arm)

Peak RSS GiB / wall s (median of 3). All arms of a row give byte-identical idXML and PIN unless marked.

| search | base | s1 (PR, off) | s2 | s4 | s8 | s16 |
|---|---:|---:|---:|---:|---:|---:|
| frozen `hfx_A2` | 2.69 / 9.3¹ | 2.68 / 9.2¹ | 2.36 / 10.9 | 1.92 / 14.2 | | 1.09 / 30.6 |
| frozen `hfx_B1` | 2.70 / 10.0 | 2.70 / 8.8 | 2.36 / 10.7 | 1.89 / 13.6 | | 1.09 / 30.7 |
| frozen `astral_B1` | 2.76 / 14.9 | 2.77 / 13.8 | 2.54 / 16.3 | 2.13 / 19.3 | | 1.47 / 35.7 |
| frozen `velos_125_R1` (`hy_td`) | 1.40 / 9.0 | 1.40 / 8.9 | 1.20 / 10.6 | 0.96 / 11.4 | | 0.77 / 19.2 |
| frozen `hfx_A2`, 1 thread | 2.64 / 33.0 | 2.64 / 33.9 | | 1.87 / 47.8 | | 1.04 / 112.6 |
| `multi3` (3 files, one call) | 4.65 / 27.6 | 4.64 / 26.9 | | 2.87 / 30.2 | | 1.90 / 46.7 |
| whole `hfx_A2` run (111k MS2) | 5.98 / 82.5 | 5.98 / 80.0 | | 4.88 / 95.0 | | 4.87 / 109.0 |
| `phos2_hye` | 7.00 / 22.2 | 7.00 / 22.5 | | 5.96 / 37.1 | | 3.62 / 94.0 |
| `phos2_E` | | > 12.5 (killed ×3) | | 11.77 / 73.3 | 9.50 / 117.6 | 7.08 / 194.3 |
| `phos2_E`, `chunk_size` (c4 / c16)² | | | | 9.35 / 78.6 | | **3.03 / 70.5** |
| frozen `hfx_A2`, `chunk_size 32000`³ | 1.40 / 13.6 | 1.41 / 13.6 | | 1.41 / 13.3 | | |
| frozen `hfx_A2`, calibration on⁴ | 4.13 / 13.0 | 4.13 / 11.6 | | 1.92 / 16.2 | | 1.11 / 32.1 |
| `mhc` (unspecific, 68k spectra)⁵ | killed (> 12.5) | killed (> 12.5) | | 9.29 / 81.9 | | 6.43 / 208.5 |
| `mhc`, `low_memory` / `chunk_size` 16 chunks⁵ | | 2.05 / 163.4 · 2.35 / 57.9 | | | | |
| andes, `phos2_E`, slice 10 / 40 / 150 Da | | | | | | 1.43 / 134.0 · 2.20 / 131.7 · 5.07 / 137.1 |

¹ From the `order` control (s1, base, s1, base per rep). In the main runs `base` came first in each rep, right after
the `s16` run, and was 1–2 s slower with 2.5× the system time (10.3 vs 8.6 s); with the order swapped the two builds
are equal. ² `chunk_size` changes the results (see below). ³ `mass_slices` is ignored with `chunk_size > 0` (as
documented): the three arms are identical to each other (and differ from the unchunked search, a pre-existing
`chunk_size` effect). ⁴ base = s1; s4 and s16 differ (documented): calibration runs on a 1/N database sample. ⁵ s4 = s16 (identical);
`low_memory` and `chunk_size` give other idXML/PIN. s16 times 174 / 303 / 209 s (the 303 s run followed a worker
restart).

Base = PR with `mass_slices 1` everywhere (byte-identical idXML/PIN): nothing changes by default.

PSMs at q ≤ 0.01 (mean of 3 seeds) and entrapment FDP, `results/percolator.jsonl`:

| search | arm(s) | PSMs | FDP |
|---|---|---:|---:|
| frozen `hfx_A2` / `hfx_B1` / `astral_B1` / `velos` | all (identical) | 4,492 / 4,362 / 2,953 / 2,977 | 1.48 / 1.57 / 2.10 / – % |
| whole `hfx_A2` | all (identical) | 63,449 | 1.71 % |
| `phos2_hye` | all (identical) | 4,621 | – |
| `phos2_E` | s4 = s8 = s16 | 4,420 | 1.21 % |
| `phos2_E` | c4 / c16 | 4,419 / 4,428 | 1.30 / 1.26 % |
| `phos2_E` | andes (10 = 40 = 150 Da, identical PIN) | 4,450 | 1.54 % |
| frozen `hfx_A2`, `chunk_size 32000` | c4 | 4,498 | 1.42 % |
| frozen `hfx_A2`, calibration | s1 / s4 / s16 | 4,341 / 4,472 / 4,499 | 1.63 / 1.49 / 1.79 % |
| `mhc` | s4 = s16 / `low_memory` / `chunk_size` 16 | 2,011 / 2,011 / 1,987 | – |

With calibration, the 1/N sample gives fewer calibration PSMs (s1 268, s4 160, s16 94) and wider windows (s1
precursor −8.3/+8.6 ppm; s4 −18.0/+20.0 ppm, fragment 14.8 ppm; s16 −19.1/+19.3, fragment 17.9 ppm). Here the wider
windows find 3–4 % more PSMs; with s16 the FDP is 0.2 points higher than with s1.

### Why the memory falls less than 1/N (`results/slice_balance.txt`)

1. Every slice digests the whole database and holds all of its peptides before `FragmentIndex::build` keeps those in
   the slice's precursor windows (19.5 M peptides per slice for `frozen` on `hye_E`, 81 M for `phos2_hye`, 162 M for
   `phos2_E`). Only the fragments shrink with the slices. The `phos2_E` search without slices is killed during this
   step, and `chunk_size` (which splits the database itself) needs 3.0 GiB with 16 chunks where 16 slices need 7.1 GiB.
   andes (#77, "mass-bounded form walk") only enumerates the forms inside a slice's mass window.
2. Slices hold equal numbers of spectra, not of peptides. The top slice spans the heavy precursors, where the
   peptidoforms are: with 4 slices it holds 43 % (`frozen`) to 68 % (`phos2_hye`) of the reachable peptides, with 16
   slices 13 % to 31 %. Peak memory follows the largest slice.
3. The time cost is one database digestion per slice: +28–65 % (s4) and +116–256 % (s16) on the `frozen` files,
   +65 / +318 % on `phos2_hye`; less where scoring dominates (`multi3` +12 / +74 %, whole run +19 / +36 %). andes's time is flat in the slice width (132–137 s for 10–150 Da).
On `mhc` the slices reach 94 % of the 134 M peptides in total, and with 4 slices the two upper slices hold 57 M and
52 M, so 4 slices still need 9.3 GiB; `low_memory` (2.05 GiB) and `chunk_size` (2.4 GiB, 3.6× faster than 16 slices)
are the better options for unspecific searches.

4. On the whole run, the spectra held in memory (all files, for the whole search) dominate beyond 4 slices
   (4.88 → 4.87 GiB from 4 to 16 slices).

### andes reference

bigbio/andes#77 (`9ed6bc0`): with the out-of-core candidate index, a file-order chunk "was the whole database (100 GB
on phospho)"; mass-ordered slices made the phospho search (PXD007653, mouse + 1:1 entrapment, 4 variable mods) run:
"a 150 Da slice peaks at ~25 GB RSS, about 0.13 GB per Da over a ~5 GB base", and 80 Da slices give 13.7 GB (andes
README, 637 s wall at 8 threads; #77 title: 169 s vs Comet 226 s, +9.6 % PSMs). On our `phos2_E` input andes shows
the same shape: 1.43 → 5.07 GiB from 10 to 150 Da at constant time and identical output.

| metric | andes | OpenMS (this PR) | verdict |
|---|---|---|---|
| results identical with slicing | yes (byte-identical PIN for 10/40/150 Da here) | byte-identical idXML and PIN in all searches without calibration (single file, 3 files, whole run, phospho), also vs. base | reproduced |
| peak memory falls with the slices | ~linear in slice width (0.13 GB/Da); search becomes feasible | −8…−14 / −23…−31 / −45…−60 % for 2/4/16 slices (`frozen`, 4 files); −38/−59 % (`multi3`, s4/s16); −15/−48 % (`phos2_hye`); −18 % (whole run); makes `phos2_E` (> 12.5 GiB) run in 11.8 GiB (s4) / 7.1 GiB (s16) | partly: direction yes, size smaller (full digestion per slice, slices balanced by spectra) |
| time | flat in slice width | +12–65 % at 4 slices, +36–320 % at 16 (1 thread: +41 / +232 %) | not reproduced (one full digestion per slice) |
| vs. ProSE's existing options | – | the single-use index restriction already gives single-file searches andes's "only reachable peptides" index; `chunk_size` (16 chunks) beats 16 slices on `phos2_E` (3.0 vs 7.1 GiB, 70 vs 194 s) with equal PSMs/FDP but changed results; unspecific `mhc`: `low_memory` 2.05 GiB / 163 s and `chunk_size` 2.35 GiB / 58 s vs. 6.43 GiB / 209 s with 16 slices | the gain over develop is only identity of results |

## Fix pushed to the PR

`2b9ee0be`: `scoreMassSlices_()` passed all ProSE parameters to every slice index, and FragmentIndex warned about each
unknown one at every build and query: 928 warning lines for one file with 16 slices (29 per index on develop). Now 0;
idXML/PIN, RSS and time unchanged (`fix_frozen`, `fix_multi3`); ProSEAlgorithm_test and the 57 TOPP_ProSE tests pass.

## Proposals (not done here)

1. Generate only the peptides of a slice's mass range (bound the mass walk in `FragmentIndex::generatePeptides`, as
   andes does), so each slice costs 1/N of the digestion in time and memory.
2. Cut slices by reachable peptides (or a fixed precursor-mass span, as andes's 10–150 Da) instead of equal spectrum
   counts.
3. Qualify the documentation: "the index takes about 1/slices of the memory" holds for the fragments only.
4. Log a line when `mass_slices` is ignored (`chunk_size > 0`, open search, `low_memory` without enzyme).
5. The search report shows no `search` time for sliced runs.

## Reproduce

```bash
# inputs (frozen bundle prose-andes-reproduction.zip, sha256 e24cbe74…; lxml + pyopenms for prepare_databases.py)
python3 reproduce.py inputs hfx_A2 --keep-raw && python3 reproduce.py inputs hfx_B1 astral_B1 velos_125_R1
python3 scripts/make_entrapment.py broad_benchmark/data/hye_td.fasta db/hye_E.fasta
gunzip -c broad_benchmark/raw/hfx_A2.mzML.gz > whole/hfx_A2_full.mzML
# builds: /opt/openms-base (74048ef5), /opt/openms-pr (cd287f54), /opt/openms-fix (2b9ee0be): bin/ lib/ share/
ProSE -write_ini prose.ini      # per build; then per dataset / space:
python3 scripts/gen_ini.py prose.ini jobs.json hfx_A2 ini/hfx_A2.ini
python3 scripts/gen_ini.py prose.ini jobs.json hfx_A2 ini/hfx_A2_phos2.ini \
  "modifications:variable=Oxidation (M),Phospho (S),Phospho (T),Phospho (Y)" modifications:variable_max_per_peptide=2
scripts/watchdog.sh &
bash scripts/queue1.sh; bash scripts/queue2.sh; bash scripts/queue3.sh     # experiments.json holds every arm
python3 scripts/bench.py scripts/experiments.json fix_frozen 3
python3 scripts/summarize.py results/*.tsv
python3 scripts/perc_eval.py runs/<exp>/<arm>_r1/<file>.pin runs/<exp>/<arm>_r1/perc
```

A single search, e.g. `phos2_E` with 16 slices:

```bash
ProSE -ini ini/hfx_A2_phos2.ini -in hfx_A2.mzML -database hye_E.fasta -out_idxml hfx_A2.idXML -out_pin hfx_A2.pin \
  -threads 4 -no_progress -Search:database:mass_slices 16
```
