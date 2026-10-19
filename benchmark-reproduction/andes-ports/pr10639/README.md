# PR #10639 — ProSE `database:mass_slices` (port of bigbio/andes#77, #120–#123)

Question: does the memory/speed gain andes reported for searching in precursor-mass slices also show up in OpenMS, and
do the results stay identical?

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
| andes, `phos2_E`, slice 10 / 40 / 150 Da | | | | | | 1.43 / 134.0 · 2.20 / 131.7 · 5.07 / 137.1 |

¹ From the `order` control (s1, base, s1, base per rep). In the main runs `base` came first in each rep, right after
the `s16` run, and was 1–2 s slower with 2.5× the system time (10.3 vs 8.6 s); with the order swapped the two builds
are equal. ² `chunk_size` changes the results (see below). ³ `mass_slices` is ignored with `chunk_size > 0` (as
documented): the three arms are identical to each other (and differ from the unchunked search, a pre-existing
`chunk_size` effect). ⁴ base = s1; s4 and s16 differ (documented): calibration runs on a 1/N database sample.

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
| vs. ProSE's existing options | – | the single-use index restriction already gives single-file searches andes's "only reachable peptides" index; `chunk_size` (16 chunks) beats 16 slices on `phos2_E` (3.0 vs 7.1 GiB, 70 vs 194 s) with equal PSMs/FDP but changed results | the gain over develop is only identity of results |

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
