# PR #10632: ProteomicsLFQ `-consensus_elution_width` on ProteoBench DDA QExactive

This folder benchmarks OpenMS PR #10632, "ProteomicsLFQ -consensus_elution_width: one elution width per linked
precursor", a port of bigbio/andes#150. It asks whether the improvement andes reported also appears in OpenMS.
The pipeline is the one #10585 used, so the option-off rows can be checked against #10585's numbers.

## Builds

| Build | Commit | Role |
|---|---|---|
| Baseline | `77cccda7` | `git merge-base origin/develop` of the PR |
| PR | `6fa67d5f` | PR head `529b2027` plus the checkpoint fix described below |

Commit `e66bfaa1` on top of `6fa67d5f` adds only a TOPP test, so its binaries are the same.

Both builds use Ubuntu 24.04 system packages, Release, `-DWITH_GUI=OFF -DWITH_THERMO_RAW=OFF`. Each `bin/` and `lib/`
was copied to `/opt/openms-<arm>`, and `env.sh` selects the build through `BIN`.

## Data

ProteoBench module "DDA quantification – precursor ions", Q Exactive HF-X (PXD028735). It has 3 + 3 runs with
2 h gradients.

- **mzML:** ProteoBench's own files, `https://proteobench.cubimed.rub.de/raws/DDA/mzml/`, the same ones #10585 used.
  Their sha256 values are in `inputs.sha256`.
  - The mzML in `archive.openms.de/.../ProteoBench_Module_2/` is a different conversion: the spectrum count is the
    same, but the bytes and size differ. It was not used.
- **FASTA:** `ProteoBenchFASTA_MixedSpecies_HYE.fasta`, sha256 `d9ac434d…`, 31,889 entries.
  - `make_db.sh` builds the target-decoy database with `DecoyDatabase -method reverse -enzyme Trypsin` and the
    `DECOY_` prefix. Result: `HYE_td.fasta`, sha256 `63410dbe…`, as in #10585.

## Pipeline (#10585's recipe)

1. **Identification** (`search_one.sh`, per run):
   - CometAdapter 2025.01: 10 ppm, isotope error 0/1, 0.02 Da fragments, Trypsin, 2 missed cleavages, 7–40 aa,
     charge 2–5. Fixed CAM; variable Ox(M) and Acetyl(Protein N-term), at most 3.
   - PSMFeatureExtractor, then PercolatorAdapter: Percolator 3.06 from OpenMS/THIRDPARTY (#10585 used 3.08), PEP,
     TDC, seed 42.
   - Filter at q ≤ 0.10, with PEP as the main score.
   - PSMs at q < 0.01: 50,906 / 52,121 / 55,935 / 51,758 / 51,605 / 56,258. These are identical to #10585.
2. **ProteomicsLFQ** (`plfq.sh`):

   ```
   -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation -feat_dir <dir>
   ```

   All other options are defaults: star alignment, B-spline, `unique_peptides`, symmetric (Gaussian) elution model.
   - **Identified only:** `-targeted_only true`
   - **MBR:** `-targeted_only false -pip_echo true`
   - **Option:** each arm runs with `-consensus_elution_width false` and with `true`. Both runs use the same IDs and
     the same `-feat_dir` checkpoints. The off run detects; the on run reuses all 6 checkpoints.
3. **Normalization:**
   - **"none":** ProteomicsLFQ's own `-out_msstats`.
   - **"rt_local":** `-out_cxml` → `ConsensusMapNormalizer -algorithm_type rt_local` → `MSstatsConverter`.
4. **Scoring** (`pb_score.py`):
   - Official `proteobench` 0.18.6, format "quantms", min 3 runs.
   - The script reproduces the public quantms submission `d01e87b9…` exactly: 33,040 / 0.2088 / 0.3866 / CV 0.2236.
   - `pb_diag.py` adds the per-species median signed ε (bias) and median |ε − species median| (spread).

`run_arms.sh` runs every arm one after the other.

## Results (`results/`)

ProteoBench (min 3 runs). The off and on rows of each arm come from the same IDs and the same `-feat_dir`
checkpoints.

- **Identified only:** the off run detected and wrote the checkpoints. Combining them again (`T_off_ckpt`) gives
  byte-identical MSstats CSVs.
- **PIP-ECHO:** the off row is `P_off_ckpt`, combined from the checkpoints like the on row. Combining from checkpoints
  is not lossless on the PIP-ECHO path, as described below.

| ProteomicsLFQ arm | normalization | consensus width | precursors | median \|ε\| | E. coli | yeast | human | CV |
|---|---|---|---|---|---|---|---|---|
| identified only | none | off | 42,253 | 0.2223 | 0.4048 | 0.2599 | 0.2108 | 0.2309 |
| identified only | none | **on** | 42,253 | **0.2179** | **0.2906** | 0.2515 | 0.2068 | 0.2320 |
| identified only | rt_local | off | 42,253 | 0.1494 | 0.2978 | 0.1823 | 0.1431 | 0.0879 |
| identified only | rt_local | **on** | 42,253 | **0.1390** | **0.1818** | 0.1341 | 0.1390 | 0.0926 |
| MBR (PIP-ECHO) | none | off | 46,415 | 0.2214 | 0.4054 | 0.2598 | 0.2087 | 0.2357 |
| MBR (PIP-ECHO) | none | **on** | 46,415 | **0.2168** | **0.3060** | 0.2511 | 0.2065 | 0.2361 |
| MBR (PIP-ECHO) | rt_local | off | 46,415 | 0.1572 | 0.3363 | 0.1994 | 0.1491 | 0.0950 |
| MBR (PIP-ECHO) | rt_local | **on** | 46,415 | **0.1466** | **0.1974** | 0.1434 | 0.1461 | 0.0985 |

Other PIP-ECHO option-off runs, for reference:

| Run | Build | Detection | precursors | median \|ε\| | E. coli | CV |
|---|---|---|---|---|---|---|
| `P_off` | PR | in memory | 46,419 | 0.2214 | 0.4090 | 0.2358 |
| `P_base` | baseline | split over 3 processes by container restarts | 46,407 | 0.2214 | 0.4061 | 0.2357 |

With rt_local normalization these two give 0.1572 / 0.3357 / 0.0951 and 0.1572 / 0.3352 / 0.0951.

Per-species bias, as median log2 A/B minus expected, PIP-ECHO, no normalization → rt_local:

| Species | Off | On |
|---|---|---|
| E. coli | +0.33 → +0.25 | +0.18 → +0.06 |
| yeast | −0.09 → −0.12 | +0.05 → −0.04 |

The spread around each species' median also shrinks, e.g. E. coli 0.39 → 0.32 without normalization. The
identified-only arm behaves the same way (`hfx_species_bias.tsv`).

**Width fits.** "Updated" means re-quantified at the consensus width; "kept" means the feature kept its own
intensity.

| Arm | Features | Valid Gaussian fit | Updated | Kept |
|---|---|---|---|---|
| Identified only | 301,098 | 267,249 | 267,249 | 33,849 |
| PIP-ECHO | 471,838 | 381,546 | 295,961 | 46,355 |

- **Identified only:** every feature with a valid fit was updated, and the kept features are exactly those
  without one.
- **PIP-ECHO:** 83,033 consensus features have a consensus width.
- `fit_stats.py` produces the per-run counts.

**Option off vs baseline build** (`results/identity_*.txt`).

- **Identified only:**
  - The MSstats CSVs (none and rt_local) are byte-identical.
  - consensusXML and mzTab are identical once lines with paths, the version, dates and random feature ids are
    removed.
- **PIP-ECHO:** the outputs differ, but the baseline itself causes the difference, not the PR:
  - **Checkpoints:** 5 of 6 runs have identical checkpoints for the two builds (`compare_ckpt.py`: rt, m/z, charge,
    intensity, quality).
  - **Run A_01:** detection in PIP-ECHO mode (`-targeted_only false`) does not repeat. Detecting A_01 again with the
    same baseline build gives 77,858 vs 77,859 features, and about 15,500 features differ in intensity (relative
    difference around 1e-5). The PR build varies in the same way.
  - **Effect on the metrics:** at most 0.0001 in median |ε| and CV, and 0.004 in E. coli |ε|.
- **Combining from checkpoints (PIP-ECHO path).** Combining the same PIP-ECHO detection from the checkpoints
  (`P_off_ckpt`) does not reproduce the in-memory result (`P_off`) exactly either.
  - **Identified only:** only the formatting of `percolator_score` changes, and the quantities are identical.
  - **Both effects** (A_01 detection and combining from checkpoints) are already on develop and are outside this PR.

With the option on, the PIP-ECHO MSstats output has the same 265,768 rows as `P_off_ckpt`; only intensities change.

## andes reference (bigbio/andes#150, quoted in PR #10632)

| ProteoBench DDA QExactive, MBR quantities (andes' raw `quantity`) | andes before → after | OpenMS PIP-ECHO, none, off → on | verdict |
|---|---|---|---|
| median \|ε\| | 0.253 → 0.225 (−0.028) | 0.2214 → 0.2168 (−0.005) | partly: same direction, about 1/6 of the size |
| E. coli \|ε\| | 0.382 → 0.276 (−0.106) | 0.4054 → 0.3060 (−0.099) | reproduced |
| CV | 0.251 → 0.240 (−0.011) | 0.2357 → 0.2361 (±0) | not reproduced |

## Second dataset: Agilent 6550 Q-TOF (PXD069532)

Source: `archive.openms.de/openms/testfiles/20260717-ProteomicsLFQ-Multi-Instrument/`. Broccoli, conditions NPS vs
OCBS, 3 + 3 runs. The FASTA has decoys and contaminants.

- **Search:** same recipe, with FragPipe's search settings (`run_agilent.sh`).
- **Evaluation:** `agilent_eval.py`.
- **Size:** only 520–760 PSMs per run at 1 % FDR, and the FragPipe reference has 670–920. 73–76 precursors are
  quantified in all 6 runs, so this set has little power.

Results:

| arm | norm | width | precursors 6/6 | CV | Pearson log2FC vs FragPipe (n) |
|---|---|---|---|---|---|
| identified only | none | off → on | 73 | 0.330 → 0.318 | 0.13 → 0.02 (66) |
| identified only | rt_local | off → on | 73 | 0.310 → 0.285 | 0.16 → 0.08 |
| PIP-ECHO | none | off → on | 76 | 0.330 → 0.321 | 0.11 → −0.01 (67) |
| PIP-ECHO | rt_local | off → on | 76 | 0.317 → 0.278 | 0.14 → 0.04 |

- The CV improves slightly with the option.
- Fold changes agree with FragPipe in neither arm: the correlation is near 0, so the conditions barely differ for
  these precursors.
- No conclusion on accuracy.

## Scripts

| File | Purpose |
|---|---|
| `env.sh` | Paths and build selection (`BIN`) |
| `dl_pb_mzml.sh` | Download ProteoBench's mzML |
| `make_db.sh` | Target-decoy database |
| `search_one.sh`, `search_all.sh` | Identification |
| `design.tsv` | Experimental design |
| `plfq.sh`, `run_arms.sh` | ProteomicsLFQ arms and normalization |
| `pb_score.py`, `pb_diag.py` | Scoring |
| `fit_stats.py` | Width-fit counts from the checkpoints |
| `compare_outputs.sh`, `compare_ckpt.py`, `determinism.sh` | Identity checks against the baseline, checkpoint comparison, repeat detection of run A_01 |
| `make_td_fasta.py` | Byte-identical #10364 database (sha256 `efbf5ac2…`); only used for the dropped ProSE arm |

Requirements:

- the OpenMS builds;
- `comet.exe` and `percolator` from https://github.com/OpenMS/THIRDPARTY (`Linux/x86_64`);
- a venv with `proteobench==0.18.6` and `pyarrow`.
