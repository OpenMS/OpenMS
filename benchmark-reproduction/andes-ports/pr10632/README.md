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

ProteoBench (min 3 runs):

| ProteomicsLFQ arm | normalization | consensus width | precursors | median \|ε\| | E. coli | yeast | human | CV |
|---|---|---|---|---|---|---|---|---|
| identified only | none | off | 42,253 | 0.2223 | 0.4048 | 0.2599 | 0.2108 | 0.2309 |
| identified only | none | **on** | 42,253 | **0.2179** | **0.2906** | 0.2515 | 0.2068 | 0.2320 |
| identified only | rt_local | off | 42,253 | 0.1494 | 0.2978 | 0.1823 | 0.1431 | 0.0879 |
| identified only | rt_local | **on** | 42,253 | **0.1390** | **0.1818** | 0.1341 | 0.1390 | 0.0926 |
| MBR (PIP-ECHO) | none | off | 46,419 | 0.2214 | 0.4090 | 0.2595 | 0.2087 | 0.2358 |
| MBR (PIP-ECHO) | none | **on** | 46,415 | **0.2168** | **0.3060** | 0.2511 | 0.2065 | 0.2361 |
| MBR (PIP-ECHO) | rt_local | off | 46,419 | 0.1572 | 0.3357 | 0.1993 | 0.1492 | 0.0951 |
| MBR (PIP-ECHO) | rt_local | **on** | 46,415 | **0.1466** | **0.1974** | 0.1434 | 0.1461 | 0.0985 |

Per-species bias, as median log2 A/B minus expected, PIP-ECHO, no normalization → rt_local:

| Species | Off | On |
|---|---|---|
| E. coli | +0.33 → +0.25 | +0.18 → +0.06 |
| yeast | −0.09 → −0.12 | +0.05 → −0.04 |

The identified-only arm behaves the same way (`hfx_species_bias.tsv`).

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

**Option off vs baseline build.** `compare_outputs.sh` compares the outputs of the two builds
(`results/identity_*.txt`).

- **Identified only:** the MSstats CSVs (none and rt_local) are byte-identical. consensusXML and mzTab are identical
  once lines with paths, the version, dates and random feature ids are removed.
- **PIP-ECHO:** not yet compared; the baseline run is in progress.

**PIP-ECHO, 4 fewer precursors with the option.** 87 and 73 of 265k MSstats rows differ. After linking,
`IDConflictResolverAlgorithm::resolve` gives an identification to the most intense of several consensus features.
New intensities change that choice in a few cases.

## andes reference (bigbio/andes#150, quoted in PR #10632)

| ProteoBench DDA QExactive, MBR quantities (andes' raw `quantity`) | andes before → after | OpenMS PIP-ECHO, none, off → on |
|---|---|---|
| median \|ε\| | 0.253 → 0.225 (−0.028) | 0.2214 → 0.2168 (−0.005) |
| E. coli \|ε\| | 0.382 → 0.276 (−0.106) | 0.4090 → 0.3060 (−0.103) |
| CV | 0.251 → 0.240 (−0.011) | 0.2358 → 0.2361 (+0.000) |

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
| `compare_outputs.sh` | Identity check against the baseline |
| `make_td_fasta.py` | Byte-identical #10364 database (sha256 `efbf5ac2…`); only used for the dropped ProSE arm |

Requirements:

- the OpenMS builds;
- `comet.exe` and `percolator` from https://github.com/OpenMS/THIRDPARTY (`Linux/x86_64`);
- a venv with `proteobench==0.18.6` and `pyarrow`.
