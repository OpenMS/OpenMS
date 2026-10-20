# OpenDIA WP-00 benchmarks

Reusable local scripts for the WP-00 baseline. They use only Bash, Python's standard library, CMake, `/usr/bin/time`, and the tools required by OpenDIA.

## Build

```bash
tools/scripts/opendia_benchmark/build.sh ../openms-build
```

An existing configured build is reconfigured in place. A missing build uses `PRESET=linux-x64-release`. The script enables OpenTIMS and Thermo RAW support, builds the `OpenDIA` target, smoke-tests the executable, and records `OpenDIA.build-manifest.txt`.

All runners accept `BUILD_DIR`, `OPEN_DIA`, and `THREADS` environment overrides. Output directories must be empty, so a previous baseline cannot be overwritten accidentally.

## Synthetic TimSim 3x3

```bash
BUILD_DIR=../openms-build THREADS=4 \
  tools/scripts/opendia_benchmark/run_timsim.sh \
  /home/sing/Documents/datasets/OpenMS-TimSim-Fixture/packages/opendia_timsim_fixture_validation_3x3 \
  "$PWD/opendia-benchmark/timsim"
```

The runner verifies the packaged raw, library, truth, and metadata checksums, extracts the six `.d` archives into the output directory, and runs separate target-only and target-plus-entrapment searches. This separation is required because changing the scoring population can change quantification. Use the optional third argument `target` or `entrapment` to run one search.

`validate_timsim.py` reports target recovery, RT/IM error, realized-versus-observed precursor/protein fold-change correlation, and final-export entrapment FDP:

```bash
python3 tools/scripts/opendia_benchmark/validate_timsim.py \
  DATASET_DIR target_only/results/OpenDIA.results.tsv \
  entrapment/results/OpenDIA.results.tsv --output validation.json
```

## ProteoBench Astral entrapment

```bash
BUILD_DIR=../openms-build THREADS=4 \
  tools/scripts/opendia_benchmark/run_astral.sh \
  /home/sing/Documents/datasets/proteobench-astral-entrapment \
  "$PWD/opendia-benchmark/astral"
```

The checked dataset currently contains the 100K sampled FASTA. Set `FASTA=/path/to/full.fasta` to run another search space with the same directory layout. Extra OpenDIA arguments can follow the output directory.

`metrics_astral.py` reports per-run and pooled IDs, diagnostic `_p_target` entrapment FDP, precursor/protein completeness, median CV for complete entities, and pairwise log2-intensity correlations. These quantification metrics measure replicate repeatability, not accuracy, because this Astral directory has no quantitative ground truth.

## Captured evidence

Each run stores:

- the shell-escaped command, source revision/status, and input/binary checksums under `meta/`;
- stdout, stderr, GNU time resource data and one-second disk samples under `logs/`;
- peak output-directory bytes and a final file-size manifest under `meta/`;
- final metrics as JSON.

Historical values are comparison observations, not hard-coded pass thresholds. Use a fresh output directory for each source revision and compare metrics explicitly.

## First recorded baseline

On October 20, 2026, Release binary `037f6c02b761eb9d6ed82a2a97b3fbd5c61bd70a1993190f6f09b1bcedd8dcdc` from the dirty sandbox tree at revision `8fc1874f9bc0668d48b381554f026ef2b920dce3` completed the synthetic 3x3 baseline with four threads:

| Search | Wall time | Peak RSS | Peak output | Target recovery | Entrapment FDP |
|---|---:|---:|---:|---:|---:|
| target-only | 65.87 s | 840 MiB | 1.48 GB | 100% | n/a |
| target+entrapment | 80.28 s | 909 MiB | 1.49 GB | n/a | 0/5,938 rows |

The target-only run produced precursor and protein realized-versus-observed fold-change correlations of 0.7922 and 0.8934. The Astral metrics script was checked against the existing 100K result: 4,304 target rows, 70 diagnostic entrapment rows, 1,815 unique target precursors, and precursor pairwise log2 correlations of 0.9874 to 0.9898. A fresh Astral search was not run during this script change because it is the long-running benchmark.
