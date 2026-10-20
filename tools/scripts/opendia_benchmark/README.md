# OpenDIA WP-00 benchmarks

Reusable local scripts for the WP-00 baseline. The benchmark runners and metric scripts use Bash, Python's standard library, CMake, `/usr/bin/time`, and the tools required by OpenDIA. PDF rendering also uses Python `markdown` and Chrome or Chromium.

## Build

```bash
tools/scripts/opendia_benchmark/build.sh ../openms-build
```

An existing configured build is reconfigured in place. A missing build uses `PRESET=linux-x64-release`. The script enables OpenTIMS and Thermo RAW support, builds the `OpenDIA` target, smoke-tests the executable, and records `OpenDIA.build-manifest.txt`.

All runners accept `BUILD_DIR`, `OPEN_DIA`, and `THREADS` environment overrides. Output directories must be empty, so a previous baseline cannot be overwritten accidentally.

Both runners use the same speed settings from `lib.sh`: `readOptions=cacheWorkingInMemory`, `batchSize=0`, `innerBatchSize=-1`, `maxConcurrentSwaths=-1`, and `outer_loop_threads=-1`. These enable automatic wave scheduling. The Astral-only limit of one SWATH map is removed. Instrument and library settings remain different because the inputs are different. These settings do not prove maximum speed, and the full Astral run can exceed this host's RAM.

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

## Final report

```bash
python3 tools/scripts/opendia_benchmark/make_report.py CURRENT_TIMSIM_DIR \
  --previous-timsim PREVIOUS_TIMSIM_DIR \
  --develop-revision DEVELOP_COMMIT \
  --astral-metrics ASTRAL_METRICS.json \
  --astral-time ASTRAL_TIME.txt \
  --output-dir tools/scripts/opendia_benchmark
```

The command writes `FINAL_REPORT.md` and five plots: run time, peak memory, identification, quantification, and FDR diagnostics. Use fresh validator JSON files so they contain `fdr_curve`. The report uses short sentences and simple terms.

FDR plots show raw entrapment shares at five q-value limits. They are not search-space-adjusted FDR estimates. Final exports exclude rejected results, so these plots cannot prove full FDR calibration. TimSim and Astral can use different score columns, which the report names. Astral quantification measures repeatability, not accuracy.

Render the report and plots to PDF with the installed Python `markdown` package and Chrome or Chromium:

```bash
python3 tools/scripts/opendia_benchmark/render_report.py \
  tools/scripts/opendia_benchmark/FINAL_REPORT.md
```

Run the small regression check with `python3 tools/scripts/opendia_benchmark/test_benchmark.py`. Its PDF failure check also needs Python `markdown`.

## Captured evidence

Each run stores:

- the shell-escaped command, source revision/status, and input/binary checksums under `meta/`;
- stdout, stderr, GNU time resource data and one-second disk samples under `logs/`;
- peak output-directory bytes and a final file-size manifest under `meta/`;
- final metrics as JSON.

Historical values are comparison observations, not hard-coded pass thresholds. Use a fresh output directory for each source revision and compare metrics explicitly.

## Current result

See [`FINAL_REPORT.md`](FINAL_REPORT.md) or [`FINAL_REPORT.pdf`](FINAL_REPORT.pdf) for the short result and plots. Astral values remain historical until a full run completes with the shared settings.
