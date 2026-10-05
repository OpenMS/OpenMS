# IdentificationData native I/O benchmark

Build against a Release OpenMS build with its dependencies available:

```bash
cmake -S tools/benchmarks/identification_data -B ../id-benchmark \
  -DOpenMS_DIR=/absolute/path/to/OpenMS-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../id-benchmark
```

Each operation runs in its own process so peak RSS from owning construction does not
contaminate the standalone streaming measurements. Outputs are operation, match rows,
seconds and a score checksum, separated by tabs. Additional dimension and size lines
start with `#`. Every operation checks record counts. `scan` projects numeric score
columns without molecule, evidence, annotation or metadata payloads; `full-scan`
decodes every match field. The checksum sums the first score column where present;
it is a consistency check, not a comparable scientific score across different runs.

For an existing search result, one invocation separates XML parsing, owning conversion,
native writing, projected scanning and full scanning:

```bash
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark \
  idxml /data/search.idXML /data/search-native.idparquet
```

The repository contains an MS-GF+ result with modifications, target/decoy annotations,
protein sequences and repeated search metadata. From the repository root:

```bash
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark \
  idxml src/tests/topp/IDFilter_missed_cleavages_input.idXML \
  ../msgf-native.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark scan ../msgf-native.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark full-scan ../msgf-native.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark load ../msgf-native.idparquet
```

`idxml` checks that importing preserves query and match counts and that both native
scans reproduce those counts and the owning score checksum. It reports input/native
bytes, Parquet file count, run/source dimensions and resident manifest text bytes.
Its scan timings are warm immediately after writing. Its process-wide peak RSS includes
XML parsing and conversion with both representations resident, even though those
objects are released before scanning. Use the separate `scan` processes above for
streaming peak RSS. The small bundled result exercises realistic record structure;
use a complete local search result for representative throughput.

Synthetic scaling examples:

```bash
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark generate sample.idparquet 1000000 1000
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark scan sample.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark full-scan sample.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark filter sample.idparquet reduced.idparquet
/usr/bin/time -v ../id-benchmark/IdentificationDataBenchmark load sample.idparquet
```

This deliberately simple synthetic workload has two score columns, repeated molecular
strings and parent evidence, and one candidate per query. It isolates storage/scan
scaling; it is not a representative estimate of biological sequence diversity, inference
cost or compression ratios. Output directories must not already exist. Report compiler,
build flags, dataset dimensions, filesystem cache state and peak RSS with timings.
The benchmark never overwrites or removes an existing output directory.

Pooled inference can be measured separately:

```bash
OMP_NUM_THREADS=1 /usr/bin/time -v ./IdentificationDataBenchmark inference 100000 10
```

This mode uses 10,000 distinct peptidoforms and 2,000 proteins with consistent
mappings across runs, and validates the contributing-run count. Its graph
is synthetic; inference cost for other ambiguity structures can differ substantially.
See `doc/openms/identification_data_validation.md` for measured Release results.

## Comparison against legacy

`IdentificationDataLegacyBenchmark` compares idXML, existing PSM Parquet, OMS and
the owning native format with equivalent modified-peptide input. Run on Linux:

```bash
python3 tools/benchmarks/identification_data/run_legacy_comparison.py \
  ../id-benchmark/IdentificationDataLegacyBenchmark ../legacy-comparison-results
```

The output directory must not exist. Configure runtime library/data paths for your
OpenMS build first. The driver measures fresh-process peak RSS, rotates read order,
and includes one warmup plus three measured reads per format. See
`doc/openms/identification_data_legacy_comparison.md` for results and caveats.

Compare serial and threaded native I/O with the existing PSM Parquet API:

```bash
python3 tools/benchmarks/identification_data/run_legacy_comparison.py \
  ../id-benchmark/IdentificationDataLegacyBenchmark ../thread-comparison \
  --formats parquet native --cases 1000000:1 1000000:1000 \
  --native-threads 1 2 4 8
```

`--cases` accepts `PSMS:RUNS` pairs. `--read-repetitions` defaults to three measured
reads after one warmup. Each native worker count gets its own regenerated dataset;
all variants must retain the same content digest. Native I/O uses a private CPU
pool per operation and does not change Arrow's global pool. The existing Parquet
API retains its own defaults. Native single-threaded operation remains the default.
The output includes the method, raw phase timings, process CPU time and peak RSS.
Process CPU/RSS includes verification and, for writes, generation and conversion;
the separate I/O wall-clock phases exclude that work.

The additional `results/memory-footprint.json` archive contains Linux/glibc memory
diagnostics, exact object sizes for this build, public-accessor storage accounting,
and the complete standalone diagnostic source and reproduction instructions. These
measurements explain owning-load memory and compare it with a non-retaining scan;
they are separate from the repeated I/O timing matrix.
