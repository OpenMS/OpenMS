# PercolatorAdapter on internal PSM parquet (idparquet)

**Issue:** [OpenMS/OpenMS#9225](https://github.com/OpenMS/OpenMS/issues/9225)
**Date:** 2026-04-28
**Status:** Design approved, awaiting implementation plan.

## Background

Issue #9225 asks PercolatorAdapter to read and write PSM data in parquet
form, with an FDR filter and a "best PSM per spectrum" option. The issue
text uses "qpx" loosely; this design uses the **internal** PSM parquet
schema (`PSMSchema`, lossless round-trip), not the cross-tool quantms
exchange schema (`QPXPSMSchema`, lossy). Both are declared in
`src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`.

The change extends beyond PercolatorAdapter: the broader pipeline —
search-engine adapters and IDMerger — also gain native parquet IO, so
the entire identification stage can run on the internal parquet format
without converting through idXML.

## Decisions

| Question | Decision |
|---|---|
| Scope | Full bundle: PercolatorAdapter + IDMerger + search-engine adapters |
| FDR filter behaviour | Output post-filter on Percolator's q-value |
| Best vs all PSMs | Default best per spectrum, flag for all passing |
| Picked protein FDR | Out of scope — applied via `ProteinFDR` chained downstream |
| Container | Directory containing 4 parquet files (PSMSchema + 3 protein tables) |
| File extension | `.idparquet`, register `FileTypes::IDPARQUET` |
| Search-engine integration | Extend each adapter's `out` valid-formats |
| Schema | `PSMSchema` (lossless), not `QPXPSMSchema` |

`QPXConverter` is not modified; it remains the producer of the
quantms-exchange parquet (`QPXPSMSchema`). The new path is parallel.

## Architecture

Three layers, narrowest at the top.

### 1. Format layer

New class `PSMArrowIO` in
`src/openms/{include,source}/OpenMS/FORMAT/PSMArrowIO.{h,cpp}`:

- `static bool exportToParquet(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&, const String& dir, bool export_all_psms = true, const ParquetWriteConfig& = {})`
- `static bool importFromParquet(const String& dir, std::vector<ProteinIdentification>&, PeptideIdentificationList&)`

Internally:

- Export delegates to `QPXFile::exportToArrow` (PSMSchema) for
  `psms.parquet`, and to `ProteinIdentificationArrowIO::export*` for
  `proteins.parquet`, `protein_groups.parquet`, `search_params.parquet`.
- Import delegates to `ProteinIdentificationArrowIO::importFromParquet`
  for the protein side, and to a new free-standing PSMSchema reader
  extracted from `ConsensusMapArrowIO::importPSMsFromArrow` (the half
  that does not depend on `ConsensusFeature` linkage). The extracted
  helper is exposed as `static bool QPXFile::importFromArrow(table, protein_ids, peptide_ids)`.
  `ConsensusMapArrowIO` is refactored to call this helper for the
  PSM-construction half so we maintain a single PSM reader.

### 2. Dispatch layer

`FileHandler::loadIdentifications` and `storeIdentifications` already
dispatch on extension to idXML / mzIdentML / oms. We add an `IDPARQUET`
branch calling `PSMArrowIO::{import,export}FromParquet`. New
`FileTypes::IDPARQUET` enum value with `.idparquet` extension and the
corresponding `nameToType` / `typeToName` mappings in
`src/openms/source/FORMAT/FileTypes.cpp`.

### 3. Tool layer

Every TOPP tool that already takes idXML/mzid as input or writes them
as output adds `"idparquet"` to its `setValidFormats_` call. The
`FileHandler` dispatch makes that a one-line change per tool.

Tools touched:

- `PercolatorAdapter` — input + output
- `IDMerger` — input + output (uses a `formats` list)
- `CometAdapter`, `SageAdapter`, `MSGFPlusAdapter`, `MSFraggerAdapter`,
  `XTandemAdapter` — output

Most of these adapters already call `FileHandler().storeIdentifications`,
but with an explicit allowed-types whitelist (typically
`{FileTypes::IDXML}`). Adding parquet support requires expanding that
whitelist to include `FileTypes::IDPARQUET` at each call site.
**Exception:** `SageAdapter` currently calls `IdXMLFile().store(...)`
directly — it needs to be migrated to
`FileHandler().storeIdentifications(...)` first (a minor refactor).
`PercolatorAdapter` already passes
`{FileTypes::IDXML, FileTypes::MZIDENTML}` and just needs `IDPARQUET`
added.

PercolatorAdapter additionally gets two new options:

- `score:fdr` (double, default `1.0`, range `[0.0, 1.0]`) — post-filter
  cutoff on Percolator's q-value before writing output.
- `keep_all_passing` (flag, default off) — keep all PSMs that pass the
  cutoff, not just the best per spectrum.

## File inventory

### New files

- `src/openms/include/OpenMS/FORMAT/PSMArrowIO.h`
- `src/openms/source/FORMAT/PSMArrowIO.cpp` (~150–250 lines)
- `src/tests/class_tests/openms/source/PSMArrowIO_test.cpp`

### Modified files (format / dispatch)

- `src/openms/include/OpenMS/FORMAT/FileTypes.h` — add `IDPARQUET`
  enum value (alongside `PARQUET`, `CHROMPARQUET`, `MOBILPARQUET`).
- `src/openms/source/FORMAT/FileTypes.cpp` — register `.idparquet`
  in `nameToType` / `typeToName` / properties.
- `src/openms/include/OpenMS/FORMAT/QPXFile.h` and `.cpp` — add
  `static bool importFromArrow(...)` extracted from
  `ConsensusMapArrowIO::importPSMsFromArrow`.
- `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp` — call the
  extracted helper to avoid duplication. Light refactor; behaviour
  preserved by the existing `ConsensusMapArrowIO` tests.
- `src/openms/include/OpenMS/FORMAT/FileHandler.h` and
  `src/openms/source/FORMAT/FileHandler.cpp` — extend
  `loadIdentifications` / `storeIdentifications` with `IDPARQUET`
  branch.
- `src/openms/includes.cmake` and `src/openms/sources.cmake` —
  register the new header/source.
- `src/tests/class_tests/openms/executables.cmake` — register the new
  test.

### Modified files (tools)

- `src/topp/PercolatorAdapter.cpp` — add `idparquet` to `in`/`out`
  valid-formats, expand the allowed-types whitelist on the
  `FileHandler().storeIdentifications(...)` call (currently
  `{FileTypes::IDXML, FileTypes::MZIDENTML}`) to include
  `FileTypes::IDPARQUET`, add `score:fdr` and `keep_all_passing`
  options, apply the post-filter step before writing the output.
- `src/topp/IDMerger.cpp` — add `idparquet` to its `formats` list and
  to the allowed-types whitelist on the `storeIdentifications` call.
- `src/topp/CometAdapter.cpp`, `MSGFPlusAdapter.cpp`,
  `MSFraggerAdapter.cpp`, `XTandemAdapter.cpp` — add `idparquet`
  to `out` valid-formats and expand the allowed-types whitelist on
  the `FileHandler().storeIdentifications(...)` call to include
  `FileTypes::IDPARQUET`.
- `src/topp/SageAdapter.cpp` — same as above, plus migrate the direct
  `IdXMLFile().store(...)` call (currently at line 889) to
  `FileHandler().storeIdentifications(...)` so it picks up the dispatch.

## Data flow

### Flow A — End-to-end pipeline

```
CometAdapter --in run1.mzML --out run1.idparquet
SageAdapter  --in run2.mzML --out run2.idparquet
IDMerger --in run1.idparquet run2.idparquet --out merged.idparquet
PercolatorAdapter --in merged.idparquet --out scored.idparquet --score:fdr 0.01
```

Inside each search-engine adapter:

1. Producer assembles `protein_ids` + `peptide_ids` in memory (unchanged).
2. `FileHandler::storeIdentifications("run1.idparquet", protein_ids, peptide_ids)`
   detects `IDPARQUET`, calls `PSMArrowIO::exportToParquet`.
3. `PSMArrowIO` creates the directory; calls `QPXFile::exportToArrow`
   (writes `psms.parquet`); calls `ProteinIdentificationArrowIO::export*`
   (writes the three protein tables).

In IDMerger:

1. `FileHandler::loadIdentifications` on each input → in-memory
   `protein_ids` + `peptide_ids` per input.
2. Existing merge logic runs (unchanged) — produces a single combined
   `protein_ids` (one entry per source run, distinguished by
   `run_identifier`) and combined `peptide_ids`.
3. Store combined result via `FileHandler::storeIdentifications`. The
   resulting `merged.idparquet/` has multiple rows in `proteins.parquet`
   and `search_params.parquet` (one per run); PSMs reference them by
   `run_identifier`.

In PercolatorAdapter:

1. `FileHandler::loadIdentifications("merged.idparquet")` → fully
   reconstructed `protein_ids` (with `extra_features`, `enzyme`,
   `charges`, fixed/var mods intact) + `peptide_ids` (with PIN-feature
   metavalues `COMET:*`, `MS:*`, … intact via `psm_metavalues`).
2. Existing `.pin` writer (`PercolatorInfile::store`) consumes
   `peptide_ids` as today.
3. Percolator runs; `.pout` parsed; q-values written back into
   `PeptideHit` scores (existing logic).
4. **New post-filter step**:
   - If `score:fdr < 1.0`: `IDFilter::filterHitsByScore(peptide_ids, cutoff)`.
     Drop empty `PeptideIdentification`s.
   - If `!keep_all_passing` (default): `IDFilter::keepBestPeptideHits(peptide_ids)`.
5. `FileHandler::storeIdentifications("scored.idparquet", …)`.

### Flow B — Backwards compatibility

The same code paths handle idXML, mzIdentML, and idparquet. Dispatch
is by extension. Mixed-extension pipelines work: e.g.
`--in merged.idparquet --out scored.idXML`.

### Flow C — `export_all_psms` semantics

`QPXFile::exportToArrow(..., export_all_psms)` already chooses "all
hits" vs "best hit only" at table-construction time. `PSMArrowIO::exportToParquet`
defaults `export_all_psms = true` — lossless round-trip is the point.
PercolatorAdapter's `keep_all_passing` flag operates on the in-memory
`peptide_ids` *before* the store call; the on-disk parquet faithfully
represents what is in memory. The two flags do not interact.

### Edge cases handled by existing infrastructure

- Multiple search engines per spectrum (post-IDMerger): `run_identifier`
  column distinguishes them; PSMSchema rows preserve order.
- Decoy PSMs: `is_decoy` column round-trips; `target_decoy` PeptideHit
  metavalue is restored on read.
- Charges / enzyme / `extra_features`: live in
  `SearchParamsSchema::metavalues` / `extra_features` field.

## Error handling

### Format layer (`PSMArrowIO`)

- **Output directory pre-flight.** Before writing, `exportToParquet`
  checks: parent directory exists; target path is either nonexistent
  or an existing directory. If target exists as a directory, the four
  canonical files are overwritten and any other sidecar files are
  left alone. If target exists as a regular file, log error, return `false`.
- **Per-table write failure.** Each `arrow::io::FileOutputStream::Open`
  / `parquet::arrow::WriteTable` call is checked. On the first failure:
  log which sub-file failed (`psms.parquet` / `proteins.parquet` /
  `protein_groups.parquet` / `search_params.parquet`), do not attempt
  the remaining writes, return `false`. Partial directories are not
  rolled back; the import path detects them.
- **Schema build failure.** `QPXFile::exportToArrow` and
  `ProteinIdentificationArrowIO::export*` already return `nullptr` on
  failure with their own logging. `PSMArrowIO` checks each return and
  propagates `false`.
- **Import: missing files.** All four files must be present. If any is
  missing, log which one, return `false`. Rationale: the directory is
  the format; an incomplete one is corrupt, not a degraded variant.
- **Import: schema mismatch.** Existing
  `ArrowSchemaValidation::validate` (subset mode) is already used
  inside `QPXFile` / `ProteinIdentificationArrowIO` /
  `ConsensusMapArrowIO`. Errors surface from there.
- **Import: empty PSM table.** Permitted; returns success with empty
  `peptide_ids`. Mirrors idXML behaviour.

### Dispatch layer (`FileHandler`)

`loadIdentifications` / `storeIdentifications` already throw
`Exception::UnableToCreateFile` / `Exception::FileNotFound` for I/O
issues and return early on type-detection failure. The `IDPARQUET`
branch follows the same convention: when `PSMArrowIO::*` returns
`false`, throw `Exception::UnableToCreateFile` (store) or
`Exception::ParseError` (load). Type detection is by extension
(`.idparquet`); mismatched/missing extensions fall through to the
existing `UNKNOWN` error path.

### Tool layer

- `PercolatorAdapter`: a filter that drops every PSM produces a
  warning (`OPENMS_LOG_WARN`) but does not fail the run — empty output
  is valid.
- `score:fdr` outside `[0.0, 1.0]` rejected at parse time
  (`setMinFloat_` / `setMaxFloat_`).
- `IDMerger`: no new error paths. Mixed-format inputs (e.g.,
  `.idXML` + `.idparquet`) already error today with "all must have the
  same type"; the message is unchanged.
- Search-engine adapters: no new error paths. Output dispatch via
  `FileHandler` inherits its error handling.

### Out of scope

- Best-effort partial reads.
- Format auto-promotion (e.g., reading a single `.parquet` file as a
  bundle).
- Retries on transient I/O.

## Testing

### Unit tests — `PSMArrowIO_test.cpp`

- *Round-trip, simple.* Construct minimal `protein_ids` + `peptide_ids`
  programmatically (one run, two PSMs, one protein hit). Export →
  import → assert structural equality on identifier, score, score
  type, sequence, charge, modifications, protein accessions,
  `target_decoy` metavalue.
- *Round-trip, full fidelity.* Load an existing idXML test fixture,
  export to idparquet, import back, compare to the original. Asserts:
  `extra_features` SearchParameters metavalue survives; PIN-feature
  PeptideHit metavalues (`COMET:deltaCn`, `MS:1002049`, …) survive;
  `enzyme` / `charges` / `fixed_modifications` round-trip.
- *Multi-run merged content.* Two `ProteinIdentification`s with distinct
  `run_identifier`s and PSMs distributed across both — round-trip
  preserves the run linkage.
- *Empty PSM table.* Export with empty `peptide_ids`, re-import yields
  empty `peptide_ids` and unchanged protein-side tables.
- *Error: missing sub-file.* Build a valid idparquet directory, delete
  `psms.parquet`, attempt import → returns `false`, error log mentions
  the missing file.
- *Error: target path is a regular file.* Pre-create a regular file at
  the target path, attempt export → returns `false`.

### Dispatch test (`FileHandler_test.cpp`)

Add one round-trip case using `loadIdentifications` /
`storeIdentifications` against an `.idparquet` path. Asserts the
dispatch picks the right branch by extension.

### TOPP tool tests

For each tool we touch, add at most one new test case that mirrors an
existing idXML case but uses idparquet. Reference outputs are
checked-in idparquet directories — small (low single-digit KB),
deterministic given the same input and `ParquetWriteConfig`:

- *PercolatorAdapter_idparquet_test* — input is a checked-in
  `.idparquet/` fixture (built once from the existing idXML fixture).
  Output is `.idparquet/`. Uses `score:fdr 0.01`. Compares against a
  reference idparquet directory file-by-file.
- *PercolatorAdapter_score_fdr_test* — option-only test on the existing
  idXML path. Feed pre-scored input where `score:fdr 0.01` should drop
  a known number of PSMs.
- *PercolatorAdapter_keep_all_passing_test* — confirms the flag retains
  multi-rank hits.
- *IDMerger_idparquet_test* — two `.idparquet/` inputs → one merged
  `.idparquet/`. Asserts row counts and that both `run_identifier`s
  appear in the proteins table.
- *CometAdapter_idparquet_test*, *SageAdapter_idparquet_test* — output
  an `.idparquet/`, smoke-test only (count PSMs, assert non-empty).
  Other engine adapters keep their existing idXML test; the parquet
  path is opt-in.

### Reference comparison

We do not diff parquet bytes (compression metadata is not deterministic
across Arrow patch versions). The test step imports the produced
directory back into memory and compares against the reference's
imported in-memory representation via the existing equality checks
used for idXML. A small CMake macro is added to
`src/tests/topp/CMakeLists.txt`:

```
add_idparquet_test(<tool> <name>)
```

It wraps the load-and-compare.

### Out of scope

- Arrow / Parquet library internals.
- Performance / file-size benchmarks.
- Cross-version `PSMSchema` migration (the schema is `@experimental`
  and has no v1/v2 story yet).

### CI

- The new tests need Arrow available. Existing parquet tests already
  require it — no new dependency.
- Reference idparquet fixtures are small and check in cleanly.
