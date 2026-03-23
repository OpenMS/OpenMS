# Replace timsrust with opentims for Bruker TimsTOF Reading

**Date:** 2026-03-23
**Status:** Approved
**Branch:** opentims

## Prerequisites

**License clarification (BLOCKING):** opentims repo-level LICENSE says MIT; individual C++ source headers say GPLv3. OpenMS is BSD-3-Clause — GPLv3 is incompatible. **Must confirm with maintainer (michalsta) that MIT applies to all files we consume before writing any code.** If GPLv3 is confirmed, fall back to Approach B (self-contained C++ reader without opentims dependency).

## Goal

Replace the timsrust Rust bridge with the opentims C++ library plus open-source calibration converters for reading Bruker TimsTOF `.d` files. Eliminate all Rust from the dependency chain. Maintain the existing public API (`BrukerTimsFile::load()`, `BrukerTimsFile::transform()`).

## Background

OpenMS currently reads Bruker TimsTOF `.d` directories via `timsrust_cpp_bridge`, which wraps the Rust `timsrust` library behind a C FFI. Platform-specific pre-built binaries are downloaded at CMake configure time. This works but introduces Rust as a foreign dependency in an otherwise pure C++ project.

[opentims](https://github.com/michalsta/opentims) is an open-source C++ library for reading TDF binary data. It provides TDF I/O (ZSTD decompression, memory-mapped frame reading) with a clean factory/strategy pattern for unit converters (TOF-to-m/z, scan-to-ion-mobility). However, opentims ships no open-source calibration — it relies on Bruker's proprietary `.so`/`.dll` for m/z and IM conversion.

timsrust's calibration math is trivially simple (~20 lines of arithmetic) and can be ported to C++ as custom converter factories registered into opentims's existing extension points.

**opentims reference version:** Pinned to commit on the `master` branch of `https://github.com/michalsta/opentims`. The exact commit hash must be determined at implementation time after license clarification. All API references in this spec are based on the opentims codebase as of 2026-03-23.

## Architecture

```
FileConverter (TOPP)
  -> FileHandler (dispatch on .d suffix)
    -> BrukerTimsFile (reader, unchanged public API)
      -> opentims TimsDataHandle / TimsFrame (TDF I/O)
      -> OpenSourceTof2MzConverter (calibration math)
      -> OpenSourceScan2ImConverter (calibration math)
      -> TdfMetadataReader (SQL queries for DIA/DDA metadata)
      -> analysis.tdf (SQLite) + analysis.tdf_bin (ZSTD-compressed frames)
```

## Components

### 1. opentims Integration (Build)

- `FetchContent_Declare(opentims ...)` in `cmake/cmake_findExternalLibs.cmake`, pinned to a specific commit hash
- Pull `src/opentims++/` only (no Python, R, Scala bindings)
- Strip opentims's bundled `sqlite/` directory: build opentims via explicit `add_library()` with a hand-maintained file list (not `add_subdirectory()`), excluding `sqlite/sqlite3.c`. Provide OpenMS's sqlite3 headers via `target_include_directories()` pointing at `src/openms/extern/SQLiteCpp/sqlite3/`. Link against OpenMS's existing `sqlite3` static library target.
- Keep opentims's bundled ZSTD (small, stateless, no conflict risk)
- Compile opentims sources directly into OpenMS
- `WITH_OPENTIMS` CMake option (default ON), replaces `WITH_TIMSRUST` entirely (no deprecation shim)
- Compile definition: `WITH_OPENTIMS=1`, header guards `#ifdef WITH_OPENTIMS`

### 2. Open-Source Calibration Converters

New files: `OpenTimsCalibration.h` / `OpenTimsCalibration.cpp` in `src/openms/{include,source}/OpenMS/FORMAT/`.

#### OpenSourceTof2MzConverter (extends opentims `Tof2MzConverter`)

Calibration model (linear in sqrt-space, matching timsrust):

```
intercept = sqrt(mz_min)
slope     = (sqrt(mz_max) - sqrt(mz_min)) / tof_max_index

convert:  mz = (intercept + slope * tof_index)^2
inverse:  tof = (sqrt(mz) - intercept) / slope
```

Parameters from `analysis.tdf` GlobalMetadata:
- `MzAcqRangeLower` -> mz_min
- `MzAcqRangeUpper` -> mz_max
- `DigitizerNumSamples` -> tof_max_index
- `AcquisitionSoftware`: if `"Bruker otofControl"`, expand mz range by +/-5.0 Da

The `frame_id` parameter is ignored (global calibration). This is a known accuracy trade-off vs Bruker's per-frame polynomial calibration. The global calibration means accuracy may drift over long acquisitions if the instrument's calibration shifts; the OLS recalibration (below) mitigates this.

#### OpenSourceScan2ImConverter (extends opentims `Scan2InvIonMobilityConverter`)

```
intercept = OneOverK0AcqRangeUpper  (im_max)
slope     = (im_min - im_max) / scan_max_index

convert:  im = intercept + slope * scan_index
```

Parameters from `analysis.tdf`:
- GlobalMetadata: `OneOverK0AcqRangeLower`, `OneOverK0AcqRangeUpper`
- Frames table: `MAX(NumScans)` -> scan_max_index

#### Factory Classes

```cpp
OpenSourceTof2MzConverterFactory : Tof2MzConverterFactory
  produce(TimsDataHandle&, pcs) -> reads tims_dir_path, constructs converter

OpenSourceScan2ImConverterFactory : Scan2InvIonMobilityConverterFactory
  produce(TimsDataHandle&, pcs) -> same pattern
```

Registered in `BrukerTimsFile::load()` before constructing any `TimsDataHandle`:
```cpp
DefaultTof2MzConverterFactory::setAsDefault<OpenSourceTof2MzConverterFactory>();
DefaultScan2InvIonMobilityConverterFactory::setAsDefault<OpenSourceScan2ImConverterFactory>();
```

The `pressure_compensation_strategy` parameter is part of opentims's `produce()` signature. Our factories accept it but do not act on it — the open-source linear model has no pressure correction. This is a silent no-op, not an error.

#### Optional OLS Recalibration

When `Config::calibrate = true`, a post-open pass refines the TOF→m/z slope and intercept:

1. For each DDA precursor, read `MonoisotopicMz` from the `Precursors` SQL table
2. Find the precursor's TOF index by locating the highest-intensity peak in the precursor's scan range within the relevant frame
3. Convert that TOF index to m/z using the current (unrefined) calibration
4. If `|converted_mz - MonoisotopicMz| < calibration_tolerance`, include the pair `(tof_index, sqrt(MonoisotopicMz))` in the regression set
5. With >= 2 pairs, compute OLS linear regression: `sqrt(mz) = intercept + slope * tof_index`
6. Update the converter's `intercept_` and `slope_` fields

Controlled by `Config::calibration_tolerance` (default 0.1 Da matching threshold). Requires >= 2 data points; if insufficient, the original calibration is kept with a warning.

### 3. TdfMetadataReader (Private Helper)

Private implementation detail inside `BrukerTimsFile.cpp`. Uses OpenMS's `SqliteConnector` or `SQLiteCpp::Database` to query `analysis.tdf`.

#### Calibration Metadata

```sql
SELECT Key, Value FROM GlobalMetadata
WHERE Key IN ('MzAcqRangeLower', 'MzAcqRangeUpper', 'DigitizerNumSamples',
              'OneOverK0AcqRangeLower', 'OneOverK0AcqRangeUpper',
              'AcquisitionSoftware')
```

#### Frame Metadata (for scan_max_index and frame counts)

```sql
SELECT MAX(NumScans) AS scan_max FROM Frames
SELECT COUNT(*) AS total_frames FROM Frames
SELECT COUNT(*) AS ms1_count FROM Frames WHERE MsMsType = 0
SELECT COUNT(*) AS ms2_count FROM Frames WHERE MsMsType != 0
```

These replace `tims_file_info()` for obtaining `num_frames`, `ms1.count`, `ms2.count`. The MS2 spectrum count for DDA is obtained from:
```sql
SELECT COUNT(DISTINCT Precursor) FROM PasefFrameMsMsInfo
```

#### DIA Mode Detection

```cpp
bool isDIA = sqliteHasNonEmptyTable(db, "DiaFrameMsMsWindows")
          || sqliteHasNonEmptyTable(db, "DiaFrameMsMsInfo");
```

#### DIA Window Query

```sql
-- Primary (newer firmware)
SELECT WindowGroup, ScanNumBegin, ScanNumEnd,
       IsolationMz, IsolationWidth, CollisionEnergy
FROM DiaFrameMsMsWindows

-- Fallback (older firmware)
SELECT Frame, ScanNumBegin, ScanNumEnd,
       IsolationMz, IsolationWidth, CollisionEnergy
FROM DiaFrameMsMsInfo
```

IM bounds per window derived by converting `ScanNumBegin`/`ScanNumEnd` via the scan-to-IM converter.

#### DDA Precursor Query

```sql
SELECT p.Id, p.MonoisotopicMz, p.Charge, p.Intensity,
       pfm.Frame, pfm.ScanNumBegin, pfm.ScanNumEnd,
       pfm.IsolationMz, pfm.IsolationWidth
FROM Precursors p
JOIN PasefFrameMsMsInfo pfm ON p.Id = pfm.Precursor
ORDER BY pfm.Frame, pfm.ScanNumBegin
```

### 4. BrukerTimsFile.cpp Rewrite

Same public API. Internal rewrite against opentims `TimsDataHandle` / `TimsFrame`.

#### API Mapping

| timsrust FFI | opentims |
|---|---|
| `tims_open()` | `TimsDataHandle(path)` |
| `tims_close()` | destructor (RAII) |
| `tims_get_frames_by_level(level)` | Filter `TimsDataHandle` frames by `TimsFrame::msms_type` field (0 = MS1, != 0 = MS2) |
| `tims_get_spectrum(idx)` | Reconstruct from frames + SQL (see DDA MS2 reconstruction below) |
| `tims_get_swath_windows()` | SQL query (see DIA window query) |
| `tims_convert_tof_to_mz_array()` | `TimsFrame::save_to_buffs()` returns m/z via registered converter |
| `tims_convert_scan_to_im_array()` | `TimsFrame::save_to_buffs()` returns inv_ion_mobility via registered converter |
| `tims_num_spectra()` | SQL: `SELECT COUNT(DISTINCT Precursor) FROM PasefFrameMsMsInfo` |
| `tims_file_info()` | SQL queries on Frames table (see Frame Metadata above) |
| `tims_free_*()` | Not needed (RAII) |

#### DDA MS2 Spectrum Reconstruction

This is the most significant new logic, replacing timsrust's pre-assembled `tims_get_spectrum()`:

1. Execute the DDA precursor query (Section 3) to get all precursor-to-frame mappings
2. Group results by precursor ID (a single precursor may span multiple PASEF frames)
3. For each precursor:
   a. For each associated frame: use the frame's `scan_offsets` array (length `num_scans + 1`) to identify peaks belonging to scans in `[ScanNumBegin, ScanNumEnd)`. `scan_offsets[s]` to `scan_offsets[s+1]` gives the peak index range for scan `s`. This filtering is done on raw integer scan indices (exact), not on converted IM values (avoids floating-point edge cases).
   b. Call `TimsFrame::save_to_buffs()` for the full frame to get m/z, intensity, and inv_ion_mobility arrays, then select only the peaks at the indices determined in step (a). (Alternatively, if opentims exposes per-scan extraction, use that directly.)
   c. Merge peaks across frames (concatenate, then sort by m/z)
   d. Compute scalar IM as intensity-weighted mean of peak IM values
   e. Build `MSSpectrum`:
      - Peak data: merged m/z + intensity arrays
      - Scalar `driftTime`: weighted mean IM (VSSC units)
      - RT: from the first contributing frame
      - MS level: 2
   f. Build `Precursor`:
      - m/z: `MonoisotopicMz` from SQL (selected ion m/z)
      - charge, intensity from SQL
      - isolation window: computed from `IsolationMz` and `IsolationWidth`
      - drift time: same as spectrum scalar IM

#### What stays the same

- MS1 frame loading: iterate frames, build CONCATENATED spectra with per-peak IM FloatDataArray
- DIA MS2 splitting: filter peaks by IM bounds per SWATH window
- IM-aware MS1 centroiding (`FrameCentroider` — pure OpenMS C++, unchanged)
- `IMSDataConsumer` streaming via `transform()`
- OpenMS metadata: nativeID, source files, CV terms, drift time units

#### Memory management simplification

The current `TimsDatasetHandle`, `TimsConfigHandle`, and `RAIICleanup` guards for C arrays are all replaced by opentims's C++ value types (`TimsDataHandle` owns everything, `TimsFrame` is copyable).

### 5. Config Changes

Removed fields (were default-off, no-ops):
- `smoothing_window` (timsrust applied during conversion; opentims has no equivalent)
- `centroiding_window` (same; OpenMS-side IM centroiding via `ms1_centroid_*` fields remains)

All other fields unchanged:
- `calibration_tolerance`, `calibrate` (drive OLS recalibration)
- `ms1_centroid_mz_ppm`, `ms1_centroid_im_pct` (OpenMS-side centroiding)
- `export_mode` (AUTO/SPECTRUM/FRAME)

### 6. TOPP Parameter Migration (User-Visible)

`FileConverter.cpp` registers a `timsrust:` TOPP parameter subsection (lines 180-204) with `#ifdef WITH_TIMSRUST`. This is a **user-visible parameter namespace change**:

- Rename subsection from `timsrust:` to `bruker:` (more generic, survives future backend changes)
- Remove `bruker:smoothing_window` and `bruker:centroiding_window` parameters
- Keep: `bruker:calibration_tolerance`, `bruker:calibrate`, `bruker:export_mode`, `bruker:ms1_centroid_mz_ppm`, `bruker:ms1_centroid_im_pct`
- Update `#ifdef WITH_TIMSRUST` to `#ifdef WITH_OPENTIMS`
- Update `getTimsConfig_()` helper (lines 226-241) accordingly

Other TOPP tools that accept `.d` input (e.g., SageAdapter) may also need parameter updates — check at implementation time.

## Known Trade-offs

### m/z Accuracy

The linear-in-sqrt calibration model is an approximation. AlphaTims documents typical errors < 0.02 Th but has observed outliers up to 6 Th. The optional OLS recalibration mitigates this. Bruker's proprietary SDK uses per-frame polynomial calibration with pressure compensation, which is more accurate. The global calibration means accuracy may degrade over long acquisitions if the instrument's calibration drifts frame-to-frame.

### Dropped Capabilities

- **TSF format**: timsrust supports TSF (non-TIMS TOF); opentims does not. Add later if needed as a separate reader.
- **In-bridge smoothing/centroiding**: Removed. Use OpenMS-side centroiding instead.
- **Pressure compensation**: Not replicated in the open-source converters. The `pressure_compensation_strategy` parameter is accepted by the factory interface but silently ignored.

## Testing Strategy

### Golden-File Regression

Before migration: generate reference TSV/mzML outputs from timsrust-based reader for all test datasets (DDA HeLa, DIA HeLa). After migration: compare opentims-based output with m/z tolerance (0.05 Th).

### Existing Unit Tests

`BrukerTimsFile_test.cpp` adapted to new internals. Same assertions on spectrum counts, MS levels, IM presence, mzML round-trip. Centroiding tests unchanged.

### New Tests

- Converter unit tests: known TOF indices -> expected m/z (hand-computed from metadata)
- SQL query tests: DIA window and precursor extraction against known `.d` files
- Factory registration: verify error without registration, success with registration
- Table fallback: verify `DiaFrameMsMsInfo` fallback when `DiaFrameMsMsWindows` absent
- Accuracy report: max/mean/95th percentile m/z deviation vs golden files

### Test Infrastructure Rename

CMake test macros rename from `TIMSRUST_*` to `OPENTIMS_*`:
- `ENABLE_TIMSRUST_TESTS` -> `ENABLE_OPENTIMS_TESTS`
- `TIMSRUST_DDA_TEST_DATA` -> `OPENTIMS_DDA_TEST_DATA`
- `TIMSRUST_DIA_TEST_DATA` -> `OPENTIMS_DIA_TEST_DATA`
- `TIMSRUST_TEST_FASTA` -> `OPENTIMS_TEST_FASTA`
- `TIMSRUST_DDA_SYMLINK` -> `OPENTIMS_DDA_SYMLINK` (in `src/tests/topp/CMakeLists.txt`)

TOPP test parameters update from `timsrust:*` to `bruker:*` in test command lines.

## File Changes

| File | Action |
|---|---|
| `CMakeLists.txt` (top-level) | Replace `WITH_TIMSRUST`/`ENABLE_TIMSRUST_TESTS` options and test data fetch with `WITH_OPENTIMS`/`ENABLE_OPENTIMS_TESTS` |
| `cmake/cmake_findExternalLibs.cmake` | Replace timsrust FetchContent with opentims FetchContent |
| `cmake/OpenMSConfig.cmake.in` | Replace `set(WITH_TIMSRUST ...)` export with `set(WITH_OPENTIMS ...)` |
| `src/openms/CMakeLists.txt` | Update source lists, link targets, compile definitions (`WITH_TIMSRUST` -> `WITH_OPENTIMS`) |
| `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` | Update includes, remove smoothing/centroiding config fields, update `#ifdef` |
| `src/openms/source/FORMAT/BrukerTimsFile.cpp` | Full rewrite against opentims API |
| `src/openms/include/OpenMS/FORMAT/OpenTimsCalibration.h` | **NEW**: converter + factory class declarations |
| `src/openms/source/FORMAT/OpenTimsCalibration.cpp` | **NEW**: converter implementation + calibration math |
| `src/openms/source/FORMAT/FileHandler.cpp` | Update `#ifdef` from `WITH_TIMSRUST` to `WITH_OPENTIMS` |
| `src/topp/FileConverter.cpp` | Rename `timsrust:` subsection to `bruker:`, remove smoothing/centroiding params, update `#ifdef` and `getTimsConfig_()` |
| `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp` | Adapt to API changes, add converter/SQL tests |
| `src/tests/topp/CMakeLists.txt` | Rename `TIMSRUST_*` variables/symlinks, update test parameter names |
