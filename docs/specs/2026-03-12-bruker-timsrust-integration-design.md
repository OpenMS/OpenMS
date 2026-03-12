# Bruker TimsTOF Integration via timsrust_cpp_bridge

**Date:** 2026-03-12
**Status:** Approved

## Overview

Add native support for reading Bruker TimsTOF `.d` files into OpenMS via the open-source [timsrust_cpp_bridge](https://github.com/OpenMS/timsrust_cpp_bridge) library. This enables direct conversion of DDA-PASEF and DIA-PASEF data to mzML, with support for raw frame-level 4D access, integrated into `FileConverter` and the `FileHandler` dispatch system. The integration is cross-platform, vendor-DLL-independent, and guarded by a `WITH_TIMSRUST` CMake flag (default ON).

## Goals

- Load Bruker `.d` directories (TDF format) into `MSExperiment` with correct IM data representation for DDA and DIA workflows
- Expose timsrust processing parameters (smoothing, centroiding, calibration) via `FileConverter`
- Support streaming conversion for large datasets via `IMSDataConsumer`
- Provide raw frame-level 4D export mode for custom analysis
- Integrate seamlessly with existing IM-aware tools (`PeakPickerIM`, `FeatureFinderMetabo`, `OpenSwathWorkflow`)

## Non-Goals

- No vendor DLL fallback path. The integration uses timsrust exclusively.
- No new TOPP tool. All functionality lives in `FileConverter` behind `#ifdef WITH_TIMSRUST` guards.
- No miniTDF support in the initial implementation (can be added later since the bridge supports it).

## Architecture

```
┌─────────────────────────────────────────────────┐
│ FileConverter (TOPP)                            │
│  #ifdef WITH_TIMSRUST                           │
│  ├─ timsrust_* parameters                       │
│  └─ auto/spectrum/frame export mode             │
├─────────────────────────────────────────────────┤
│ FileHandler                                     │
│  ├─ detectType(.d dir) → BRUKER_TDF            │
│  └─ loadExperiment() → delegates to ↓           │
├─────────────────────────────────────────────────┤
│ BrukerTimsFile                                  │
│  ├─ load(path, MSExperiment, config)            │
│  ├─ transform(path, Interfaces::IMSDataConsumer*)│
│  ├─ loadDDA()   ─┐                             │
│  ├─ loadDIA()    ├─ internal mapping logic      │
│  └─ loadFrames() ┘                             │
├─────────────────────────────────────────────────┤
│ timsrust_cpp_bridge (C ABI)                     │
│  ├─ tims_open[_with_config]()                   │
│  ├─ tims_get_spectrum()      ← DDA MS2         │
│  ├─ tims_get_frame()         ← MS1 + DIA MS2   │
│  ├─ tims_get_swath_windows() ← DIA detection    │
│  ├─ tims_convert_tof_to_mz_array()             │
│  ├─ tims_convert_scan_to_im_array()            │
│  └─ tims_file_info()         ← metadata         │
└─────────────────────────────────────────────────┘
```

## CMake Integration

### Option Flags

Top-level `CMakeLists.txt`:
```cmake
option(WITH_TIMSRUST "Build with Bruker TimsTOF .d file support via timsrust" ON)
option(ENABLE_TIMSRUST_TESTS "Download Bruker test data for timsrust integration tests" OFF)
```

### Library Acquisition

In `cmake/cmake_findExternalLibs.cmake`, use `FetchContent` to download the platform-specific release archive from GitHub at configure time. Note: this introduces `FetchContent` to the core OpenMS build for the first time (previously only used in `src/pyOpenMS/CMakeLists.txt` for nanobind).

- Version pinned via a `TIMSRUST_VERSION` CMake variable for easy bumping
- Platform detection selects the correct archive (linux-x86_64, linux-aarch64, macos-arm64, windows-x86_64)
- After extraction, `find_package(timsrust_cpp_bridge REQUIRED)` locates the CMake config files
- Link `timsrust_cpp_bridge::timsrust_cpp_bridge` as a private dependency in `src/openms/CMakeLists.txt`

### Compile-Time Flag

In `src/openms/include/OpenMS/config.h.in`:
```c
#cmakedefine WITH_TIMSRUST 1
```

### Test Data

When `ENABLE_TIMSRUST_TESTS=ON`, FetchContent downloads `DDA_HeLa_50ng_5_6min.d.zip` and `DIA_HeLa_50ng_5_6min.d.zip` from the `test-data-v1` GitHub release.

## File Type Registration

### FileTypes.h

Add `BRUKER_TDF` to the `FileTypes::Type` enum.

### FileTypes.cpp

Add `TypeNameBinding` entry in the `type_with_annotation__` array:
- Extension name: `"d"` (matches the `.d` directory suffix)
- Description: `"Bruker TDF"`
- `FileProperties` flags: `PROVIDES_EXPERIMENT`, `READABLE`

Add `typeToMZML` mapping for `BRUKER_TDF` → `"Bruker TDF format"` (matching the PSI-MS term in `share/OpenMS/CV/psi-ms.obo`).

### FileHandler.cpp

Detection logic — the existing flow is file-oriented; directories need special handling:

1. `getTypeByFileName()` already handles `.d` via the `TypeNameBinding` with extension `"d"` → returns `BRUKER_TDF` automatically
2. In `getType()`: once `getTypeByFileName()` returns `BRUKER_TDF`, short-circuit before calling `getTypeByContent()` (which opens the path as a file and would fail for directories). Optionally validate that the directory contains `analysis.tdf` or `analysis.tdf_bin` to distinguish from other `.d` directories.
3. In `loadExperiment()`: skip `computeFileHash()` for directory inputs (it opens the path as a file and returns empty string for directories). Handle `File::basename()` for paths with trailing slashes (returns empty string). Set source file metadata without hash.

Dispatch: `loadExperiment()` delegates to `BrukerTimsFile::load()` when type is `BRUKER_TDF`.

## Core Reader: BrukerTimsFile

### Header: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`
### Source: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

Both files conditionally compiled via `sources.cmake` when `WITH_TIMSRUST` is ON.

### Public API

```cpp
class OPENMS_DLLAPI BrukerTimsFile
{
public:
  struct Config
  {
    uint32_t smoothing_window = 0;         // 0 = timsrust default; bin count
    uint32_t centroiding_window = 0;       // 0 = timsrust default; bin count
    double calibration_tolerance = 0.0;    // 0 = timsrust default; m/z tolerance
    bool calibrate = false;                // off by default (timsrust 0.4.2 caveat)
    enum ExportMode { AUTO, SPECTRUM, FRAME };
    ExportMode export_mode = AUTO;
  };

  /// Load entire .d directory into MSExperiment
  void load(const String& path, MSExperiment& exp, const Config& config = {});

  /// Streaming: read .d and feed spectra to consumer
  void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config = {});
};
```

### DDA Mapping (`loadDDA`)

**MS1 frames:**
1. `tims_get_frames_by_level(handle, 1)` → fetch all MS1 frames
2. For each frame: `tims_convert_tof_to_mz_array()` + `tims_convert_scan_to_im_array()`
3. Build one `MSSpectrum` per frame in **CONCATENATED** format:
   - Peak data: m/z (from converted TOF, double precision) + intensity (uint32 → double)
   - IM `FloatDataArray` configured via `IMDataConverter::setIMUnit(array, DriftTimeUnit::VSSC)`, which sets the array name to `"raw inverse reduced ion mobility array"` with CV term `MS:1003008`. This is required: `MSSpectrum::containsIMData()` delegates to `IMDataConverter::getIMUnit()` which inspects array names, not CV accessions. Using the wrong name (e.g., `"Ion Mobility"`) would cause misclassification as MILLISECOND.
   - `DriftTimeUnit::VSSC`
   - RT from `tims_frame.rt_seconds`
   - MS level = 1

**Scan-offset to IM mapping** (core algorithm for frame-level path):
```
for scan_idx in 0..num_scans:
    start = scan_offsets[scan_idx]
    end   = scan_offsets[scan_idx + 1]
    im    = tims_convert_scan_to_im(handle, scan_idx)
    for peak in start..end:
        // peak has tof_indices[peak] → m/z, intensities[peak], and im
```

**MS2 spectra:**
1. Iterate `tims_get_spectrum(handle, i)` for i in 0..`tims_num_spectra()`
2. One `MSSpectrum` per precursor:
   - Peak data: m/z + intensity arrays from `tims_spectrum` (note: bridge returns `float*` for both — must be widened to `double` for `Peak1D`)
   - Scalar `driftTime` = `tims_spectrum.im` (1/K₀, VSSC)
   - `Precursor` populated with:
     - m/z = `tims_spectrum.precursor_mz`
     - charge = `tims_spectrum.charge`
     - intensity = `tims_spectrum.precursor_intensity`
     - drift time = `tims_spectrum.im`
     - isolation window = `tims_spectrum.isolation_mz` ± `tims_spectrum.isolation_width / 2`
   - RT from `tims_spectrum.rt_seconds`
   - MS level = 2

### DIA Mapping (`loadDIA`)

**MS1 frames:** Same as DDA.

**MS2 frames:**
1. `tims_get_swath_windows(handle)` → get all DIA isolation windows with m/z and IM bounds
2. `tims_get_frames_by_level(handle, 2)` → fetch all MS2 frames
3. For each frame: convert TOF→m/z and scan→IM
4. Split peaks by SWATH window: assign each peak to the window whose IM bounds contain the peak's IM value
5. One `MSSpectrum` per window per frame in **CONCATENATED** format:
   - Per-peak IM in `FloatDataArray` configured via `IMDataConverter::setIMUnit(..., DriftTimeUnit::VSSC)` (sets name `"raw inverse reduced ion mobility array"`, CV `MS:1003008`)
   - Isolation window metadata from `tims_swath_window` (mz_lower, mz_upper, mz_center)
   - RT from frame
   - MS level = 2

This preserves within-window IM resolution needed by `OpenSwathWorkflow` for 2D extraction and `IonMobilityScoring`.

### Raw Frame Export (`loadFrames`)

- All frames via `tims_get_frames_by_level()` for each MS level
- Convert TOF→m/z and scan→IM
- One CONCATENATED `MSSpectrum` per frame with per-peak IM
- No precursor grouping or SWATH window splitting
- Pure 4D dump (RT × IM × m/z × intensity)

### DDA vs DIA Auto-Detection

```cpp
tims_get_swath_windows(handle, &count, &windows);
if (count > 0) → DIA path
else → DDA path
```

### Streaming (`transform`)

1. Open dataset, call `tims_file_info()` for counts
2. Compute expected spectrum count based on export mode:
   - SPECTRUM/DDA: `tims_file_info().num_frames` (MS1) + `tims_num_spectra()` (MS2)
   - DIA: `tims_file_info().num_frames` (MS1) + num_ms2_frames × num_swath_windows (MS2)
   - FRAME: `tims_file_info().num_frames` (all levels)
3. `consumer->setExpectedSize(computed_count, 0)`
4. `consumer->setExperimentalSettings(...)` — note: `tims_file_info()` only provides counts and ranges (RT/m/z/IM/intensity), not instrument identity. Instrument metadata must be populated from what's available (RT range, m/z range, etc.); instrument name/model are not available from the bridge.
5. Iterate frames/spectra in RT order, convert, call `consumer->consumeSpectrum()` for each
6. Enables constant-memory mzML conversion via `PlainMSDataWritingConsumer`

### Export Mode Semantics

| Mode | DDA dataset | DIA dataset |
|------|-------------|-------------|
| AUTO | DDA path (spectrum reader for MS2) | DIA path (frame reader + SWATH split for MS2) |
| SPECTRUM | DDA-style: spectrum reader for MS2, no window splitting | DDA-style even for DIA: spectrum reader for MS2, no window splitting (loses within-window IM resolution) |
| FRAME | Raw frame dump, no precursor grouping | Raw frame dump, no SWATH window splitting |

### Spectrum Ordering

Output spectra are interleaved by RT across MS levels, matching the conventional mzML output from TimsTOF instruments (as produced by MSConvert). MS1 and MS2 spectra from the same RT region appear adjacent, sorted by RT.

### Memory Management

The bridge has strict ownership semantics:
- `tims_get_spectrum()` returns pointers into internal buffers that become invalid on the next call — all data must be copied into `MSSpectrum` before the next bridge call
- `tims_get_frames_by_level()`, `tims_get_swath_windows()` return caller-owned arrays that require explicit `tims_free_*` calls

The implementation must use RAII wrappers (scope guards or `unique_ptr` with custom deleters) for `tims_dataset*`, `tims_frame*` arrays, and `tims_swath_window*` arrays to prevent leaks on exception paths.

### Common Properties

- All IM values in 1/K₀ (VSSC) units — matching TimsTOF native representation
- `DriftTimeUnit::VSSC` set on all spectra
- IM `FloatDataArray` entries must be configured via `IMDataConverter::setIMUnit(array, DriftTimeUnit::VSSC)` — sets name `"raw inverse reduced ion mobility array"` with CV `MS:1003008`
- Spectra sorted by RT after loading, interleaved across MS levels
- `ExperimentalSettings` populated from `tims_file_info()` where available (RT/m/z/IM ranges, peak counts). Instrument name/model not available from bridge — left empty.
- Error handling: check `timsffi_status` return codes, throw `Exception::FileNotReadable` or `Exception::ParseError` with message from `tims_get_last_error()`

## FileConverter Integration

In `src/topp/FileConverter.cpp`, under `#ifdef WITH_TIMSRUST`:

### Additional Parameters

```
-timsrust_smoothing_window <int>                    (default: 0 = library default; bin count)
-timsrust_centroiding_window <int>                  (default: 0 = library default; bin count)
-timsrust_calibration_tolerance <float>             (default: 0 = library default; m/z tolerance)
-timsrust_calibrate <string: "true","false">        (default: "false"; modeled as string toggle per TOPP conventions)
-timsrust_export_mode <string: "auto","spectrum","frame">  (default: "auto")
```

Note: `timsrust_calibrate` uses a `"true"/"false"` string restriction rather than a raw bool, following the existing TOPP parameter pattern (see e.g., `force` parameter).

### Input Format Registration

Add `"d"` to FileConverter's input format list (currently at line ~152). This is required for TOPP parameter validation.

### Processing

- When input type is `BRUKER_TDF`, build `BrukerTimsFile::Config` from parameters
- Low-memory mode (`-process_lowmemory`): currently only allows mzML/mzXML input (throws `IllegalArgument` otherwise). Must extend the low-memory branch to accept `BRUKER_TDF` as input, routing to `BrukerTimsFile::transform()` with a writing consumer.
- Normal mode: use `BrukerTimsFile::load()` then write

## Vendor DLL Independence

This integration uses timsrust exclusively. It does **not** use or depend on `timsdata.dll` (the proprietary Bruker SDK found in `THIRDPARTY/pwiz-bin/`). The TOF→m/z and scan→IM conversions are performed by timsrust's open-source implementation reading the TDF SQLite calibration data.

No vendor DLL fallback path is provided. If numerical accuracy against the vendor SDK becomes a concern, that is addressed upstream in timsrust, not by branching the OpenMS reader.

## Testing

### Unit Tests (always compiled when `WITH_TIMSRUST=ON`)

File: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`

- Error handling: invalid path, missing TDF files, bridge error codes
- DDA/DIA auto-detection logic
- SWATH window parsing
- `FileTypes::BRUKER_TDF` detection in `FileHandler`
- CMake linkage verification

### Integration Tests (gated behind `ENABLE_TIMSRUST_TESTS=ON`)

Download real HeLa datasets from `test-data-v1` release:

- **DDA test**: Load `DDA_HeLa_50ng_5_6min.d` → verify spectrum count, RT range, MS levels, precursor metadata, MS1 CONCATENATED with IM FloatDataArray, MS2 with scalar drift times
- **DIA test**: Load `DIA_HeLa_50ng_5_6min.d` → verify SWATH window count and bounds, MS2 per-peak IM arrays, isolation window metadata
- **Round-trip test**: Load `.d` → write mzML → reload → compare spectrum count, RT range, m/z range, IM range
- **Streaming test**: Compare `load()` vs `transform()` output for consistency
- **Frame-level test**: Load in frame mode → verify all frames present, per-peak IM arrays, no precursor grouping

## Files Changed

| Action | File |
|--------|------|
| Create | `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` |
| Create | `src/openms/source/FORMAT/BrukerTimsFile.cpp` |
| Create | `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp` |
| Modify | `CMakeLists.txt` — add `WITH_TIMSRUST` and `ENABLE_TIMSRUST_TESTS` options |
| Modify | `cmake/cmake_findExternalLibs.cmake` — FetchContent + find_package |
| Modify | `src/openms/include/OpenMS/config.h.in` — add `#cmakedefine WITH_TIMSRUST 1` |
| Modify | `src/openms/CMakeLists.txt` — conditional linking + public compile definition |
| Modify | `cmake/OpenMSConfig.cmake.in` — export `WITH_TIMSRUST` for downstream consumers |
| Modify | `src/openms/include/OpenMS/FORMAT/FileTypes.h` — add `BRUKER_TDF` enum value |
| Modify | `src/openms/source/FORMAT/FileTypes.cpp` — add `TypeNameBinding`, `typeToMZML` mapping |
| Modify | `src/openms/source/FORMAT/FileHandler.cpp` — directory-aware detection + dispatch |
| Modify | `src/openms/source/FORMAT/sources.cmake` — conditional `BrukerTimsFile.cpp` |
| Modify | `src/openms/include/OpenMS/FORMAT/sources.cmake` — conditional `BrukerTimsFile.h` |
| Modify | `src/topp/FileConverter.cpp` — timsrust parameters + BRUKER_TDF input |
| Modify | `src/tests/class_tests/openms/executables.cmake` — add `BrukerTimsFile_test` |
| Modify | `src/tests/class_tests/openms/source/FileTypes_test.cpp` — update type count |
