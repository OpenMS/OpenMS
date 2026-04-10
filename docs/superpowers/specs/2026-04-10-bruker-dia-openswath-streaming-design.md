# Bruker DIA → OpenSwath Consumer-Based Streaming

## Goal

Replace the current Bruker DIA → OpenSwath loading path (which materializes the entire dataset as an MSExperiment before partitioning) with a consumer-based streaming path that feeds spectra one-at-a-time to a `FullSwathFileConsumer`. This enables disk-cached loading for large diaPASEF datasets and eliminates the 2× memory duplication of the current approach.

## Problem

The current path in `SwathFile::loadBrukerTdf()`:

1. `BrukerTimsFile::load()` → builds full MSExperiment (all ~8.5K spectra, ~50M peaks, ~8 GB RAM)
2. `countScansInSwath_()` → re-scans all spectra to discover SWATH windows (redundant: windows are already known from SQL)
3. `RegularSwathFileConsumer` → partitions spectra into per-window PeakMaps (second copy in RAM)
4. Returns `SpectrumAccessOpenMS` per SwathMap (pointing to the per-window PeakMaps)

For a production diaPASEF run (60-min gradient, ~100K frames), this would require tens of GB of RAM. The mzML path solves this with `CachedSwathFileConsumer` (stream to `.mzML.cached` on disk, keep only metadata in memory), but Bruker loading doesn't support it.

## Design

### New methods on `BrukerTimsFile`

Two new public methods split the work into metadata (cheap) and streaming (heavy):

**Step 1 — `readDIAMetadata()`**: Read SWATH boundaries and spectrum counts from SQL without touching peak data. The caller uses this to construct the appropriate consumer.

```cpp
/// Metadata needed to construct a FullSwathFileConsumer for DIA streaming.
struct DIAStreamingMetadata
{
    std::vector<OpenSwath::SwathMap> boundaries;  // one per DIAWindow
    int nr_ms1_spectra = 0;
    std::vector<int> nr_ms2_spectra;              // per-window spectrum counts
};

/// Read DIA metadata from a Bruker .d directory (SQL only, no peak data).
/// Also populates exp_settings with instrument/source metadata.
DIAStreamingMetadata readDIAMetadata(
    const String& path,
    ExperimentalSettings& exp_settings,
    const Config& config = Config());
```

**Step 2 — `loadDIAStreaming()`**: Stream all spectra to a pre-constructed consumer. The consumer must have been initialized with the boundaries from step 1.

```cpp
/// Stream DIA spectra to a consumer one-at-a-time.
/// MS2 spectra are always raw (no aggregation/denoising/centroiding).
/// MS1 centroiding respects the Config settings.
void loadDIAStreaming(
    const String& path,
    FullSwathFileConsumer& consumer,
    const Config& config = Config());
```

**MS2 is always raw**: The `dia_ms2_n_neighbors`, `dia_ms2_min_support`, and `dia_ms2_centroid` config parameters are ignored. OpenSwath operates on raw frame-level data for chromatogram extraction. Aggregation/denoising is only relevant for the mzML export path (FileConverter, PeakPickerIM).

**MS1 centroiding**: Respects `ms1_centroid_mz_ppm` and `ms1_centroid_im_pct` from Config, same as the existing path.

### Modified: `SwathFile::loadBrukerTdf()`

Extended signature to accept `readoptions` and `tmp` (cache directory), mirroring `loadMzML()`:

```cpp
std::vector<OpenSwath::SwathMap> loadBrukerTdf(
    const String& file,
    const String& tmp,
    std::shared_ptr<ExperimentalSettings>& exp_meta,
    const String& readoptions);
```

Implementation:
```cpp
BrukerTimsFile bruker;
bruker.setLogType(this->getLogType());

// Step 1: metadata (SQL only, no peak data)
ExperimentalSettings settings;
auto meta = bruker.readDIAMetadata(file, settings);
auto exp_meta_ptr = std::make_shared<PeakMap>(settings);
exp_meta = exp_meta_ptr;

// Step 2: construct consumer based on readoptions
std::shared_ptr<FullSwathFileConsumer> consumer;
if (readoptions == "normal")
    consumer = std::make_shared<RegularSwathFileConsumer>(meta.boundaries);
else if (readoptions == "cache")
    consumer = std::make_shared<CachedSwathFileConsumer>(
        meta.boundaries, tmp, File::getUniqueName(),
        meta.nr_ms1_spectra, meta.nr_ms2_spectra);
else
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unknown or unsupported readoption '" + readoptions + "' for Bruker TDF loading");
consumer->setExperimentalSettings(settings);

// Step 3: stream spectra
bruker.loadDIAStreaming(file, *consumer);

// Step 4: finalize
std::vector<OpenSwath::SwathMap> swath_maps;
consumer->retrieveSwathMaps(swath_maps);
return swath_maps;
```

Note: `readDIAMetadata()` and `loadDIAStreaming()` each open their own `TimsDataHandle` internally. Opening the handle is cheap (calibration setup + SQL metadata); the expensive part is reading peak data. If profiling shows this matters, a future optimization could share the handle between the two calls via an internal cache or a combined method.

The old signature `loadBrukerTdf(file, exp_meta)` becomes a convenience wrapper that calls the new one with `readoptions="normal"`.

### Modified: `OpenSwathBase::loadSwathFiles_()`

Pass `readoptions` and `tmp` through to the Bruker loading path (currently hardcoded to in-memory):

```cpp
// Before:
auto maps = swath_file.loadBrukerTdf(f, exp_meta);

// After:
auto maps = swath_file.loadBrukerTdf(f, tmp, exp_meta, readoptions);
```

## Data Flow

### Current (all in memory)

```
.d → BrukerTimsFile::load()
   → MSExperiment (8.5K spectra, 50M peaks, ~8 GB)
   → countScansInSwath_() [redundant scan]
   → RegularSwathFileConsumer [second copy]
   → SpectrumAccessOpenMS per SwathMap
```

### New: readoptions="cache" (streaming to disk)

```
.d → readDIAWindows() + readFrameToWindowGroupMapping()      [SQL, ~1ms]
   → SwathMap boundaries from DIAWindow structs               [no discovery pass]
   → CachedSwathFileConsumer(boundaries, tmp_dir)
   → iterate per-WindowGroup/frame:
       build 1 MSSpectrum → consumer.consumeSpectrum()        [O(1 spectrum) memory]
         → MSDataCachedConsumer writes to .mzML.cached        [data on disk]
         → metadata-only spectrum kept in PeakMap              [lightweight]
   → consumer.retrieveSwathMaps()
   → SpectrumAccessOpenMSCached per SwathMap                  [on-demand disk reads]
```

### New: readoptions="normal" (in memory, but single pass)

```
.d → readDIAWindows() + readFrameToWindowGroupMapping()
   → SwathMap boundaries
   → RegularSwathFileConsumer(boundaries)
   → iterate per-WindowGroup/frame:
       build 1 MSSpectrum → consumer.consumeSpectrum()
         → spectrum appended to per-window PeakMap
   → SpectrumAccessOpenMS per SwathMap
```

Still in-memory, but avoids the intermediate full MSExperiment and the countScansInSwath_ discovery pass. One copy of data, not two.

## Implementation Details

### Building SwathMap boundaries from DIAWindow structs

The consumer needs boundaries before consuming spectra. We build them directly from SQL:

```cpp
// In loadDIAToSwathMaps():
auto windows = readDIAWindows(db, *handle.scan2inv_ion_mobility_converter);
std::vector<OpenSwath::SwathMap> boundaries;
for (const auto& w : windows)
{
    boundaries.emplace_back(
        w.mz_center - w.mz_width / 2.0,  // lower
        w.mz_center + w.mz_width / 2.0,  // upper
        w.mz_center,                       // center
        w.im_lower,                        // imLower
        w.im_upper,                        // imUpper
        false);                            // ms1 = false
}
```

### Spectrum count pre-computation

`CachedSwathFileConsumer` accepts `nr_ms1_spectra` and `nr_ms2_spectra` (per-window counts) for pre-sizing cache files. Both are computable from SQL without reading peak data:

```cpp
// MS1: count frames with msms_type == 0
int nr_ms1_spectra = 0;
for (uint32_t fid = handle.min_frame_id(); fid <= handle.max_frame_id(); ++fid)
    if (handle.has_frame(fid) && handle.get_frame(fid).msms_type == 0)
        ++nr_ms1_spectra;

// MS2: for each DIAWindow, count frames belonging to its WindowGroup
auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
std::vector<int> nr_ms2_spectra;
nr_ms2_spectra.reserve(windows.size());
for (const auto& w : windows)
{
    auto it = group_to_frames.find(w.window_group);
    nr_ms2_spectra.push_back(it != group_to_frames.end()
        ? static_cast<int>(it->second.size()) : 0);
}
```

### Per-WindowGroup iteration (MS2 raw path)

Identical to the corrected raw path in `loadDIA_()`, but calls `consumer.consumeSpectrum(spec)` instead of `exp.addSpectrum(std::move(spec))`. The spectrum is consumed (data cleared by CachedSwathFileConsumer after writing to disk), then goes out of scope.

### MS1 streaming

Same as existing `loadDIA_()` MS1 loop: iterate MS1 frames, optionally centroid, build spectrum, feed to consumer as MS level 1.

### ExperimentalSettings

The consumer needs `ExperimentalSettings` via `setExperimentalSettings()`. The existing `BrukerTimsFile::load()` sets a `SourceFile` entry on the MSExperiment (source file name, path, type, nativeID format). It does NOT currently read instrument model or contact info from SQL tables.

Factor the `SourceFile` setup into a private `loadExperimentalSettings_(TimsDataHandle&, ExperimentalSettings&)` helper used by both `load()` and `readDIAMetadata()`. No peak data is read — only the source file metadata already set by the existing code.

### Thread safety

`lightClone()` on the resulting `SpectrumAccessOpenMSCached` creates independent file handles for parallel extraction threads. This works the same as for mzML-cached files.

## Edge Cases

- **DDA files**: `readDIAMetadata()` calls `isDIA_()` internally. If the file is DDA (no DIA windows), it throws `Exception::InvalidParameter` — callers should only pass DIA files to this path. `SwathFile::loadBrukerTdf()` can check file type before calling.
- **Empty DIA windows**: If `readDIAWindows()` returns empty, `readDIAMetadata()` returns empty boundaries. The consumer produces no SwathMaps. Same behavior as existing code (BrukerTimsFile.cpp line 1578-1582 logs a warning).
- **Old TDF format**: `readDIAWindows()` and `readFrameToWindowGroupMapping()` already handle the old format (no `DiaFrameMsMsWindows` table) via geometry deduplication. No changes needed.
- **No MS1 frames**: `nr_ms1_spectra` is 0. The consumer handles `ms1_map_ == nullptr` gracefully (SwathFileConsumer.h line 130).
- **Spectrum ordering**: The streaming path emits all MS1 first, then all MS2 by WindowGroup — unlike mzML which interleaves. The consumer routes correctly regardless of order; this only affects cache file write patterns (all MS1 data is contiguous in the MS1 cache file).

## What doesn't change

- `loadDIA_()` for the mzML export path (FileConverter, PeakPickerIM) — unchanged, still supports aggregation/denoising/centroiding
- OpenSwathWorkflow's downstream processing — it already handles `SpectrumAccessOpenMSCached` from the mzML path
- `load_into_memory` option — already wraps any `ISpectrumAccess` in `SpectrumAccessOpenMSInMemory`; works transparently with the new cached Bruker path

## Files to Modify

| File | Change |
|------|--------|
| `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` | Add `DIAStreamingMetadata` struct, `readDIAMetadata()`, `loadDIAStreaming()` declarations |
| `src/openms/source/FORMAT/BrukerTimsFile.cpp` | Implement `readDIAMetadata()` and `loadDIAStreaming()`; factor out `loadExperimentalSettings_()` |
| `src/openms/include/OpenMS/FORMAT/SwathFile.h` | Extend `loadBrukerTdf()` signature with `tmp`, `readoptions` |
| `src/openms/source/FORMAT/SwathFile.cpp` | Implement new `loadBrukerTdf()` with consumer selection; old signature becomes wrapper |
| `src/openms/source/APPLICATIONS/OpenSwathBase.cpp` | Pass `readoptions` and `tmp` to Bruker path (both call sites at lines 135 and 181) |
| `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp` | Add test sections for `readDIAMetadata()` and `loadDIAStreaming()` |

## Testing

- Unit test: Load DIA test data via `loadDIAToSwathMaps()` with `RegularSwathFileConsumer`, verify SwathMap count and spectrum counts match the existing `loadBrukerTdf()` output
- Unit test: Load with `CachedSwathFileConsumer`, verify `.mzML.cached` files are created and `SpectrumAccessOpenMSCached` returns correct spectra
- Verify no MS2 aggregation/denoising occurs regardless of Config settings
- Verify MS1 centroiding works when enabled
