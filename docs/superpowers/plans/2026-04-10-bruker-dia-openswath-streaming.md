# Bruker DIA → OpenSwath Consumer-Based Streaming Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stream Bruker diaPASEF spectra directly to OpenSwath's `FullSwathFileConsumer` without an intermediate MSExperiment, enabling disk-cached loading for large datasets.

**Architecture:** Two new public methods on `BrukerTimsFile` — `readDIAMetadata()` (SQL-only, returns SwathMap boundaries + counts) and `loadDIAStreaming()` (streams spectra one-at-a-time to a consumer). `SwathFile::loadBrukerTdf()` is extended with `readoptions`/`tmp` parameters and wired through `OpenSwathBase`.

**Tech Stack:** C++ (OpenMS), SQLite (TDF metadata), opentims++ (frame decompression)

**Spec:** `docs/superpowers/specs/2026-04-10-bruker-dia-openswath-streaming-design.md`

---

### Task 1: Factor out `loadExperimentalSettings_()` helper

Both `readDIAMetadata()` and the existing `load()` need to populate `SourceFile` metadata. Extract it into a reusable private method.

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h:91-102`
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp:1046-1056`

- [ ] **Step 1: Add private declaration to BrukerTimsFile.h**

In `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`, add after line 102 (`bool isDIA_(...)`):

```cpp
    /// Populate SourceFile metadata from the .d path (no peak data read)
    void loadExperimentalSettings_(const String& path, ExperimentalSettings& settings);
```

- [ ] **Step 2: Implement the helper in BrukerTimsFile.cpp**

In `src/openms/source/FORMAT/BrukerTimsFile.cpp`, add a new method (before the `load()` method, around line 1020):

```cpp
  void BrukerTimsFile::loadExperimentalSettings_(const String& path, ExperimentalSettings& settings)
  {
    SourceFile sf;
    sf.setNameOfFile(File::basename(path));
    sf.setPathToFile(File::path(path));
    sf.setFileType("Bruker TDF");
    sf.setNativeIDType("scan number only nativeID format");
    sf.setNativeIDTypeAccession("MS:1000776");
    settings.getSourceFiles().push_back(sf);
  }
```

- [ ] **Step 3: Replace inline SourceFile code in `load()` with helper call**

In `src/openms/source/FORMAT/BrukerTimsFile.cpp`, replace lines 1046-1056 (the SourceFile block after loading):

```cpp
    // Before:
    SourceFile sf;
    sf.setNameOfFile(File::basename(path));
    sf.setPathToFile(File::path(path));
    sf.setFileType("Bruker TDF");
    sf.setNativeIDType("scan number only nativeID format");
    sf.setNativeIDTypeAccession("MS:1000776");
    exp.getSourceFiles().push_back(sf);

    // After:
    loadExperimentalSettings_(path, exp);
```

- [ ] **Step 4: Build and verify existing tests still pass**

```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc)
```

Expected: Build succeeds. The refactoring is behavior-preserving.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h \
        src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "refactor: extract loadExperimentalSettings_() helper from BrukerTimsFile::load()"
```

---

### Task 2: Add `DIAStreamingMetadata` struct and `readDIAMetadata()` method

Read SWATH boundaries and spectrum counts from SQL metadata without touching peak data.

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h:78-82`
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Add includes and forward declarations to BrukerTimsFile.h**

In `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`, add after line 13 (`#include <string>`):

```cpp
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
```

- [ ] **Step 2: Add `DIAStreamingMetadata` struct and `readDIAMetadata()` declaration**

In `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`, add after line 77 (end of Config struct), before the `load()` methods:

```cpp
    /// Metadata for constructing a FullSwathFileConsumer for DIA streaming.
    struct DIAStreamingMetadata
    {
      std::vector<OpenSwath::SwathMap> boundaries;  ///< one SwathMap per DIAWindow (MS2 only)
      int nr_ms1_spectra = 0;                       ///< number of MS1 frames
      std::vector<int> nr_ms2_spectra;              ///< per-window spectrum counts (parallel to boundaries)
    };

    /// Read DIA SWATH boundaries and spectrum counts from a .d directory (SQL only, no peak data).
    /// Also populates exp_settings with source file metadata.
    DIAStreamingMetadata readDIAMetadata(const String& path, ExperimentalSettings& exp_settings,
                                         const Config& config = Config());
```

- [ ] **Step 3: Implement `readDIAMetadata()` in BrukerTimsFile.cpp**

Add the implementation. This opens a `TimsDataHandle` for IM calibration, reads DIAWindows and frame-to-WindowGroup mapping from SQL, computes boundaries and counts.

```cpp
  BrukerTimsFile::DIAStreamingMetadata BrukerTimsFile::readDIAMetadata(
      const String& path, ExperimentalSettings& exp_settings, const Config& config)
  {
    auto handle = openTimsDataHandle(path, config);
    String tdf_path = path + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    if (!isDIA(db))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "readDIAMetadata() requires a DIA dataset, but '" + path + "' appears to be DDA.");
    }

    loadExperimentalSettings_(path, exp_settings);

    // Read DIA windows (with IM conversion via handle's calibration)
    auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);

    DIAStreamingMetadata meta;
    if (windows.empty()) { return meta; }

    // Build SwathMap boundaries from DIAWindow structs
    meta.boundaries.reserve(windows.size());
    for (const auto& w : windows)
    {
      meta.boundaries.emplace_back(
        w.mz_center - w.mz_width / 2.0,  // lower
        w.mz_center + w.mz_width / 2.0,  // upper
        w.mz_center,                       // center
        w.im_lower,                        // imLower
        w.im_upper,                        // imUpper
        false);                            // ms1 = false
    }

    // Count MS1 frames
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ++meta.nr_ms1_spectra;
    }

    // Count MS2 spectra per window (= frames per WindowGroup)
    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
    meta.nr_ms2_spectra.reserve(windows.size());
    for (const auto& w : windows)
    {
      auto it = group_to_frames.find(w.window_group);
      meta.nr_ms2_spectra.push_back(
        it != group_to_frames.end() ? static_cast<int>(it->second.size()) : 0);
    }

    return meta;
  }
```

Note: `isDIA()` is the free function in the anonymous namespace (not `isDIA_()` which takes a path). Check which variant is used in context — the free function takes a `SQLite::Database&`. If only `isDIA_(tdf_path)` is available as a member, use that instead and open the DB after.

- [ ] **Step 4: Build to verify compilation**

```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc)
```

Expected: Build succeeds. No tests call the new method yet.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h \
        src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: add BrukerTimsFile::readDIAMetadata() for SQL-only boundary/count extraction"
```

---

### Task 3: Implement `loadDIAStreaming()` method

Stream MS1 + raw MS2 spectra to a `FullSwathFileConsumer` one-at-a-time.

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Add forward declaration and include for FullSwathFileConsumer**

In `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`, add a forward declaration after line 16 (`class TimsDataHandle;`):

```cpp
class FullSwathFileConsumer;
```

Then add the method declaration after `readDIAMetadata()`:

```cpp
    /// Stream DIA spectra to a consumer one-at-a-time without accumulating.
    /// MS2 spectra are always raw (no aggregation/denoising/centroiding).
    /// MS1 centroiding respects the Config settings.
    void loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer,
                          const Config& config = Config());
```

- [ ] **Step 2: Implement `loadDIAStreaming()` in BrukerTimsFile.cpp**

Add `#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>` to the includes at the top of BrukerTimsFile.cpp.

The implementation follows the same per-WindowGroup iteration as the raw path in `loadDIA_()`, but calls `consumer.consumeSpectrum(spec)` instead of `exp.addSpectrum()`:

```cpp
  void BrukerTimsFile::loadDIAStreaming(
      const String& path, FullSwathFileConsumer& consumer, const Config& config)
  {
    auto handle = openTimsDataHandle(path, config);
    String tdf_path = handle->get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    // --- MS1 frames ---
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;

    std::vector<uint32_t> ms1_frame_ids;
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ms1_frame_ids.push_back(fid);
    }

    startProgress(0, ms1_frame_ids.size(), "Streaming DIA-PASEF MS1 frames");
    for (size_t i = 0; i < ms1_frame_ids.size(); ++i)
    {
      TimsFrame& frame = handle->get_frame(ms1_frame_ids[i]);
      MSSpectrum spec;
      if (do_centroid)
        centroidMS1Frame(frame, spec, config, centroider);
      else
        frameToSpectrum(frame, spec, 1);
      consumer.consumeSpectrum(spec);
      setProgress(i);
    }
    endProgress();

    // --- MS2 frames: raw per-WindowGroup iteration (no aggregation) ---
    auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);
    if (windows.empty()) { return; }

    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);

    std::map<int, std::vector<const DIAWindow*>> group_to_windows;
    for (const auto& w : windows)
      group_to_windows[w.window_group].push_back(&w);

    Size total_work = 0;
    for (const auto& [group, frames] : group_to_frames)
    {
      auto wit = group_to_windows.find(group);
      if (wit != group_to_windows.end())
        total_work += frames.size() * wit->second.size();
    }

    startProgress(0, total_work, "Streaming DIA-PASEF MS2 frames");
    Size progress_count = 0;

    for (const auto& [group, frame_ids] : group_to_frames)
    {
      auto wit = group_to_windows.find(group);
      if (wit == group_to_windows.end()) continue;
      const auto& dia_windows = wit->second;

      for (const DIAWindow* win : dia_windows)
      {
        for (size_t i = 0; i < frame_ids.size(); ++i)
        {
          setProgress(progress_count++);
          TimsFrame& frame = handle->get_frame(frame_ids[i]);
          if (frame.num_peaks == 0) continue;

          std::vector<uint32_t> scan_ids(frame.num_peaks);
          std::vector<uint32_t> intensities(frame.num_peaks);
          std::vector<double> mzs(frame.num_peaks);
          std::vector<double> inv_ion_mobilities(frame.num_peaks);

          frame.save_to_buffs(nullptr, scan_ids.data(), nullptr, intensities.data(),
                              mzs.data(), inv_ion_mobilities.data(), nullptr);

          MSSpectrum spec;
          spec.setRT(frame.time);
          spec.setMSLevel(2);
          spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
          spec.setNativeID("frame=" + String(frame.id) + " windowGroup=" + String(win->window_group));

          Precursor prec;
          prec.setMZ(win->mz_center);
          prec.setIsolationWindowLowerOffset(win->mz_width / 2.0);
          prec.setIsolationWindowUpperOffset(win->mz_width / 2.0);
          spec.setPrecursors({prec});

          spec.setMetaValue("ion mobility lower limit", win->im_lower);
          spec.setMetaValue("ion mobility upper limit", win->im_upper);

          DataArrays::FloatDataArray im_array;
          IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

          for (uint32_t p = 0; p < frame.num_peaks; ++p)
          {
            if (inv_ion_mobilities[p] >= win->im_lower && inv_ion_mobilities[p] <= win->im_upper)
            {
              Peak1D peak;
              peak.setMZ(mzs[p]);
              peak.setIntensity(static_cast<double>(intensities[p]));
              spec.push_back(peak);
              im_array.push_back(static_cast<float>(inv_ion_mobilities[p]));
            }
          }

          if (!spec.empty())
          {
            spec.getFloatDataArrays().push_back(std::move(im_array));
            spec.setIMPeakType(IMPeakType::IM_PROFILE);
            consumer.consumeSpectrum(spec);
          }
        }
      }
    }
    endProgress();
  }
```

- [ ] **Step 3: Build to verify compilation**

```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc)
```

Expected: Build succeeds.

- [ ] **Step 4: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h \
        src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: add BrukerTimsFile::loadDIAStreaming() for consumer-based DIA loading"
```

---

### Task 4: Extend `SwathFile::loadBrukerTdf()` with readoptions support

Add the new overload that accepts `readoptions` and `tmp`, and make the old signature a wrapper.

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/SwathFile.h:84-94`
- Modify: `src/openms/source/FORMAT/SwathFile.cpp:288-326`

- [ ] **Step 1: Add new overload declaration to SwathFile.h**

In `src/openms/include/OpenMS/FORMAT/SwathFile.h`, replace lines 84-94:

```cpp
#ifdef WITH_OPENTIMS
    /**
      @brief Loads a Swath run from a Bruker .d (TDF) directory

      @param[in] file Path to a Bruker .d (TDF) directory
      @param[in] tmp Temporary directory (for cached data)
      @param[in,out] exp_meta Will be filled with ExperimentalSettings metadata
      @param[in] readoptions How spectra are accessed: "normal" (in-memory) or "cache" (disk-cached)
      @return Vector of SwathMap structures representing the loaded Swath maps
    */
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const String& file,
                                                    const String& tmp,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                    const String& readoptions);

    /// @brief Convenience overload: loads Bruker TDF in-memory (readoptions="normal")
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const String& file,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta);
#endif
```

- [ ] **Step 2: Add include for BrukerTimsFile.h in SwathFile.cpp**

In `src/openms/source/FORMAT/SwathFile.cpp`, add inside the `#ifdef WITH_OPENTIMS` block near the top includes:

```cpp
#include <OpenMS/FORMAT/BrukerTimsFile.h>
```

(Verify this include isn't already present.)

- [ ] **Step 3: Implement new `loadBrukerTdf()` overload**

Replace the existing `loadBrukerTdf()` implementation in SwathFile.cpp (lines 288-326) with both overloads:

```cpp
#ifdef WITH_OPENTIMS
  std::vector<OpenSwath::SwathMap> SwathFile::loadBrukerTdf(
    const String& file,
    const String& tmp,
    std::shared_ptr<ExperimentalSettings>& exp_meta,
    const String& readoptions)
  {
    OPENMS_LOG_INFO << "Loading Bruker TDF file " << file
                    << " using readoptions " << readoptions << '\n';
    startProgress(0, 1, "Loading Bruker TDF file " + file);

    BrukerTimsFile bruker_reader;
    bruker_reader.setLogType(this->getLogType());

    // Step 1: metadata from SQL (no peak data)
    ExperimentalSettings settings;
    auto meta = bruker_reader.readDIAMetadata(file, settings);
    auto exp_meta_ptr = std::make_shared<PeakMap>(settings);
    exp_meta = exp_meta_ptr;

    OPENMS_LOG_INFO << "Bruker TDF: " << meta.boundaries.size()
                    << " SWATH windows and " << meta.nr_ms1_spectra
                    << " MS1 spectra" << std::endl;

    // Step 2: construct consumer
    String tmp_fname = tmp.hasSuffix('/') ? File::getUniqueName() : "";
    std::shared_ptr<FullSwathFileConsumer> consumer;
    if (readoptions == "normal")
    {
      consumer = std::make_shared<RegularSwathFileConsumer>(meta.boundaries);
    }
    else if (readoptions == "cache")
    {
      consumer = std::make_shared<CachedSwathFileConsumer>(
        meta.boundaries, tmp, tmp_fname, meta.nr_ms1_spectra, meta.nr_ms2_spectra);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unknown or unsupported readoption '" + readoptions + "' for Bruker TDF loading");
    }
    consumer->setExperimentalSettings(settings);

    // Step 3: stream spectra to consumer
    bruker_reader.loadDIAStreaming(file, *consumer);

    // Step 4: finalize and retrieve SwathMaps
    std::vector<OpenSwath::SwathMap> swath_maps;
    consumer->retrieveSwathMaps(swath_maps);

    endProgress();
    return swath_maps;
  }

  std::vector<OpenSwath::SwathMap> SwathFile::loadBrukerTdf(
    const String& file,
    std::shared_ptr<ExperimentalSettings>& exp_meta)
  {
    return loadBrukerTdf(file, File::getTempDirectory(), exp_meta, "normal");
  }
#endif
```

- [ ] **Step 4: Build**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc)
```

Expected: Build succeeds.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/SwathFile.h \
        src/openms/source/FORMAT/SwathFile.cpp
git commit -m "feat: extend SwathFile::loadBrukerTdf() with readoptions/cache support"
```

---

### Task 5: Wire `readoptions` and `tmp` through OpenSwathBase

Update both call sites in `loadSwathFiles_()` to pass `readoptions` and `tmp` to the Bruker path.

**Files:**
- Modify: `src/openms/source/APPLICATIONS/OpenSwathBase.cpp:132-141,178-184`

- [ ] **Step 1: Update multi-file Bruker call site (line 135)**

In `src/openms/source/APPLICATIONS/OpenSwathBase.cpp`, replace the block at lines 132-141:

```cpp
#ifdef WITH_OPENTIMS
        else if (in_file_type == FileTypes::BRUKER_TDF)
        {
          auto maps = swath_file.loadBrukerTdf(f, tmp, exp_meta, readoptions);
          for (auto & m : maps)
          {
            swath_maps.push_back(m);
            swath_map_sources.push_back(f);
          }
        }
#endif
```

- [ ] **Step 2: Update single-file Bruker call site (line 181)**

Replace the block at lines 178-184:

```cpp
#ifdef WITH_OPENTIMS
      else if (in_file_type == FileTypes::BRUKER_TDF)
      {
        swath_maps = swath_file.loadBrukerTdf(file_list[0], tmp, exp_meta, readoptions);
        swath_map_sources.clear();
        for (Size i = 0; i < swath_maps.size(); ++i) swath_map_sources.push_back(file_list[0]);
      }
#endif
```

- [ ] **Step 3: Build**

```bash
cmake --build OpenMS-build -j$(nproc)
```

Expected: Build succeeds.

- [ ] **Step 4: Commit**

```bash
git add src/openms/source/APPLICATIONS/OpenSwathBase.cpp
git commit -m "feat: pass readoptions and tmp to Bruker TDF loading in OpenSwathBase"
```

---

### Task 6: Add tests for `readDIAMetadata()` and `loadDIAStreaming()`

Add test sections to the existing BrukerTimsFile test file.

**Files:**
- Modify: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`

- [ ] **Step 1: Add include for SwathFileConsumer**

At the top of `BrukerTimsFile_test.cpp`, add after the existing includes (around line 18):

```cpp
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
```

- [ ] **Step 2: Add `readDIAMetadata()` test section**

Insert before `#endif // OPENTIMS_DIA_TEST_DATA` (after the last DIA test section):

```cpp
START_SECTION(DIA readDIAMetadata test)
{
  BrukerTimsFile f;
  ExperimentalSettings settings;
  auto meta = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings);

  // Should have SWATH windows
  TEST_NOT_EQUAL(meta.boundaries.size(), 0);

  // Should have MS1 frames
  TEST_NOT_EQUAL(meta.nr_ms1_spectra, 0);

  // nr_ms2_spectra should have same size as boundaries
  TEST_EQUAL(meta.nr_ms2_spectra.size(), meta.boundaries.size());

  // Each window should have spectra
  for (Size i = 0; i < meta.nr_ms2_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(meta.nr_ms2_spectra[i], 0);
  }

  // Boundaries should have valid m/z and IM ranges
  for (const auto& b : meta.boundaries)
  {
    TEST_EQUAL(b.ms1, false);
    TEST_EQUAL(b.center > 0, true);
    TEST_EQUAL(b.lower < b.upper, true);
    TEST_EQUAL(b.imLower >= 0, true);
    TEST_EQUAL(b.imUpper > b.imLower, true);
  }

  // ExperimentalSettings should have source file
  TEST_EQUAL(settings.getSourceFiles().size(), 1);

  STATUS("readDIAMetadata: " << meta.boundaries.size() << " windows, "
         << meta.nr_ms1_spectra << " MS1 spectra");
}
END_SECTION
```

- [ ] **Step 3: Add `loadDIAStreaming()` test section**

Insert after the previous test section:

```cpp
START_SECTION(DIA loadDIAStreaming test)
{
  BrukerTimsFile f;

  // Get metadata first
  ExperimentalSettings settings;
  auto meta = f.readDIAMetadata(OPENTIMS_DIA_TEST_DATA, settings);

  // Create in-memory consumer with known boundaries
  RegularSwathFileConsumer consumer(meta.boundaries);
  consumer.setExperimentalSettings(settings);

  // Stream spectra
  f.loadDIAStreaming(OPENTIMS_DIA_TEST_DATA, consumer);

  // Retrieve SwathMaps
  std::vector<OpenSwath::SwathMap> swath_maps;
  consumer.retrieveSwathMaps(swath_maps);

  // Should have MS1 map + MS2 maps
  Size ms1_count = 0, ms2_count = 0;
  for (const auto& m : swath_maps)
  {
    if (m.ms1) ++ms1_count;
    else ++ms2_count;
  }
  TEST_EQUAL(ms1_count, 1);
  TEST_EQUAL(ms2_count, meta.boundaries.size());

  // MS1 map should have spectra
  for (const auto& m : swath_maps)
  {
    if (m.ms1)
    {
      TEST_EQUAL(m.sptr->getNrSpectra(), meta.nr_ms1_spectra);
      break;
    }
  }

  // Each MS2 map should have spectra with IM data
  for (const auto& m : swath_maps)
  {
    if (!m.ms1 && m.sptr->getNrSpectra() > 0)
    {
      auto spec = m.sptr->getSpectrumById(0);
      TEST_NOT_EQUAL(spec->getMZArray()->data.size(), 0);
      // Should have drift time array (per-peak IM)
      TEST_NOT_EQUAL(spec->getDriftTimeArray(), boost::shared_ptr<OpenSwath::BinaryDataArray>());
      break;
    }
  }

  STATUS("loadDIAStreaming: " << swath_maps.size() << " SwathMaps "
         << "(" << ms1_count << " MS1, " << ms2_count << " MS2)");
}
END_SECTION
```

- [ ] **Step 4: Build and run test**

```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc)
/path/to/OpenMS-build/src/tests/class_tests/bin/BrukerTimsFile_test
```

Expected: All tests PASS (including the new sections). The test will take ~8-10 minutes due to DIA data loading.

- [ ] **Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp
git commit -m "test: add tests for BrukerTimsFile::readDIAMetadata() and loadDIAStreaming()"
```
