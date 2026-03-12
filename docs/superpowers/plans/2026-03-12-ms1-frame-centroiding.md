# MS1 Frame Centroiding Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add config-driven IM-dimension centroiding for MS1 frames in BrukerTimsFile, adapted from Sage's fastcentroid_frame algorithm.

**Architecture:** A `FrameCentroider` struct (hidden in the .cpp anonymous namespace) performs DBSCAN-like centroiding on converted frame peaks. Load methods (`loadDDA_`, `loadDIA_`, `loadFrames_`) instantiate the centroider when config enables it and apply it to MS1 frames before building `MSSpectrum` objects. FileConverter exposes the two tolerance parameters as TOPP options.

**Tech Stack:** C++17, timsrust_cpp_bridge FFI, OpenMS MSExperiment/MSSpectrum API

**Spec:** `docs/superpowers/specs/2026-03-12-ms1-frame-centroiding-design.md`

---

## Chunk 1: Foundation — Config, expandScanOffsets, FrameCentroider, refactored call sites

All changes in this chunk are committed together to keep the build green at every commit.

### Task 1: Add config fields, new helpers, and refactor existing code

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h:53-82`
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Add config fields to header**

In `BrukerTimsFile.h`, add two fields to `struct Config` after the `calibrate` field (line 58), before the `ExportMode` enum:

```cpp
      float ms1_centroid_mz_ppm = 0.0f; ///< MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0). Adapted from Sage (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
      float ms1_centroid_im_pct = 0.0f;  ///< MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)
```

- [ ] **Step 2: Update private method signatures in header**

In the same header, update the private method declarations (lines 76-82) to accept `const Config&`:

```cpp
    void loadDDA_(void* handle, MSExperiment& exp, const Config& config);
    void loadDIA_(void* handle, MSExperiment& exp, const Config& config);
    void loadFrames_(void* handle, MSExperiment& exp, const Config& config);
```

`frameToSpectrum_()` and `isDIA_()` signatures remain unchanged.

- [ ] **Step 3: Add expandScanOffsets template function to .cpp**

Insert after line 56 (after `TimsConfigHandle` struct, before `getTimsError`), inside `namespace OpenMS`:

```cpp
  // Helper: expand run-length-encoded scan offsets to per-peak IM values.
  // Translated from Sage's PeakBuffer::expand_mobility_iter().
  // T = float for frameToSpectrum_() / centroiding path,
  // T = double for loadDIA_() SWATH splitting (preserves precision for window boundary comparisons).
  template <typename T>
  static void expandScanOffsets(const uint64_t* scan_offsets, uint32_t num_scans,
                                const double* scan_im, uint32_t num_peaks,
                                std::vector<T>& out_im)
  {
    out_im.resize(num_peaks);
    for (uint32_t scan_idx = 0; scan_idx < num_scans; ++scan_idx)
    {
      uint64_t start = scan_offsets[scan_idx];
      uint64_t end = scan_offsets[scan_idx + 1];
      for (uint64_t p = start; p < end; ++p)
      {
        out_im[p] = static_cast<T>(scan_im[scan_idx]);
      }
    }
  }
```

- [ ] **Step 4: Add ImsPeak struct and MAX_CENTROID_PEAKS constant**

Insert right after `expandScanOffsets`:

```cpp
  // IM-aware peak for centroiding. Adapted from Sage's ImsPeak.
  struct ImsPeak
  {
    float mz;
    float intensity;
    float im;
  };

  static constexpr size_t MAX_CENTROID_PEAKS = 10000;
```

- [ ] **Step 5: Add FrameCentroider struct**

Insert right after `MAX_CENTROID_PEAKS`. This is a direct C++ translation of Sage's `PeakBuffer`:

```cpp
  // Reusable buffer for IM-dimension centroiding of a single frame.
  // Translated from Sage's PeakBuffer (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
  // Buffers persist across frames within a single load() call for memory reuse.
  struct FrameCentroider
  {
    std::vector<ImsPeak> peaks;
    std::vector<size_t> order;       // indices sorted by descending intensity
    std::vector<ImsPeak> agg_buff;   // centroided output

    void clear()
    {
      peaks.clear();
      order.clear();
      agg_buff.clear();
    }

    // Load already-converted frame data into the buffer.
    // mz_values: converted m/z (double, from tims_convert_tof_to_mz_array)
    // intensities: raw uint32_t from tims_frame::intensities (cast to float internally)
    // im_values: expanded per-peak IM (float, from expandScanOffsets)
    void loadFrame(const double* mz_values, const uint32_t* intensities,
                   const float* im_values, uint32_t count)
    {
      clear();
      peaks.reserve(count);
      for (uint32_t i = 0; i < count; ++i)
      {
        peaks.push_back({static_cast<float>(mz_values[i]),
                         static_cast<float>(intensities[i]),
                         im_values[i]});
      }

      // Sort by m/z for binary-search neighbor finding
      std::sort(peaks.begin(), peaks.end(),
                [](const ImsPeak& a, const ImsPeak& b) { return a.mz < b.mz; });

      // Build intensity-descending index
      order.resize(count);
      std::iota(order.begin(), order.end(), size_t(0));
      std::sort(order.begin(), order.end(),
                [this](size_t a, size_t b)
                {
                  return peaks[b].intensity < peaks[a].intensity; // descending
                });
    }

    // Centroid the loaded frame by collapsing the IM dimension.
    // Iterates peaks in descending intensity order. For each apex peak, aggregates
    // all unconsumed neighbors within m/z (ppm) and IM (percent) tolerances.
    // Output is capped at MAX_CENTROID_PEAKS entries.
    void centroid(float mz_ppm, float im_pct,
                  std::vector<double>& out_mz,
                  std::vector<double>& out_intensity,
                  std::vector<float>& out_im)
    {
      agg_buff.clear();
      agg_buff.reserve(std::min(peaks.size(), MAX_CENTROID_PEAKS + 1));

      const float utol = mz_ppm / 1e6f;
      const float im_tol_frac = im_pct / 100.0f;
      size_t global_consumed = 0;

      for (size_t idx : order)
      {
        if (peaks[idx].intensity <= 0.0f) continue;  // already consumed

        if (agg_buff.size() > MAX_CENTROID_PEAKS)
        {
          if (peaks[idx].intensity > 200.0f)
          {
            OPENMS_LOG_DEBUG << "FrameCentroider: reached MAX_CENTROID_PEAKS at index "
                             << idx << "/" << peaks.size()
                             << " intensity=" << peaks[idx].intensity << std::endl;
          }
          break;
        }

        const float mz = peaks[idx].mz;
        const float im = peaks[idx].im;
        const float da_tol = mz * utol;
        const float left_mz = mz - da_tol;
        const float right_mz = mz + da_tol;

        // Binary search for m/z neighbor range
        auto it_start = std::lower_bound(peaks.begin(), peaks.end(), left_mz,
                          [](const ImsPeak& p, float val) { return p.mz < val; });
        auto it_end = std::upper_bound(peaks.begin(), peaks.end(), right_mz,
                          [](float val, const ImsPeak& p) { return val < p.mz; });

        const float abs_im_tol = im * im_tol_frac;
        const float left_im = im - abs_im_tol;
        const float right_im = im + abs_im_tol;

        float curr_intensity = 0.0f;
        size_t num_consumed = 0;

        for (auto it = it_start; it != it_end; ++it)
        {
          if (it->intensity <= 0.0f) continue;  // already consumed by earlier apex
          if (it->im >= left_im && it->im <= right_im)
          {
            curr_intensity += it->intensity;
            it->intensity = -1.0f;  // mark consumed
            ++num_consumed;
          }
        }

        agg_buff.push_back({mz, curr_intensity, im});
        global_consumed += num_consumed;

        if (global_consumed == peaks.size()) break;  // all peaks assigned
      }

      // Sort output by m/z
      std::sort(agg_buff.begin(), agg_buff.end(),
                [](const ImsPeak& a, const ImsPeak& b) { return a.mz < b.mz; });

      // Write to output vectors
      out_mz.resize(agg_buff.size());
      out_intensity.resize(agg_buff.size());
      out_im.resize(agg_buff.size());
      for (size_t i = 0; i < agg_buff.size(); ++i)
      {
        out_mz[i] = static_cast<double>(agg_buff[i].mz);
        out_intensity[i] = static_cast<double>(agg_buff[i].intensity);
        out_im[i] = agg_buff[i].im;
      }
    }
  };
```

- [ ] **Step 6: Update load() call sites**

In `BrukerTimsFile::load()` (around line 131-142), change each call to forward config:
```cpp
      loadFrames_(ds.ptr, exp);       →  loadFrames_(ds.ptr, exp, config);
      loadDDA_(ds.ptr, exp);          →  loadDDA_(ds.ptr, exp, config);
      loadDIA_(ds.ptr, exp);          →  loadDIA_(ds.ptr, exp, config);
```

- [ ] **Step 7: Update transform() call sites**

In `BrukerTimsFile::transform()` (around lines 204-209), apply the same three changes:
```cpp
      loadFrames_(ds.ptr, exp);       →  loadFrames_(ds.ptr, exp, config);
      loadDIA_(ds.ptr, exp);          →  loadDIA_(ds.ptr, exp, config);
      loadDDA_(ds.ptr, exp);          →  loadDDA_(ds.ptr, exp, config);
```

- [ ] **Step 8: Update loadDDA_ signature**

Change line 218:
```cpp
  void BrukerTimsFile::loadDDA_(void* handle, MSExperiment& exp)
```
to:
```cpp
  void BrukerTimsFile::loadDDA_(void* handle, MSExperiment& exp, const Config& config)
```

(The `config` parameter is unused in this task — centroiding integration comes in Task 2.)

- [ ] **Step 9: Update loadDIA_ signature and refactor scan-offset expansion**

Change line 300:
```cpp
  void BrukerTimsFile::loadDIA_(void* handle, MSExperiment& exp)
```
to:
```cpp
  void BrukerTimsFile::loadDIA_(void* handle, MSExperiment& exp, const Config& config)
```

Replace the inline scan-offset expansion (lines 367-377):
```cpp
      // Build per-peak IM from scan offsets
      std::vector<double> im_values(frame.num_peaks);
      for (uint32_t scan_idx = 0; scan_idx < frame.num_scans; ++scan_idx)
      {
        uint64_t start = frame.scan_offsets[scan_idx];
        uint64_t end = frame.scan_offsets[scan_idx + 1];
        for (uint64_t p = start; p < end; ++p)
        {
          im_values[p] = scan_im[scan_idx];
        }
      }
```
with:
```cpp
      // Build per-peak IM from scan offsets
      std::vector<double> im_values;
      expandScanOffsets(frame.scan_offsets, frame.num_scans, scan_im.data(), frame.num_peaks, im_values);
```

- [ ] **Step 10: Update loadFrames_ signature**

Change line 422:
```cpp
  void BrukerTimsFile::loadFrames_(void* handle, MSExperiment& exp)
```
to:
```cpp
  void BrukerTimsFile::loadFrames_(void* handle, MSExperiment& exp, const Config& config)
```

- [ ] **Step 11: Refactor frameToSpectrum_ to use expandScanOffsets**

In `frameToSpectrum_()` (around lines 482-492), replace the inline scan-offset expansion:
```cpp
    // Build per-peak IM values from scan offsets
    std::vector<float> im_values(frame->num_peaks);
    for (uint32_t scan_idx = 0; scan_idx < frame->num_scans; ++scan_idx)
    {
      uint64_t start = frame->scan_offsets[scan_idx];
      uint64_t end = frame->scan_offsets[scan_idx + 1];
      for (uint64_t p = start; p < end; ++p)
      {
        im_values[p] = static_cast<float>(scan_im[scan_idx]);
      }
    }
```
with:
```cpp
    // Build per-peak IM values from scan offsets
    std::vector<float> im_values;
    expandScanOffsets(frame->scan_offsets, frame->num_scans, scan_im.data(), frame->num_peaks, im_values);
```

- [ ] **Step 12: Build and verify**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -20`
Expected: Clean compilation, no errors. The new helpers are defined but not yet called from centroiding logic.

- [ ] **Step 13: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat(BrukerTimsFile): add FrameCentroider, expandScanOffsets, thread Config

Adds FrameCentroider (adapted from Sage's PeakBuffer, Lazear 2023,
doi:10.1021/acs.jproteome.3c00486) and expandScanOffsets<T>() utility.
Threads Config through loadDDA_/loadDIA_/loadFrames_ from load() and
transform(). Replaces duplicated inline scan-offset expansion."
```

---

## Chunk 2: Centroiding Integration and FileConverter

### Task 2: Integrate FrameCentroider into load methods for MS1 frames

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Add partial-config validation helper**

Insert after the `FrameCentroider` struct definition (before `getTimsError`):

```cpp
  // Helper: check if MS1 centroiding is enabled, warn on partial config
  static bool isCentroidingEnabled(const BrukerTimsFile::Config& config)
  {
    bool has_mz = config.ms1_centroid_mz_ppm > 0.0f;
    bool has_im = config.ms1_centroid_im_pct > 0.0f;
    if (has_mz != has_im)
    {
      OPENMS_LOG_WARN << "Warning: ms1_centroid_mz_ppm and ms1_centroid_im_pct must both be > 0 "
                      << "to enable MS1 centroiding. Centroiding disabled." << std::endl;
      return false;
    }
    return has_mz && has_im;
  }
```

- [ ] **Step 2: Add centroidMS1Frame_ helper free function**

Insert after `isCentroidingEnabled`. This function handles the full centroiding pipeline for a single MS1 frame, avoiding code duplication across the three load methods:

```cpp
  // Centroid a single MS1 frame: convert TOF->mz, expand scan offsets, run centroider,
  // build MSSpectrum with centroided peaks and IM FloatDataArray.
  // Note: This duplicates the TOF-to-mz and scan-to-IM conversion from frameToSpectrum_()
  // because centroiding needs the intermediate arrays before they are packed into MSSpectrum.
  static void centroidMS1Frame_(tims_dataset* ds, const tims_frame& frame,
                                MSSpectrum& spec, const BrukerTimsFile::Config& config,
                                FrameCentroider& centroider)
  {
    spec.clear(true);
    spec.setRT(frame.rt_seconds);
    spec.setMSLevel(frame.ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    spec.setType(SpectrumSettings::CENTROID);

    if (frame.num_peaks == 0) return;

    // Batch-convert TOF indices to m/z
    std::vector<double> mz_values(frame.num_peaks);
    timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame.tof_indices, frame.num_peaks, mz_values.data());
    if (mz_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "TOF to m/z conversion failed: " + getTimsError(ds));

    // Batch-convert scan indices to IM
    std::vector<uint32_t> scan_indices(frame.num_scans);
    std::iota(scan_indices.begin(), scan_indices.end(), 0u);
    std::vector<double> scan_im(frame.num_scans);
    timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame.num_scans, scan_im.data());
    if (im_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Scan to IM conversion failed: " + getTimsError(ds));

    // Expand scan offsets to per-peak IM
    std::vector<float> im_values;
    expandScanOffsets(frame.scan_offsets, frame.num_scans, scan_im.data(), frame.num_peaks, im_values);

    // Run centroiding
    centroider.loadFrame(mz_values.data(), frame.intensities, im_values.data(), frame.num_peaks);

    std::vector<double> cent_mz, cent_intensity;
    std::vector<float> cent_im;
    centroider.centroid(config.ms1_centroid_mz_ppm, config.ms1_centroid_im_pct,
                        cent_mz, cent_intensity, cent_im);

    // Build MSSpectrum from centroided output
    spec.reserve(cent_mz.size());
    for (size_t i = 0; i < cent_mz.size(); ++i)
    {
      Peak1D peak;
      peak.setMZ(cent_mz[i]);
      peak.setIntensity(cent_intensity[i]);
      spec.push_back(peak);
    }

    // Set IM float data array
    DataArrays::FloatDataArray im_array;
    im_array.assign(cent_im.begin(), cent_im.end());
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    spec.getFloatDataArrays().push_back(std::move(im_array));
  }
```

- [ ] **Step 3: Integrate centroiding into loadDDA_()**

In `loadDDA_()`, after the RAII guard and before the MS1 frame loop (around line 237), add centroider setup:

Replace the MS1 frame loop (lines 239-245):
```cpp
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
```
with:
```cpp
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame_(ds, ms1_frames[i], spec, config, centroider);
      }
      else
      {
        frameToSpectrum_(handle, &ms1_frames[i], spec);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
```

- [ ] **Step 4: Integrate centroiding into loadDIA_()**

In `loadDIA_()`, replace the MS1 frame loop (lines 314-322):
```cpp
    startProgress(0, ms1_count, "Loading DIA-PASEF MS1 frames");
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
    endProgress();
```
with:
```cpp
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    startProgress(0, ms1_count, "Loading DIA-PASEF MS1 frames");
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame_(ds, ms1_frames[i], spec, config, centroider);
      }
      else
      {
        frameToSpectrum_(handle, &ms1_frames[i], spec);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
    endProgress();
```

- [ ] **Step 5: Integrate centroiding into loadFrames_()**

In `loadFrames_()`, replace the inner frame loop (lines 442-450):
```cpp
      startProgress(0, count, String("Loading MS") + String(level) + " frames");
      for (unsigned int i = 0; i < count; ++i)
      {
        MSSpectrum spec;
        frameToSpectrum_(handle, &frames[i], spec);
        exp.addSpectrum(std::move(spec));
        setProgress(i);
      }
      endProgress();
```
with:
```cpp
      bool do_centroid = (level == 1) && isCentroidingEnabled(config);
      FrameCentroider centroider;
      startProgress(0, count, String("Loading MS") + String(level) + " frames");
      for (unsigned int i = 0; i < count; ++i)
      {
        MSSpectrum spec;
        if (do_centroid)
        {
          centroidMS1Frame_(ds, frames[i], spec, config, centroider);
        }
        else
        {
          frameToSpectrum_(handle, &frames[i], spec);
        }
        exp.addSpectrum(std::move(spec));
        setProgress(i);
      }
      endProgress();
```

- [ ] **Step 6: Build and verify**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -20`
Expected: Clean compilation, no errors.

- [ ] **Step 7: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat(BrukerTimsFile): integrate FrameCentroider into MS1 frame loading

When ms1_centroid_mz_ppm and ms1_centroid_im_pct are both > 0, MS1 frames
are centroided across the IM dimension before building MSSpectrum objects.
Centroided spectra are marked with SpectrumSettings::CENTROID type.
Algorithm adapted from Sage (Lazear 2023, doi:10.1021/acs.jproteome.3c00486)."
```

---

### Task 3: Add FileConverter TOPP parameters

**Files:**
- Modify: `src/topp/FileConverter.cpp:182-197` (parameter registration)
- Modify: `src/topp/FileConverter.cpp:219-231` (getTimsConfig_)

- [ ] **Step 1: Register TOPP parameters**

In `registerOptionsAndFlags_()`, after the `timsrust:export_mode` registration (line 196), add before the `#endif`:

```cpp
    registerDoubleOption_("timsrust:ms1_centroid_mz_ppm", "<float>", 0.0,
      "MS1 frame IM-centroiding m/z tolerance in ppm. Collapses the ion mobility dimension "
      "by aggregating neighboring peaks. Both this and ms1_centroid_im_pct must be > 0 to enable. "
      "Suggested value: 5.0. Algorithm from Sage (Lazear 2023).", false, true);
    setMinFloat_("timsrust:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("timsrust:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
    setMinFloat_("timsrust:ms1_centroid_im_pct", 0.0);
```

- [ ] **Step 2: Wire into getTimsConfig_()**

In `getTimsConfig_()`, after the export_mode handling (line 229), add before `return c;`:

```cpp
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("timsrust:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("timsrust:ms1_centroid_im_pct"));
```

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build --target FileConverter -j$(nproc) 2>&1 | tail -20`
Expected: Clean compilation.

- [ ] **Step 4: Verify parameter appears in help**

Run: `OpenMS-build/bin/FileConverter --helphelp 2>&1 | grep ms1_centroid`
Expected: Both `ms1_centroid_mz_ppm` and `ms1_centroid_im_pct` appear in output.

- [ ] **Step 5: Commit**

```bash
git add src/topp/FileConverter.cpp
git commit -m "feat(FileConverter): expose MS1 centroiding params as TOPP options

Adds timsrust:ms1_centroid_mz_ppm and timsrust:ms1_centroid_im_pct
parameters for controlling IM-dimension centroiding of MS1 frames."
```

---

## Chunk 3: Tests

### Task 4: Add integration tests for MS1 centroiding

**Files:**
- Modify: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`

- [ ] **Step 1: Add DDA centroiding test**

After the existing `DDA frame-level loading test` section (before `#endif // TIMSRUST_DDA_TEST_DATA` around line 153), add:

```cpp
START_SECTION(DDA MS1 centroiding test)
{
  BrukerTimsFile f;

  // Load without centroiding (baseline)
  MSExperiment exp_raw;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_raw);

  // Load with centroiding enabled
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_cent, cfg);

  // Count MS1 and MS2 spectra and total MS1 peaks in both
  Size raw_ms1_peaks = 0, cent_ms1_peaks = 0;
  Size raw_ms1_count = 0, cent_ms1_count = 0;
  Size raw_ms2_count = 0, cent_ms2_count = 0;

  for (const auto& spec : exp_raw)
  {
    if (spec.getMSLevel() == 1) { ++raw_ms1_count; raw_ms1_peaks += spec.size(); }
    else { ++raw_ms2_count; }
  }
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 1) { ++cent_ms1_count; cent_ms1_peaks += spec.size(); }
    else { ++cent_ms2_count; }
  }

  // MS1 frame count should be the same
  TEST_EQUAL(raw_ms1_count, cent_ms1_count);

  // MS2 count should be unaffected by centroiding
  TEST_EQUAL(raw_ms2_count, cent_ms2_count);

  // Centroided MS1 should have fewer total peaks
  TEST_EQUAL(cent_ms1_peaks < raw_ms1_peaks, true);

  // Centroided MS1 spectra should still have IM data, be marked CENTROID,
  // have IM array length == peak count, and be sorted by m/z
  for (const auto& spec : exp_cent)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_EQUAL(spec.getType(), SpectrumSettings::CENTROID);

      // IM array length must match peak count
      const auto& fda = spec.getFloatDataArrays();
      TEST_EQUAL(fda.empty(), false);
      TEST_EQUAL(fda[0].size(), spec.size());

      // Peaks must be sorted by m/z (centroider sorts output)
      bool mz_sorted = true;
      for (Size j = 1; j < spec.size(); ++j)
      {
        if (spec[j].getMZ() < spec[j-1].getMZ()) { mz_sorted = false; break; }
      }
      TEST_EQUAL(mz_sorted, true);
      break;
    }
  }
}
END_SECTION

START_SECTION(DDA partial centroiding config test)
{
  // Setting only one tolerance should NOT enable centroiding
  BrukerTimsFile f;

  BrukerTimsFile::Config cfg_partial;
  cfg_partial.ms1_centroid_mz_ppm = 5.0f;  // set
  cfg_partial.ms1_centroid_im_pct = 0.0f;  // NOT set
  MSExperiment exp_partial;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_partial, cfg_partial);

  // Load without centroiding for comparison
  MSExperiment exp_raw;
  f.load(TIMSRUST_DDA_TEST_DATA, exp_raw);

  // MS1 peak counts should be identical (centroiding was NOT applied)
  Size partial_ms1_peaks = 0, raw_ms1_peaks = 0;
  for (const auto& spec : exp_partial)
    if (spec.getMSLevel() == 1) partial_ms1_peaks += spec.size();
  for (const auto& spec : exp_raw)
    if (spec.getMSLevel() == 1) raw_ms1_peaks += spec.size();

  TEST_EQUAL(partial_ms1_peaks, raw_ms1_peaks);
}
END_SECTION
```

- [ ] **Step 2: Add DIA centroiding test (if DIA test data available)**

After the existing `DIA round-trip test` section (before `#endif // TIMSRUST_DIA_TEST_DATA` around line 209), add:

```cpp
START_SECTION(DIA MS1 centroiding test)
{
  BrukerTimsFile f;

  // Load without centroiding
  MSExperiment exp_raw;
  f.load(TIMSRUST_DIA_TEST_DATA, exp_raw);

  // Load with centroiding
  BrukerTimsFile::Config cfg;
  cfg.ms1_centroid_mz_ppm = 5.0f;
  cfg.ms1_centroid_im_pct = 3.0f;
  MSExperiment exp_cent;
  f.load(TIMSRUST_DIA_TEST_DATA, exp_cent, cfg);

  // Count MS1 peaks
  Size raw_ms1_peaks = 0, cent_ms1_peaks = 0;
  for (const auto& spec : exp_raw)
    if (spec.getMSLevel() == 1) raw_ms1_peaks += spec.size();
  for (const auto& spec : exp_cent)
    if (spec.getMSLevel() == 1) cent_ms1_peaks += spec.size();

  // Centroided should have fewer peaks
  TEST_EQUAL(cent_ms1_peaks < raw_ms1_peaks, true);
}
END_SECTION
```

- [ ] **Step 3: Build test**

Run: `cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc) 2>&1 | tail -20`
Expected: Clean compilation.

- [ ] **Step 4: Run test (if test data available)**

Run: `ctest --test-dir OpenMS-build -R BrukerTimsFile_test -V 2>&1 | tail -40`
Expected: All tests pass (or skip gracefully if test data macros are not defined).

- [ ] **Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp
git commit -m "test(BrukerTimsFile): add integration tests for MS1 IM-centroiding

Verifies centroiding reduces MS1 peak count, leaves MS2 unaffected,
preserves IM FloatDataArray, and marks spectra as CENTROID type."
```
