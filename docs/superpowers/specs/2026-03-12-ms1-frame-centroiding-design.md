# MS1 Frame Centroiding for BrukerTimsFile

## Summary

Integrate Sage's `PeakBuffer` / `fastcentroid_frame()` algorithm into `BrukerTimsFile` as a config-driven load-time option. When enabled, MS1 frames are centroided across the ion mobility dimension, dramatically reducing peak count (e.g. 400k raw peaks to ~10k centroided peaks) while preserving per-peak IM values in the output `FloatDataArray`.

## Motivation

TimsTOF MS1 frames contain hundreds of thousands of raw peaks spread across the m/z and ion mobility dimensions. Downstream workflows (LFQ, feature finding) benefit from reduced peak counts, and the IM-dimension centroiding performed by Sage has proven effective across timsTOF Ultra and HT2 datasets. Currently OpenMS stores raw per-peak IM data but never centroids across the IM dimension.

## Reference Implementation and Attribution

The centroiding algorithm is adapted from **Sage** (Lazear, 2023), an open-source proteomics search engine:

- **Repository**: https://github.com/lazear/sage
- **Specific file**: [`crates/sage-cloudpath/src/tdf.rs`](https://github.com/lazear/sage/blob/064106c86fa2f89de04b1cd4aafe37d4a9ea3afb/crates/sage-cloudpath/src/tdf.rs) (commit `064106c`)
- **License**: MIT
- **Citation**: Lazear, M.R. (2023). Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale. *Journal of Proteome Research*, 22(11), 3652-3659. https://doi.org/10.1021/acs.jproteome.3c00486

Key components from Sage's `tdf.rs` that are translated to C++:

- `PeakBuffer` — reusable per-thread buffer holding `ImsPeak` structs (mz, intensity, im) → `FrameCentroider`
- `PeakBuffer::with_frame()` — converts raw frame data, sorts peaks by m/z, builds intensity-descending index → `FrameCentroider::loadFrame()`
- `PeakBuffer::fastcentroid_frame()` — DBSCAN-like centroiding algorithm → `FrameCentroider::centroid()`
- `PeakBuffer::expand_mobility_iter()` — expands run-length-encoded scan offsets to per-peak IM values → `expandScanOffsets()`

## Design

### 1. Config Changes (`BrukerTimsFile.h`)

Add two fields to `BrukerTimsFile::Config`:

```cpp
float ms1_centroid_mz_ppm = 0.0f;  ///< MS1 centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0)
float ms1_centroid_im_pct = 0.0f;  ///< MS1 centroiding IM tolerance in percent (0 = disabled, suggested: 3.0)
```

Centroiding is enabled when both values are > 0.

### 2. `FrameCentroider` Struct (private, in .cpp)

A C++ translation of Sage's `PeakBuffer`, defined in the anonymous namespace of `BrukerTimsFile.cpp`:

```cpp
struct ImsPeak
{
  float mz;
  float intensity;
  float im;
};

static constexpr size_t MAX_CENTROID_PEAKS = 10000;

struct FrameCentroider
{
  std::vector<ImsPeak> peaks;
  std::vector<size_t> order;       // indices sorted by descending intensity
  std::vector<ImsPeak> agg_buff;   // centroided output

  void clear();
  void loadFrame(const double* mz, const float* intensity, const float* im, uint32_t count);
  void centroid(float mz_ppm, float im_pct,
                std::vector<double>& out_mz, std::vector<double>& out_intensity,
                std::vector<float>& out_im);
};
```

**Algorithm** (translated from Sage's `fastcentroid_frame`):

1. `loadFrame()`: Convert input arrays (already-converted m/z doubles, not raw TOF indices) into `ImsPeak` structs, sort by m/z, build `order` vector sorted by descending intensity
2. `centroid()`: For each peak index in `order` (intensity-descending):
   - Skip if already consumed (`peaks[idx].intensity <= 0`)
   - If `agg_buff.size() > MAX_CENTROID_PEAKS`: log a debug message (`OPENMS_LOG_DEBUG`) if current intensity > 200, then break (note: `>` not `>=`, matching Sage — allows up to MAX+1 entries)
   - Compute m/z tolerance: `da_tol = mz * (mz_ppm / 1e6)`
   - Find m/z neighbor range: `std::lower_bound` for `mz - da_tol` (left bound), `std::upper_bound` for `mz + da_tol` (right bound)
   - Compute IM tolerance: `abs_im_tol = im * (im_pct / 100.0)`
   - Inner neighbor loop over the m/z range: for each neighbor, **skip if already consumed** (`intensity <= 0`), then check IM within `[im - abs_im_tol, im + abs_im_tol]`. If within tolerance: accumulate intensity, mark consumed (set intensity = -1), increment consumed counter
   - Emit centroided peak: apex m/z, summed intensity, apex IM
   - Break early if total consumed count equals total peak count
3. Sort `agg_buff` by m/z using `std::sort` (unstable sort is fine — centroided peaks have unique m/z from unique apexes), write to output vectors

All buffers persist across frames within a single `load()` call for memory reuse.

### 3. Scan-Offset Expansion Utility

Extract the duplicated inline scan-offset expansion (currently in both `frameToSpectrum_()` and `loadDIA_()`) into a free function:

```cpp
static void expandScanOffsets(const uint64_t* scan_offsets, uint32_t num_scans,
                              const double* scan_im, uint32_t num_peaks,
                              std::vector<float>& out_im);
```

This iterates `scan_offsets` windows, assigning the corresponding `scan_im[scan_idx]` value to each peak in the range `[scan_offsets[scan_idx], scan_offsets[scan_idx+1])`.

### 4. Integration into Frame Loading

#### Updated private method signatures

```cpp
void frameToSpectrum_(void* handle, const void* frame, MSSpectrum& spec,
                      const Config& config, FrameCentroider* centroider) const;
void loadDDA_(void* handle, MSExperiment& exp, const Config& config);
void loadDIA_(void* handle, MSExperiment& exp, const Config& config);
void loadFrames_(void* handle, MSExperiment& exp, const Config& config);
```

#### Flow in `frameToSpectrum_()`

1. Batch-convert TOF indices to m/z (existing)
2. Batch-convert scan indices to IM (existing)
3. Expand scan offsets to per-peak IM via `expandScanOffsets()` (refactored)
4. **If** `centroider != nullptr` **and** `frame->ms_level == 1`:
   - `centroider->loadFrame(mz, intensity, im, count)`
   - `centroider->centroid(config.ms1_centroid_mz_ppm, config.ms1_centroid_im_pct, out_mz, out_intensity, out_im)`
   - Use centroided output for MSSpectrum peaks and FloatDataArray
5. **Else**: use raw converted data (existing behavior)

#### Instantiation in load methods

Each of `loadDDA_()`, `loadDIA_()`, `loadFrames_()`:
- Check if centroiding is enabled (`config.ms1_centroid_mz_ppm > 0 && config.ms1_centroid_im_pct > 0`)
- If so, create one `FrameCentroider` instance and pass pointer to `frameToSpectrum_()`
- If not, pass `nullptr` (no centroiding)

#### Call site updates in `load()` and `transform()`

Both `load()` and `transform()` currently call `loadDDA_(ds.ptr, exp)` etc. without config. These call sites must be updated to forward the config:
- `load()`: `loadDDA_(ds.ptr, exp, config)`, `loadDIA_(ds.ptr, exp, config)`, `loadFrames_(ds.ptr, exp, config)`
- `transform()`: same — centroiding applies through the `transform()` path as well, since it uses the same `loadDDA_/loadDIA_/loadFrames_` methods internally

### 5. FileConverter TOPP Parameters

Add under the existing `timsrust:` subsection in `FileConverter.cpp`:

```cpp
registerDoubleOption_("timsrust:ms1_centroid_mz_ppm", "<float>", 0.0,
  "MS1 frame centroiding m/z tolerance in ppm. Collapses the ion mobility dimension "
  "by aggregating peaks within this m/z and IM tolerance. 0 = disabled, suggested value: 5.0",
  false, true);
setMinFloat_("timsrust:ms1_centroid_mz_ppm", 0.0);

registerDoubleOption_("timsrust:ms1_centroid_im_pct", "<float>", 0.0,
  "MS1 frame centroiding IM tolerance in percent. Used together with ms1_centroid_mz_ppm. "
  "0 = disabled, suggested value: 3.0",
  false, true);
setMinFloat_("timsrust:ms1_centroid_im_pct", 0.0);
```

Wire into `getTimsConfig_()`:

```cpp
c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("timsrust:ms1_centroid_mz_ppm"));
c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("timsrust:ms1_centroid_im_pct"));
```

### 6. `loadDIA_()` Refactoring

The scan-offset expansion in `loadDIA_()` (for MS2 SWATH splitting) is replaced with a call to `expandScanOffsets()`. The existing `std::vector<double> im_values` is kept as-is in `loadDIA_()` (the SWATH window IM bounds from `tims_swath_window` are likely double, so double comparisons are appropriate). The `expandScanOffsets()` utility outputs `std::vector<float>` for `frameToSpectrum_()` use; `loadDIA_()` will use a separate overload or keep its own inline expansion with double output.

### 7. Relationship Between TOF and IM Centroiding

In AUTO/SPECTRUM modes, timsrust already applies TOF-domain smoothing and centroiding during spectrum assembly. The IM centroiding described here is a separate, complementary step that collapses the ion mobility dimension within each frame. In FRAME mode, no TOF processing is applied by timsrust, so the IM centroiding operates on fully raw data. Both scenarios are valid: the IM centroiding reduces per-frame peak count regardless of whether TOF centroiding was applied first.

## Files Changed

| File | Change |
|---|---|
| `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` | Config fields, updated private method signatures |
| `src/openms/source/FORMAT/BrukerTimsFile.cpp` | `FrameCentroider`, `expandScanOffsets`, integration into load path |
| `src/topp/FileConverter.cpp` | TOPP parameters + wiring |

## Defaults and Behavior

- **Disabled by default**: both ppm and pct are 0.0, preserving current raw output
- **Suggested values**: 5.0 ppm m/z, 3.0% IM (matching Sage's proven defaults)
- **Applies to**: MS1 frames only, in all export modes (AUTO, SPECTRUM, FRAME)
- **Output format**: same as current — MSSpectrum with FloatDataArray for per-peak IM, just with fewer peaks
- **MAX_CENTROID_PEAKS = 10000**: matches Sage; in practice rarely reached for single frames

## Test Plan

1. **Unit tests for `FrameCentroider`** (in `BrukerTimsFile_test.cpp`):
   - Synthetic peak data with known expected centroid output (e.g., two clusters that should merge)
   - Edge cases: empty frame, single peak, all peaks within tolerance (single centroid output)
   - MAX_CENTROID_PEAKS limit: verify output is capped and remaining peaks are not lost silently
   - Verify that already-consumed peaks are not double-counted in later iterations

2. **Integration tests** (require `TIMSRUST_DDA_TEST_DATA` / `TIMSRUST_DIA_TEST_DATA`):
   - Load with centroiding disabled: verify output matches current behavior
   - Load with centroiding enabled (5.0 ppm, 3.0%): verify MS1 peak count is reduced
   - Verify MS2 spectra are unaffected by centroiding config
   - Verify centroided MS1 spectra still have valid FloatDataArray with per-peak IM

3. **FileConverter smoke test**:
   - Verify new TOPP parameters are accepted and passed through to Config
