# TIMS Calibration: Tiered Scan-to-1/K0 Conversion for Bruker TimsTOF Data

**Date:** 2026-03-25
**Status:** Draft
**Scope:** `BrukerTimsFile` reader — scan index → 1/K0 conversion accuracy

## Problem

OpenMS currently uses a simple linear model to convert TIMS scan indices to inverse
reduced ion mobility (1/K0) values:

```
1/K0 = OneOverK0AcqRangeUpper + (OneOverK0AcqRangeLower - OneOverK0AcqRangeUpper) / MAX(NumScans) × scan_index
```

This is a frame-independent approximation read from the `GlobalMetadata` table in
`analysis.tdf`. Every other major open-source tool (OpenTIMS, AlphaTIMS, timsrust)
uses the same linear model. The Bruker proprietary SDK provides higher accuracy via
per-frame calibration and pressure compensation, but no open-source tool exploits the
calibration coefficients stored in the `TimsCalibration` table.

## Solution

Implement a **tiered calibration strategy** in `BrukerTimsFile`:

1. **Bruker SDK** — highest accuracy, per-frame, with optional pressure compensation.
   Requires `timsdata.dll`/`libtimsdata.so` at runtime.
2. **Rational function model** — first open-source per-frame calibration using the
   `TimsCalibration` table coefficients. No SDK required.
3. **Linear model** — current behavior, used as ultimate fallback.

The strategy is configured via `BrukerTimsFile::Config` with an `AUTO` default that
tries SDK → rational → linear.

## Background: Landscape Survey

| Tool | Bruker SDK? | Open-source model | Per-frame? | Pressure comp? |
|---|---|---|---|---|
| OpenMS (current) | No | Linear (GlobalMetadata) | No | No |
| OpenTIMS | Optional | Linear (GlobalMetadata) | Only via SDK | No |
| AlphaTIMS | Preferred (bundled) | Linear fallback | No (uses frame 1 globally) | No |
| timsrust | No | Linear (GlobalMetadata) | No | No |
| ProteoWizard | Required | None | Via SDK | No (legacy `tims_open`) |
| diapysef | Optional | Rational (TimsCalibration) | No (single row) | No |
| TIMSCONVERT/pyTDFSDK | Yes | None | Via SDK | Yes (only tool) |

**Key finding:** No open-source tool reads per-frame calibration from the
`TimsCalibration` table. diapysef reads the table but ignores `ModelType` and does
not map frames to calibration segments. Only TIMSCONVERT exposes pressure
compensation.

## Design

### Config Additions

```cpp
struct Config
{
  // ... existing fields ...

  /// Strategy for converting TIMS scan indices to 1/K0 values.
  /// AUTO (default): tries Bruker SDK → rational → linear.
  enum class TimsCalibrationStrategy { AUTO, BRUKER_SDK, RATIONAL, LINEAR };
  TimsCalibrationStrategy tims_calibration_strategy = TimsCalibrationStrategy::AUTO;

  /// Pressure compensation strategy (only effective with BRUKER_SDK calibration).
  /// Corrects for ambient gas pressure drift during acquisition.
  /// Ignored when using RATIONAL or LINEAR strategies.
  enum class PressureCompensation { NONE, GLOBAL, PER_FRAME };
  PressureCompensation pressure_compensation = PressureCompensation::NONE;

  /// Path to Bruker SDK library (timsdata.dll / libtimsdata.so).
  /// Empty string (default) = discover from OPENMS_BRUKER_SDK_PATH env var.
  std::string bruker_sdk_path;
};
```

### RationalScan2ImConverter

New class in OpenMS inheriting from opentims's `Scan2InvIonMobilityConverter`.

**Physical model (ModelType=2):**

The calibration has two stages reflecting the TIMS device physics:

1. **Scan index → TIMS elution voltage:**
   ```
   V = c2 + ((c3 - c2) / c1) × (scan - c4 - c0)
   ```
   where c0+c4 is the scan offset, c1 is the scan span, and c2/c3 are voltage
   ramp endpoints.

2. **Voltage → inverse reduced ion mobility:**
   ```
   1/K0 = 1 / (c6 + c7 / V)
   ```
   where c6 and c7 are empirical voltage-to-mobility calibration constants,
   determined during instrument calibration with reference ions.

Combined:
```
1/K0 = 1 / (c6 + c7 / (c2 + ((c3 - c2) / c1) × (scan - c4 - c0)))
```

**Analytic inverse** (for `inverse_convert`):
```
scan = c0 + c4 + (c1 / (c3 - c2)) × (c7 / (1/y - c6) - c2)
```
where `y = 1/K0`. Singularity guards needed for:
- `c3 == c2` (division by zero in scan-to-voltage slope)
- `y == 1/c6` (division by zero when voltage term vanishes)
- `V == 0` in the forward direction (when the inner affine term evaluates to zero;
  physically unreachable for valid scan indices since c2/c3 are voltage endpoints
  typically in the hundreds, but defensive code should guard against it)

**Rounding in `inverse_convert`:** The method returns `uint32_t* scans`. Following
the convention in opentims's `OpenSourceScan2ImConverter`, results are rounded via
`static_cast<uint32_t>(val + 0.5)` with negative values clamped to 0.

**Coefficients C0–C9:** Stored in `TimsCalibration` table. C5, C8, C9 are present
in the table but unused in the only known model type (ModelType=2). No public
documentation assigns them meaning. They are stored in the struct for forward
compatibility but not referenced in the formula.

**Per-frame calibration:** The `Frames` table has a `TimsCalibration` column
(integer foreign key → `TimsCalibration.Id`). Most datasets have a single calibration
row, but the schema supports multiple segments (e.g., mid-run recalibration). The
converter stores a `vector<uint32_t>` mapping frame ID → calibration ID and an
`unordered_map<uint32_t, Coefficients>` for the calibration data. This is the first
open-source implementation to use per-frame calibration — all other tools either
ignore frame ID (linear model) or delegate to the Bruker SDK.

**ModelType handling:** Only ModelType=2 is supported (the only value observed in
public datasets). Unknown model types cause a warning log and fallback to the linear
converter.

**Thread safety:** The converter is immutable after construction — all calibration
data is loaded in the constructor and never modified. This makes `convert()` and
`inverse_convert()` safe for concurrent calls, which matters because opentims's
`TimsFrame::save_to_buffs()` calls the converter and may be invoked from multiple
threads.

**Destructor:** No custom destructor needed. The class holds only value types and
STL containers (`unordered_map`, `vector`), which are cleaned up by the default
destructor. The base class `Scan2InvIonMobilityConverter` has a virtual destructor,
so polymorphic deletion via `unique_ptr<Scan2InvIonMobilityConverter>` is safe.

```cpp
class RationalScan2ImConverter : public Scan2InvIonMobilityConverter
{
public:
  struct Coefficients
  {
    double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9;
    // c5, c8, c9: stored but unused in ModelType=2. No public documentation
    // assigns them meaning. Retained for forward compatibility.
  };

  /// @param calibrations  Map of calibration ID → coefficients (from TimsCalibration table)
  /// @param frame_to_cal  Calibration ID for each frame, indexed by frame_id (1-based)
  RationalScan2ImConverter(
    std::unordered_map<uint32_t, Coefficients> calibrations,
    std::vector<uint32_t> frame_to_cal);

  void convert(uint32_t frame_id, double* inv_ion_mobilities,
               const double* scans, uint32_t size) override;
  void convert(uint32_t frame_id, double* inv_ion_mobilities,
               const uint32_t* scans, uint32_t size) override;
  void inverse_convert(uint32_t frame_id, uint32_t* scans,
                       const double* inv_ion_mobilities, uint32_t size) override;
  /// Returns e.g. "RationalScan2ImConverter (2 calibration segments)"
  std::string description() const override;

private:
  std::unordered_map<uint32_t, Coefficients> calibrations_;
  std::vector<uint32_t> frame_to_cal_;  // indexed by frame_id (1-based, index 0 unused)

  /// Look up calibration for a frame. Bounds checks:
  /// - frame_id == 0: uses first calibration (should not occur in valid data)
  /// - frame_id > frame_to_cal_.size(): uses first calibration + OPENMS_LOG_WARN
  /// - NULL TimsCalibration FK in DB: handled at construction time by
  ///   tryCreateRationalConverter(), which falls back to linear if any frame
  ///   has a missing or invalid FK.
  const Coefficients& getCalibration(uint32_t frame_id) const;

  /// V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
  /// 1/K0 = 1.0 / (c6 + c7 / V)
  static double applyFormula(const Coefficients& c, double scan);

  /// scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1.0/inv_k0 - c6) - c2)
  static double invertFormula(const Coefficients& c, double inv_k0);
};
```

### Tiered Fallback in `openTimsDataHandle()`

The existing file-local helper is updated with a new signature:

```cpp
// Was: static std::unique_ptr<TimsDataHandle> openTimsDataHandle(const String& path)
// Now:
static std::unique_ptr<TimsDataHandle> openTimsDataHandle(
    const String& path, const Config& config)
```

**Call sites to update:** `load()` (line ~644), `transform()` (line ~690), and any
overloads that call through to `openTimsDataHandle`. Overloads without a `Config`
parameter pass a default-constructed `Config{}`.

The fallback chain:

1. **Create handle with linear converters** (safe baseline, always succeeds):
   ```cpp
   auto& tof_factory = OpenSourceTof2MzConverterFactory::instance();
   auto& im_factory = OpenSourceScan2ImConverterFactory::instance();
   auto handle = std::make_unique<TimsDataHandle>(path_string,
       NoPressureCompensation, &tof_factory, &im_factory);
   ```

2. **If strategy is AUTO or BRUKER_SDK:** try loading Bruker SDK from
   `config.bruker_sdk_path` or `OPENMS_BRUKER_SDK_PATH` env var. If found,
   map the pressure compensation enum and replace the converter:
   ```cpp
   // Map OpenMS enum → opentims enum
   // PressureCompensation::NONE   → NoPressureCompensation (0)
   // PressureCompensation::GLOBAL → AnalyisGlobalPressureCompensation (1)
   // PressureCompensation::PER_FRAME → PerFramePressureCompensationWithMissingReference (3)
   //   (uses WithMissingReference variant — more robust, does not fail if a
   //    reference frame is missing. The non-missing variant (2) is not exposed.)
   auto pcs = mapPressureCompensation(config.pressure_compensation);

   // BrukerScan2InvIonMobilityConverter needs the TimsDataHandle& and lib path
   handle->scan2inv_ion_mobility_converter =
       BrukerScan2InvIonMobilityConverterFactory::instance(sdk_path)
           .produce(*handle, pcs);
   ```
   If BRUKER_SDK was explicitly requested but SDK not found, throw
   `Exception::FileNotReadable`.

3. **If strategy is AUTO or RATIONAL:** try reading `TimsCalibration` table from
   SQLite. If present with all ModelType=2 rows, create `RationalScan2ImConverter`
   and replace the converter via direct assignment to the public `unique_ptr`:
   ```cpp
   auto converter = tryCreateRationalConverter(path);
   if (converter)
   {
     handle->scan2inv_ion_mobility_converter = std::move(converter);
   }
   ```
   If RATIONAL was explicitly requested but table missing or unknown ModelType,
   throw `Exception::FileNotReadable`.

4. **Fall through to linear** (already set from step 1).

The `TimsDataHandle::scan2inv_ion_mobility_converter` is a **public**
`std::unique_ptr` member (opentims.h line 212). Direct `std::move` assignment
replaces the converter. Do NOT use the private `set_converter()` method (line 231).

The TOF-to-m/z converter is unchanged in all cases (stays linear-in-sqrt with
optional OLS recalibration for DDA).

### Logging

All fallback transitions are logged to aid debugging:

| Event | Level | Example message |
|---|---|---|
| SDK converter selected | INFO | `"TIMS calibration: Bruker SDK (pressure_comp=global)"` |
| Rational converter selected | INFO | `"TIMS calibration: rational (TimsCalibration table, 1 segment)"` |
| Linear fallback (no table) | INFO | `"TIMS calibration: linear (GlobalMetadata)"` |
| SDK not found during AUTO | DEBUG | `"Bruker SDK not found at ..., trying rational"` |
| Unknown ModelType in table | WARN | `"TimsCalibration ModelType=3 unsupported, falling back to linear"` |
| Pressure comp requested without SDK | WARN | `"Pressure compensation requires BRUKER_SDK strategy, ignoring"` |
| Multiple calibration segments | INFO | `"TimsCalibration: 3 calibration segments for N frames"` |

### SQLite Queries

**New query — TimsCalibration table:**
```sql
SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9
FROM TimsCalibration
```

**Modified query — Frames table** (add TimsCalibration column):
```sql
SELECT Id, ..., TimsCalibration FROM Frames
```

Both queries are wrapped in try/catch so that missing tables gracefully fall back
to the linear converter. Note: the `TimsCalibration` column in the `Frames` table
may also be absent in very old TDF versions. The Frames query should be a separate
try/catch from the main Frames metadata query so that a missing column does not
break existing frame loading. If the column is absent, treat it the same as a
missing `TimsCalibration` table (fall back to linear).

### Limitations vs Bruker SDK

The rational function model does NOT replicate:

- **Pressure compensation:** Correcting for ambient gas pressure drift during
  acquisition. This is SDK-only (available when BRUKER_SDK strategy is active).
- **Recalibrated-state awareness:** The SDK can apply post-acquisition recalibration.
  The rational model uses raw calibration coefficients from acquisition time.
- **Unknown ModelType variants:** Only ModelType=2 is implemented. Other types
  (if they exist) fall back to linear.

## File Layout

**New files:**
- `src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h`
- `src/openms/source/FORMAT/RationalScan2ImConverter.cpp`

**Modified files:**
- `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` — Config additions
- `src/openms/source/FORMAT/BrukerTimsFile.cpp` — tiered `openTimsDataHandle()`

**No changes to:**
- opentims (FetchContent dependency, untouched)
- CMake build (new files compile as part of OpenMS, opentims headers already available)
- Downstream consumers (FileHandler, search engines) — backward compatible

## Backward Compatibility

`Config` defaults to `AUTO` strategy, `NONE` pressure compensation, empty SDK path.
With no env var set and a TDF file containing a `TimsCalibration` table (virtually
all files), the behavior silently upgrades from linear → rational. For files without
the table, it stays linear. No API break.

## Testing

- **Unit:** `RationalScan2ImConverter` with known coefficients → expected 1/K0 values.
  Round-trip test via `inverse_convert`. Singularity edge cases.
- **Unit:** `tryCreateRationalConverter` with mock SQLite — table present, table
  missing, unknown ModelType, missing `TimsCalibration` column in Frames.
- **Integration:** Existing Bruker test data (if available) should produce results
  close to the current linear path for simple single-calibration datasets, with
  improved accuracy at scan range edges.

**Worked validation example** (from public opentims test data, ModelType=2):

```
C0=1, C1=917, C2=213.5998, C3=75.81729, C4=33, C5=1,
C6=-0.009065829, C7=135.4364, C8=13.32608, C9=1663.341

For scan=500:
  V = 213.5998 + ((75.81729 - 213.5998) / 917) * (500 - 33 - 1)
    = 213.5998 + (-137.78251 / 917) * 466
    = 213.5998 + (-0.150255) * 466
    = 213.5998 - 70.019
    = 143.581

  1/K0 = 1 / (-0.009065829 + 135.4364 / 143.581)
       = 1 / (-0.009065829 + 0.94325)
       = 1 / 0.93419
       = 1.07045

Round-trip: invertFormula(coeffs, 1.07045) should return ≈ 500.0
```

This example should be included as a unit test assertion with appropriate floating
point tolerance (~1e-6).

## References

- diapysef rational formula: https://github.com/Roestlab/dia-pasef/blob/master/src/diapysef/diapysef/timsdata.py
- OpenTIMS converter architecture: https://github.com/michalsta/opentims
- TIMSCONVERT pressure compensation: https://github.com/gtluu/timsconvert
- Bruker tdf-sdk release notes (pressure comp in v2.21.0, CCS in v2.7.0)
- TIMS calibration background: https://pmc.ncbi.nlm.nih.gov/articles/PMC6618043/
