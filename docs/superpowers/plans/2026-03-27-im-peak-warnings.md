# IM_PEAK Format Warnings for TOPP Tools

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add IM_PEAK format detection and warnings/errors to four TOPP tools that cannot properly handle per-peak ion mobility data.

**Architecture:** After loading the MSExperiment, check for IM_PEAK format using `IMTypes::determineIMFormat()`. If detected, report IMPeakType in the message. Three tools error out (`INCOMPATIBLE_INPUT_DATA`), one warns and continues.

**Tech Stack:** C++, OpenMS IMTypes/IMPeakType API, TOPP tool framework

---

### Task 1: Resampler — error on IM_PEAK

**Files:**
- Modify: `src/topp/Resampler.cpp`

- [ ] **Step 1: Add include**

Add after line 12 (`#include <OpenMS/APPLICATIONS/TOPPBase.h>`):

```cpp
#include <OpenMS/IONMOBILITY/IMTypes.h>
```

- [ ] **Step 2: Add IM_PEAK check after loading**

Insert after line 102 (`FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);`):

```cpp
    // Check for unsupported per-peak ion mobility data
    for (const auto& spec : exp)
    {
      IMFormat im_format = IMTypes::determineIMFormat(spec);
      if (im_format == IMFormat::IM_PEAK)
      {
        OPENMS_LOG_ERROR << "Error: Input contains per-peak ion mobility data (IM_PEAK, "
                         << imPeakTypeToString(spec.getIMPeakType())
                         << ") which is not supported by Resampler. "
                         << "Preprocess with IonMobilityBinning or PeakPickerIM first." << std::endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
```

- [ ] **Step 3: Build and verify existing tests still pass**

Run: `cmake --build OpenMS-build --target Resampler -j$(nproc)`
Run: `ctest --test-dir OpenMS-build -R TOPP_Resampler -V`
Expected: All existing Resampler tests pass (non-IM data is unaffected).

- [ ] **Step 4: Commit**

```bash
git add src/topp/Resampler.cpp
git commit -m "feat(Resampler): error on unsupported IM_PEAK input data"
```

---

### Task 2: FeatureFinderCentroided — error on IM_PEAK

**Files:**
- Modify: `src/topp/FeatureFinderCentroided.cpp`

- [ ] **Step 1: Add include**

Add after line 17 (`#include <OpenMS/IONMOBILITY/IMDataConverter.h>`):

```cpp
#include <OpenMS/IONMOBILITY/IMTypes.h>
```

- [ ] **Step 2: Add IM_PEAK check after loading**

Insert after line 191 (the empty spectra check block, before the `// determine type of spectral data` comment on line 193). The check goes after `exp.updateRanges()` and the empty check, but before FAIMS splitting:

```cpp
    // Check for unsupported per-peak ion mobility data
    for (const auto& spec : exp)
    {
      IMFormat im_format = IMTypes::determineIMFormat(spec);
      if (im_format == IMFormat::IM_PEAK)
      {
        OPENMS_LOG_ERROR << "Error: Input contains per-peak ion mobility data (IM_PEAK, "
                         << imPeakTypeToString(spec.getIMPeakType())
                         << ") which is not supported by FeatureFinderCentroided. "
                         << "Preprocess with IonMobilityBinning or PeakPickerIM first." << std::endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
```

- [ ] **Step 3: Build and verify existing tests still pass**

Run: `cmake --build OpenMS-build --target FeatureFinderCentroided -j$(nproc)`
Run: `ctest --test-dir OpenMS-build -R TOPP_FeatureFinderCentroided -V`
Expected: All existing tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/topp/FeatureFinderCentroided.cpp
git commit -m "feat(FeatureFinderCentroided): error on unsupported IM_PEAK input data"
```

---

### Task 3: FeatureFinderMultiplex — error on IM_PEAK

**Files:**
- Modify: `src/topp/FeatureFinderMultiplex.cpp`

Note: This file already includes `<OpenMS/IONMOBILITY/IMTypes.h>` (line 17).

- [ ] **Step 1: Add IM_PEAK check after loading**

Insert after line 256 (`file.loadExperiment(in_, exp, {FileTypes::MZML}, log_type_);`), before the `// Prepare algorithm parameters` comment on line 258:

```cpp
    // Check for unsupported per-peak ion mobility data
    for (const auto& spec : exp)
    {
      IMFormat im_format = IMTypes::determineIMFormat(spec);
      if (im_format == IMFormat::IM_PEAK)
      {
        OPENMS_LOG_ERROR << "Error: Input contains per-peak ion mobility data (IM_PEAK, "
                         << imPeakTypeToString(spec.getIMPeakType())
                         << ") which is not supported by FeatureFinderMultiplex. "
                         << "Preprocess with IonMobilityBinning or PeakPickerIM first." << std::endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
```

- [ ] **Step 2: Build and verify existing tests still pass**

Run: `cmake --build OpenMS-build --target FeatureFinderMultiplex -j$(nproc)`
Run: `ctest --test-dir OpenMS-build -R TOPP_FeatureFinderMultiplex -V`
Expected: All existing tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/topp/FeatureFinderMultiplex.cpp
git commit -m "feat(FeatureFinderMultiplex): error on unsupported IM_PEAK input data"
```

---

### Task 4: PeakPickerHiRes — warn on IM_PEAK

**Files:**
- Modify: `src/topp/PeakPickerHiRes.cpp`

- [ ] **Step 1: Add include**

Add after line 14 (`#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>`):

```cpp
#include <OpenMS/IONMOBILITY/IMTypes.h>
```

- [ ] **Step 2: Add IM_PEAK warning after loading (in-memory path)**

Insert after line 212 (`FileHandler().loadExperiment(in, ms_exp_raw, {FileTypes::MZML}, log_type_);`), before the empty check on line 214:

```cpp
    // Warn about per-peak ion mobility data (PeakPickerHiRes picks m/z only)
    for (const auto& spec : ms_exp_raw)
    {
      IMFormat im_format = IMTypes::determineIMFormat(spec);
      if (im_format == IMFormat::IM_PEAK)
      {
        OPENMS_LOG_WARN << "Warning: Input contains per-peak ion mobility data (IM_PEAK, "
                        << imPeakTypeToString(spec.getIMPeakType())
                        << "). PeakPickerHiRes picks in m/z only and reports intensity-weighted "
                        << "mean ion mobility. This produces incorrect results on unbinned data. "
                        << "Consider IonMobilityBinning or PeakPickerIM first." << std::endl;
        break; // warn once
      }
    }
```

- [ ] **Step 3: Build and verify existing tests still pass**

Run: `cmake --build OpenMS-build --target PeakPickerHiRes -j$(nproc)`
Run: `ctest --test-dir OpenMS-build -R TOPP_PeakPickerHiRes -V`
Expected: All existing tests pass (non-IM data triggers no warning).

- [ ] **Step 4: Commit**

```bash
git add src/topp/PeakPickerHiRes.cpp
git commit -m "feat(PeakPickerHiRes): warn on IM_PEAK input data (m/z-only picking)"
```
