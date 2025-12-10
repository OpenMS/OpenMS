# MSSpectrum AoS to SoA Refactoring Plan

## Executive Summary

This document outlines a comprehensive plan to refactor `MSSpectrum` from an Array-of-Structures (AoS) layout using `Peak1D` to a Structure-of-Arrays (SoA) layout with separate m/z, intensity, and optional ion mobility arrays.

**Current State:**
```cpp
class MSSpectrum : private std::vector<Peak1D>, ... {
    // Peak1D = { DPosition<1> position_, float intensity_ }
    FloatDataArrays float_data_arrays_;  // IM hidden here
};
```

**Target State:**
```cpp
class MSSpectrum : public RangeManagerContainer<...>, public SpectrumSettings {
    // SoA storage with shared_ptr for Python zero-copy safety (see Appendix D)
    std::shared_ptr<std::vector<double>> mz_;
    std::shared_ptr<std::vector<float>> intensity_;
    std::shared_ptr<std::vector<float>> ion_mobility_;  // First-class citizen (optional)
    DriftTimeUnit ion_mobility_unit_;
    FloatDataArrays float_data_arrays_;    // Other metadata only
};
```

**Key Design Decision:** Arrays are wrapped in `shared_ptr` to enable safe zero-copy Python/NumPy views (see [Appendix D](#appendix-d-shared_ptr-wrapping-for-data-arrays)).

**Scope:** This refactoring targets `MSSpectrum` only. `MSChromatogram` (which has a similar AoS pattern with `ChromatogramPeak`) is **out of scope** for this initial refactoring but can follow the same pattern in a subsequent effort.

---

## Table of Contents

1. [Motivation and Goals](#1-motivation-and-goals)
2. [Why Not Use Proxy Iterators](#2-why-not-use-proxy-iterators)
3. [Proposed API Design](#3-proposed-api-design)
4. [Ion Mobility Integration](#4-ion-mobility-integration)
5. [Impact Analysis](#5-impact-analysis)
6. [Migration Strategy](#6-migration-strategy)
7. [Implementation Phases](#7-implementation-phases)
8. [Testing Strategy](#8-testing-strategy)
9. [pyOpenMS Bindings](#9-pyopenms-bindings)
10. [File I/O Considerations](#10-file-io-considerations)
11. [TOPP Tools Migration](#11-topp-tools-migration)
12. [Risk Assessment](#12-risk-assessment)
13. [Success Metrics](#13-success-metrics)
14. [Zero-Copy NumPy Integration](#14-zero-copy-numpy-integration)
15. [Appendix A: File List](#appendix-a-file-list)
16. [Appendix B: API Quick Reference](#appendix-b-api-quick-reference)
17. [Appendix C: Future Extensibility](#appendix-c-future-extensibility---adding-new-dimensions)
18. [Appendix D: shared_ptr Wrapping](#appendix-d-shared_ptr-wrapping-for-data-arrays)

---

## 1. Motivation and Goals

### 1.1 Performance Benefits

| Operation | AoS (Current) | SoA (Target) | Expected Improvement |
|-----------|---------------|--------------|----------------------|
| Sequential m/z scan | Cache misses (16-byte stride) | Cache-friendly (8-byte stride) | 2-4x |
| Sequential intensity scan | Cache misses | Cache-friendly (4-byte stride) | 2-4x |
| TIC calculation | Scattered access | SIMD vectorizable | 4-8x |
| Binary search on m/z | Good | Excellent (contiguous doubles) | 1.5-2x |
| Sorting | Index shuffle + Peak1D moves | Index shuffle + array permutation | 1.2-1.5x |

### 1.2 mzML Format Alignment

**Critical insight:** The mzML file format **already stores data as SoA** - m/z and intensity are separate `binaryDataArray` elements.

**Current inefficient data flow:**
```
mzML file (SoA) → MzMLHandler → Peak1D vector (AoS) → pyOpenMS → NumPy (SoA)
                              ↑                      ↑
                        Converts to AoS        Converts back to SoA
```

**After refactoring:**
```
mzML file (SoA) → MzMLHandler → SoA arrays → pyOpenMS → NumPy (SoA)
                              ↑                      ↑
                        Direct copy            Zero-copy view!
```

This eliminates **two unnecessary format conversions**, providing significant I/O performance improvement.

### 1.3 Design Goals

1. **Performance**: Enable SIMD vectorization and improve cache utilization
2. **File Format Alignment**: Match mzML's native SoA layout to eliminate conversion overhead
3. **Ion Mobility First-Class**: Elevate IM from hidden FloatDataArray to dedicated member
4. **API Clarity**: Explicit, honest API that doesn't pretend to be something it's not
5. **Gradual Migration**: Provide compatibility layer for existing code
6. **pyOpenMS Improvement**: Simplify Python bindings with zero-copy capability

### 1.4 Non-Goals

- Maintaining proxy/iterator facade that mimics `std::vector<Peak1D>`
- Binary compatibility with existing compiled code
- Zero migration effort for downstream code

### 1.5 Target Version and Timeline

- **Target Release:** OpenMS 4.0 (major version bump due to API breaking changes)
- **Deprecation Period:** Compatibility methods available for 2 releases before removal
- **MSChromatogram:** Out of scope for this refactoring; will follow in subsequent effort

---

## 2. Why Not Use Proxy Iterators

A proxy-based approach was considered but rejected due to fundamental C++ semantic issues:

### 2.1 Showstopper: Pointer Arithmetic in Existing Code

`MSSpectrum.h:344-351` contains:
```cpp
template<class Predicate>
bool isSorted(const Predicate& lambda) const {
    auto value_2_index_wrapper = [this, &lambda](const value_type& value1, const value_type& value2) {
        const Size index1 = (&value1) - (&this->front());  // POINTER ARITHMETIC!
        const Size index2 = (&value2) - (&this->front());
        return lambda(index1, index2);
    };
    return std::is_sorted(this->begin(), this->end(), value_2_index_wrapper);
}
```

This code explicitly relies on contiguous `Peak1D` storage.

### 2.2 Exposed `data()` Method

`MSSpectrum.h:148` exposes:
```cpp
using ContainerType::data;  // Returns Peak1D*
```

Any code using `data()` for direct memory access would break.

### 2.3 Reference Binding Issues

Found throughout codebase:
```cpp
for (Peak1D& pk : ms) {           // Won't compile with proxy
    pk.setIntensity(pk.getIntensity() * factor);
}
```

### 2.4 STL Algorithm Compatibility

`std::swap(*it1, *it2)` would swap proxy objects, not underlying data.

### 2.5 Conclusion

A clean, explicit API is preferable to a leaky abstraction that creates subtle bugs.

---

## 3. Proposed API Design

### 3.1 Core Data Structure

```cpp
class MSSpectrum final :
    public RangeManagerContainer<RangeMZ, RangeIntensity, RangeMobility>,
    public SpectrumSettings
{
private:
    // Core SoA arrays with shared_ptr for Python zero-copy safety
    // Always initialized to non-null empty vectors
    std::shared_ptr<std::vector<double>> mz_ = std::make_shared<std::vector<double>>();
    std::shared_ptr<std::vector<float>> intensity_ = std::make_shared<std::vector<float>>();

    // Optional ion mobility (per-peak, for CONCATENATED/CENTROIDED format)
    // nullptr when not present
    std::shared_ptr<std::vector<float>> ion_mobility_;
    DriftTimeUnit ion_mobility_unit_ = DriftTimeUnit::NONE;

    // Per-spectrum metadata (existing)
    double retention_time_ = -1;
    double drift_time_ = IMTypes::DRIFTTIME_NOT_SET;  // For MULTIPLE_SPECTRA format
    DriftTimeUnit drift_time_unit_ = DriftTimeUnit::NONE;
    UInt ms_level_ = 1;
    String name_;

    // Additional metadata arrays (arbitrary annotations, NOT ion mobility)
    FloatDataArrays float_data_arrays_;
    StringDataArrays string_data_arrays_;
    IntegerDataArrays integer_data_arrays_;

    // Copy-on-Write helper
    void ensureUniqueMZ() {
        if (mz_.use_count() > 1) {
            mz_ = std::make_shared<std::vector<double>>(*mz_);
        }
    }
    void ensureUniqueIntensity() {
        if (intensity_.use_count() > 1) {
            intensity_ = std::make_shared<std::vector<float>>(*intensity_);
        }
    }
};

// Verify Peak1D size assumption for documentation accuracy
static_assert(sizeof(Peak1D) == 16, "Peak1D size assumption violated");
```

### 3.2 Primary API (New - Preferred)

```cpp
// ═══════════════════════════════════════════════════════════════
// DIRECT ARRAY ACCESS - Preferred for performance-critical code
// ═══════════════════════════════════════════════════════════════

// M/Z array (always present, never null)
const std::vector<double>& getMZArray() const noexcept { return *mz_; }
std::vector<double>& getMZArray() noexcept { ensureUniqueMZ(); return *mz_; }

// Intensity array (always present, never null)
const std::vector<float>& getIntensityArray() const noexcept { return *intensity_; }
std::vector<float>& getIntensityArray() noexcept { ensureUniqueIntensity(); return *intensity_; }

// Ion mobility array (optional)
bool hasIonMobilityArray() const noexcept;
const std::vector<float>& getIonMobilityArray() const;  // throws if empty
std::vector<float>& getIonMobilityArray();
DriftTimeUnit getIonMobilityUnit() const noexcept;

// ═══════════════════════════════════════════════════════════════
// SHARED_PTR ACCESS - For Python/Cython zero-copy bindings
// ═══════════════════════════════════════════════════════════════

// Returns shared_ptr for safe zero-copy Python views
std::shared_ptr<const std::vector<double>> getMZArrayPtr() const noexcept { return mz_; }
std::shared_ptr<const std::vector<float>> getIntensityArrayPtr() const noexcept { return intensity_; }
std::shared_ptr<const std::vector<float>> getIonMobilityArrayPtr() const;  // nullptr if no IM

// ═══════════════════════════════════════════════════════════════
// INDEXED ACCESS - Convenience methods
// ═══════════════════════════════════════════════════════════════

double getMZ(Size i) const;
float getIntensity(Size i) const;
void setMZ(Size i, double mz);
void setIntensity(Size i, float intensity);

float getIonMobility(Size i) const;  // requires hasIonMobilityArray()
void setIonMobility(Size i, float im);
std::optional<float> maybeGetIonMobility(Size i) const;

// ═══════════════════════════════════════════════════════════════
// BULK OPERATIONS
// ═══════════════════════════════════════════════════════════════

void setData(std::vector<double> mz, std::vector<float> intensity);
void setData(std::vector<double> mz, std::vector<float> intensity,
             std::vector<float> ion_mobility, DriftTimeUnit im_unit);

void setIonMobilityArray(std::vector<float> im, DriftTimeUnit unit);
void clearIonMobility();
void initializeIonMobility(DriftTimeUnit unit);  // Creates empty array

// ═══════════════════════════════════════════════════════════════
// ADDING PEAKS
// ═══════════════════════════════════════════════════════════════

void addPeak(double mz, float intensity);
void addPeak(double mz, float intensity, float ion_mobility);

// ═══════════════════════════════════════════════════════════════
// CONTAINER OPERATIONS
// ═══════════════════════════════════════════════════════════════

Size size() const noexcept;
bool empty() const noexcept;
void reserve(Size n);
void resize(Size n);
void clear(bool clear_meta_data = false);
MSSpectrum& select(const std::vector<Size>& indices);
```

### 3.3 Compatibility API (Deprecated)

```cpp
// ═══════════════════════════════════════════════════════════════
// DEPRECATED - For gradual migration only
// ═══════════════════════════════════════════════════════════════

[[deprecated("Use getMZ(i)/getIntensity(i) for better performance")]]
Peak1D getPeak(Size i) const;

[[deprecated("Use setMZ(i)/setIntensity(i) for better performance")]]
void setPeak(Size i, const Peak1D& peak);

[[deprecated("Consider using array access instead")]]
std::vector<Peak1D> toAoS() const;

void fromAoS(const std::vector<Peak1D>& peaks);

// Legacy IM methods
[[deprecated("Use hasIonMobilityArray() instead")]]
bool containsIMData() const;

[[deprecated("Use getIonMobilityArray() and getIonMobilityUnit()")]]
std::pair<Size, DriftTimeUnit> getIMData() const;
```

### 3.4 Removed/Changed Methods

| Old Method | Status | Replacement |
|------------|--------|-------------|
| `operator[]` returning `Peak1D&` | Removed | `getMZ(i)`, `getIntensity(i)` |
| `begin()`, `end()` (Peak1D iterators) | Removed | Index-based loops or array access |
| `data()` returning `Peak1D*` | Removed | `getMZArray().data()`, `getIntensityArray().data()` |
| `push_back(Peak1D)` | Deprecated | `addPeak(mz, intensity)` |
| `emplace_back(mz, intensity)` | Changed | `addPeak(mz, intensity)` |

---

## 4. Ion Mobility Integration

### 4.1 Current State Problems

1. **Hidden in FloatDataArrays**: IM array identified by CV term name lookup
2. **O(n) Lookup**: `getIMData()` scans all float arrays
3. **No Type Safety**: IM is just another generic float array
4. **Synchronization Risk**: Easy to accidentally break peak/IM alignment

### 4.2 New Design

Ion mobility becomes a **first-class citizen**:

```cpp
// Dedicated members
std::vector<float> ion_mobility_;
DriftTimeUnit ion_mobility_unit_ = DriftTimeUnit::NONE;

// O(1) access
bool hasIonMobilityArray() const noexcept {
    return !ion_mobility_.empty();
}

const std::vector<float>& getIonMobilityArray() const {
    if (ion_mobility_.empty()) {
        throw Exception::MissingInformation(...);
    }
    return ion_mobility_;
}
```

### 4.3 IMFormat Detection

```cpp
IMFormat getIMFormat() const {
    bool has_array = hasIonMobilityArray();
    bool has_single = (drift_time_ != IMTypes::DRIFTTIME_NOT_SET);

    if (has_array && has_single) {
        throw Exception::InvalidValue("Both IM array and single drift time");
    }
    if (has_array) return IMFormat::CONCATENATED;  // or CENTROIDED
    if (has_single) return IMFormat::MULTIPLE_SPECTRA;
    return IMFormat::NONE;
}
```

### 4.4 Sorting with Ion Mobility

All three core arrays are permuted together:

```cpp
void MSSpectrum::sortByPosition() {
    if (isSorted()) return;

    std::vector<Size> indices(size());
    std::iota(indices.begin(), indices.end(), 0);
    std::stable_sort(indices.begin(), indices.end(),
        [this](Size a, Size b) { return mz_[a] < mz_[b]; });

    applyPermutation(indices);
}

void MSSpectrum::applyPermutation(const std::vector<Size>& indices) {
    mz_ = permute(mz_, indices);
    intensity_ = permute(intensity_, indices);
    if (hasIonMobilityArray()) {
        ion_mobility_ = permute(ion_mobility_, indices);
    }
    // Also permute float_data_arrays_, string_data_arrays_, etc.
}
```

---

## 5. Impact Analysis

### 5.1 Library Code (src/openms/)

**Files referencing MSSpectrum or Peak1D: ~121 files** (verified via grep)

| Access Pattern | Files | Migration Effort |
|----------------|-------|------------------|
| Direct Peak1D manipulation (`push_back`, `emplace_back`) | ~8 | High - requires API change |
| Indexed access (`spec[i].getMZ()`) | ~29 | Medium - change to `spec.getMZ(i)` |
| Range-based for loops (`for (Peak1D& p : spec)`) | ~9 | Medium - change to index loops |
| Read-only/metadata access | ~75 | Low - may not need changes |

| Category | Files | Effort |
|----------|-------|--------|
| ANALYSIS/ | ~35 | High |
| FORMAT/ | ~25 | High (I/O critical) |
| PROCESSING/ | ~20 | Medium |
| FEATUREFINDER/ | ~12 | High (performance) |
| KERNEL/ | ~8 | Core changes |
| Other | ~21 | Medium |

**Note:** Many files only reference MSSpectrum for type declarations or metadata access and won't require significant changes.

### 5.2 TOPP Tools (src/topp/)

**19 tools use MSSpectrum directly:**

| Tool | Peak Access Pattern | IM Usage | Priority |
|------|---------------------|----------|----------|
| MapNormalizer | `for (Peak1D& pk : ms)` | No | High |
| AssayGeneratorMetabo | Range-based + modification | No | High |
| PhosphoScoring | `findNearest()` + indexed | No | High |
| PeakPickerIM | IM-specialized | Yes | Critical |
| FileFilter | Filter operations | Yes | Medium |
| FileInfo | Statistics | Yes | Low |
| CometAdapter | Metadata only | Yes | Low |
| Others (12) | Various | Mixed | Medium |

### 5.3 Test Files

**Primary test files requiring updates:**

| File | Focus | Lines |
|------|-------|-------|
| `MSSpectrum_test.cpp` | Core functionality | ~800 |
| `Peak1D_test.cpp` | Peak type | ~350 |
| `MSExperiment_test.cpp` | Container | ~500 |
| `PeakPickerHiRes_test.cpp` | IM integration | ~300 |
| `PeakPickerIM_test.cpp` | IM-specific | ~400 |
| `IMDataConverter_test.cpp` | IM conversion | ~250 |

### 5.4 pyOpenMS Bindings

**Current binding patterns (src/pyOpenMS/):**

```python
# Current implementation in MSSpectrum.pyx
def get_peaks(self):
    for i in range(n):
        mzs[i] = deref(spec_)[i].getMZ()       # Iterates Peak1D
        intensities[i] = deref(spec_)[i].getIntensity()
    return mzs, intensities

def set_peaks(self, peaks):
    for i in range(N):
        p.setMZ(<double>mz)
        p.setIntensity(<float>intensity)
        spec_.push_back(p)                     # Creates Peak1D
```

**After SoA, bindings become simpler and faster:**

```python
# New implementation
def get_peaks(self):
    # Direct array copy - no Peak1D iteration
    return np.array(spec_.getMZArray()), np.array(spec_.getIntensityArray())

def set_peaks(self, peaks):
    spec_.setData(mzs, intensities)            # Bulk operation
```

---

## 6. Migration Strategy

### 6.1 Overview

```
Phase 1: Add SoA API     Phase 2: Deprecate AoS    Phase 3: Remove AoS
     │                         │                         │
     ▼                         ▼                         ▼
┌─────────────┐          ┌─────────────┐          ┌─────────────┐
│ SoA + AoS   │          │ SoA primary │          │ SoA only    │
│ both work   │    →     │ AoS warns   │    →     │ clean API   │
└─────────────┘          └─────────────┘          └─────────────┘
```

### 6.2 Code Migration Patterns

**Pattern 1: Indexed access**
```cpp
// Before
for (Size i = 0; i < spec.size(); ++i) {
    double mz = spec[i].getMZ();
    float intensity = spec[i].getIntensity();
}

// After
const auto& mz_arr = spec.getMZArray();
const auto& int_arr = spec.getIntensityArray();
for (Size i = 0; i < spec.size(); ++i) {
    double mz = mz_arr[i];
    float intensity = int_arr[i];
}

// Or simpler
for (Size i = 0; i < spec.size(); ++i) {
    double mz = spec.getMZ(i);
    float intensity = spec.getIntensity(i);
}
```

**Pattern 2: Range-based for (reading)**
```cpp
// Before
for (const Peak1D& peak : spectrum) {
    process(peak.getMZ(), peak.getIntensity());
}

// After
const auto& mz = spectrum.getMZArray();
const auto& intensity = spectrum.getIntensityArray();
for (Size i = 0; i < spectrum.size(); ++i) {
    process(mz[i], intensity[i]);
}
```

**Pattern 3: Range-based for (modifying)**
```cpp
// Before
for (Peak1D& peak : spectrum) {
    peak.setIntensity(peak.getIntensity() * factor);
}

// After
auto& intensity = spectrum.getIntensityArray();
for (float& i : intensity) {
    i *= factor;
}
// Or SIMD-friendly:
std::transform(intensity.begin(), intensity.end(), intensity.begin(),
               [factor](float i) { return i * factor; });
```

**Pattern 4: Iterator-based search (MZBegin/MZEnd)**

**Important:** Current `MZBegin()`/`MZEnd()` return `Iterator` (Peak1D pointer). After refactoring, these methods will be **deprecated** in favor of index-based alternatives.

```cpp
// Before - returns Peak1D iterator
auto it = spectrum.MZBegin(500.0);
while (it != spectrum.MZEnd(600.0)) {
    process(it->getMZ(), it->getIntensity());
    ++it;
}

// After (Option A) - use new index-based methods
Size start = spectrum.findFirstMZ(500.0);
Size end = spectrum.findLastMZ(600.0);
const auto& mz = spectrum.getMZArray();
const auto& intensity = spectrum.getIntensityArray();
for (Size i = start; i < end; ++i) {
    process(mz[i], intensity[i]);
}

// After (Option B) - use range-based helper
auto [start, end] = spectrum.findMZRange(500.0, 600.0);
for (Size i = start; i < end; ++i) {
    process(spectrum.getMZ(i), spectrum.getIntensity(i));
}
```

**Migration note for MZBegin/MZEnd users:**

| Old Method | New Method | Return Type Change |
|------------|------------|-------------------|
| `MZBegin(mz)` | `findFirstMZ(mz)` | `Iterator` → `Size` |
| `MZEnd(mz)` | `findLastMZ(mz)` | `Iterator` → `Size` |
| `MZBegin(begin, mz, end)` | Use `std::lower_bound` on `getMZArray()` | Iterator → Size |

The old iterator-based methods will remain available (deprecated) during the transition period but will be removed in a future major release.

**Pattern 5: Adding peaks**
```cpp
// Before
Peak1D p;
p.setMZ(500.0);
p.setIntensity(1000.0f);
spectrum.push_back(p);

// After
spectrum.addPeak(500.0, 1000.0f);
```

**Pattern 6: Ion mobility access**
```cpp
// Before
if (spectrum.containsIMData()) {
    auto [im_idx, unit] = spectrum.getIMData();
    const auto& im_array = spectrum.getFloatDataArrays()[im_idx];
    for (Size i = 0; i < spectrum.size(); ++i) {
        float im = im_array[i];
    }
}

// After
if (spectrum.hasIonMobilityArray()) {
    const auto& im_array = spectrum.getIonMobilityArray();
    DriftTimeUnit unit = spectrum.getIonMobilityUnit();
    for (Size i = 0; i < spectrum.size(); ++i) {
        float im = im_array[i];
    }
}
```

---

## 7. Implementation Phases

> **Note:** Time estimates removed as per project guidelines. Each phase should be completed and tested before proceeding to the next.

### Phase 0: Preparation

**0.1 Baseline Measurements**
- [ ] Create performance benchmarks for current AoS implementation
- [ ] Document current memory usage per spectrum
- [ ] Record mzML load/save times for reference datasets
- [ ] Identify all code paths using `MZBegin()`/`MZEnd()` iterators

**0.2 Feature Branch Setup**
- [ ] Create feature branch `feature/msspectrum-soa`
- [ ] Set up CI to run full test suite on each commit
- [ ] Document rollback procedure

### Phase 1: Core Infrastructure (Non-Breaking)

**Strategy:** Add new SoA storage and API **alongside** existing AoS. Both work simultaneously. The AoS storage (`std::vector<Peak1D>` inheritance) remains authoritative during this phase.

**1.1 Add SoA storage (dual-storage mode)**
- [ ] Add `shared_ptr<vector<T>>` members: `mz_`, `intensity_`, `ion_mobility_`
- [ ] Keep `std::vector<Peak1D>` inheritance (remains authoritative)
- [ ] Implement **one-way sync**: mutations to Peak1D auto-update SoA arrays
- [ ] Add `static_assert(sizeof(Peak1D) == 16)` to verify assumptions

**1.2 Implement new API methods (read from SoA)**
- [ ] `getMZArray()`, `getIntensityArray()` - return `const` references to SoA
- [ ] `getMZArrayPtr()`, `getIntensityArrayPtr()` - return `shared_ptr` for Python
- [ ] `getMZ(i)`, `getIntensity(i)` - read from SoA arrays
- [ ] `hasIonMobilityArray()`, `getIonMobilityArray()`, `getIonMobilityUnit()`
- [ ] `addPeak()` overloads - write to both AoS and SoA

**1.3 Sync mechanism**
```cpp
// During Phase 1: AoS is authoritative, SoA is derived
void push_back(const Peak1D& p) {
    ContainerType::push_back(p);  // AoS
    mz_->push_back(p.getMZ());    // Sync to SoA
    intensity_->push_back(p.getIntensity());
}
```

**1.4 Testing**
- [ ] All existing tests pass unchanged
- [ ] New tests verify SoA arrays match AoS content
- [ ] Verify `getMZArray()` returns correct data

### Phase 2: Internal Migration

**Strategy:** Migrate internal MSSpectrum methods to use SoA storage directly. AoS becomes a facade.

**2.1 Update internal methods to use SoA**
- [ ] `sortByPosition()` - use index permutation on SoA arrays
- [ ] `sortByIntensity()` - same pattern
- [ ] `sortByIonMobility()` - same pattern
- [ ] `select()` - apply to all SoA arrays
- [ ] `updateRanges()` - iterate SoA arrays directly
- [ ] `calculateTIC()` - sum intensity SoA array
- [ ] `getBasePeak()` - find max in intensity SoA array

**2.2 Add new index-based search API**
- [ ] `Size findFirstMZ(double mz) const` - binary search on `mz_`
- [ ] `Size findLastMZ(double mz) const`
- [ ] `std::pair<Size, Size> findMZRange(double mz_min, double mz_max) const`
- [ ] `findNearest()` - update to use SoA
- [ ] `findHighestInWindow()` - use SoA array access

**2.3 Deprecate iterator-based search**
- [ ] Add `[[deprecated]]` to `MZBegin()`, `MZEnd()` with message pointing to new methods
- [ ] Iterator methods internally call new index-based methods

**2.4 Testing**
- [ ] All existing tests still pass
- [ ] Performance benchmarks show improvement
- [ ] New index-based search tests added

### Phase 3: Flip Authority (SoA becomes primary)

**Strategy:** SoA becomes the authoritative storage. AoS inheritance removed. Compatibility layer provides Peak1D access.

**3.1 Remove std::vector<Peak1D> inheritance**
- [ ] Remove `private std::vector<Peak1D>` base class
- [ ] Remove all `using ContainerType::*` declarations
- [ ] Remove `data()` returning `Peak1D*`
- [ ] Update constructors to only initialize SoA

**3.2 Add compatibility layer**
```cpp
[[deprecated("Use getMZ(i)/getIntensity(i) instead")]]
Peak1D getPeak(Size i) const {
    return Peak1D((*mz_)[i], (*intensity_)[i]);
}

[[deprecated("Use setMZ(i)/setIntensity(i) instead")]]
void setPeak(Size i, const Peak1D& peak) {
    ensureUniqueMZ();
    ensureUniqueIntensity();
    (*mz_)[i] = peak.getMZ();
    (*intensity_)[i] = peak.getIntensity();
}

[[deprecated("Use getMZArray()/getIntensityArray() instead")]]
std::vector<Peak1D> toAoS() const;

void fromAoS(const std::vector<Peak1D>& peaks);  // Not deprecated - useful for migration
```

**3.3 Testing**
- [ ] All tests updated to use new API
- [ ] Deprecated methods still work correctly
- [ ] No performance regression

### Phase 4: File I/O Updates

**Strategy:** Update file handlers to read/write directly to SoA arrays, eliminating conversion overhead.

**4.1 MzMLHandler (highest priority)**
- [ ] Update `populateSpectraWithData_()` to use bulk `setData()`
- [ ] Leverage mzML's native SoA format for direct array copy
- [ ] Update ion mobility handling to use dedicated array
- [ ] Update serialization to write from SoA arrays directly

**4.2 MzXMLHandler**
- [ ] Update peak population loop
- [ ] Use `addPeak()` instead of `push_back(Peak1D)`

**4.3 MzDataHandler**
- [ ] Similar updates

**4.4 CachedMzMLHandler, MzMLSqliteHandler**
- [ ] Update for SoA layout

**4.5 Testing**
- [ ] mzML roundtrip produces identical output
- [ ] Performance improvement measured vs baseline

### Phase 5: pyOpenMS Bindings

**Strategy:** Update Cython bindings to use new SoA API and enable zero-copy.

**5.1 Update MSSpectrum.pyx**
- [ ] Simplify `get_peaks()` to use direct array copy (not Peak1D iteration)
- [ ] Add `get_peaks_view()` for zero-copy access using `getMZArrayPtr()`
- [ ] Simplify `set_peaks()` to use bulk `setData()`
- [ ] Add `edit_peaks()` context manager for zero-copy modification
- [ ] Update ion mobility methods to use new API

**5.2 Update MSSpectrum.pxd**
- [ ] Add declarations for new methods (`getMZArrayPtr`, etc.)
- [ ] Add `shared_ptr` declarations

**5.3 Update tests**
- [ ] `test_MSSpectrum.py` - verify new API
- [ ] Add zero-copy benchmarks
- [ ] `test_tutorial.py`

### Phase 6: Library Migration

**Strategy:** Migrate OpenMS library code to use new API. Prioritize by performance impact and complexity.

**6.1 High-priority modules (performance critical)**
- [ ] FEATUREFINDER/MassTraceDetection.cpp
- [ ] ANALYSIS/ID/ (many files)
- [ ] PROCESSING/CENTROIDING/ (IM integration)

**6.2 Medium-priority modules**
- [ ] ANALYSIS/OPENSWATH/
- [ ] ANALYSIS/TARGETED/
- [ ] COMPARISON/

**6.3 Low-priority modules**
- [ ] Remaining files (mostly metadata access)

### Phase 7: TOPP Tools Migration

**7.1 Critical tools (IM or performance)**
- [ ] PeakPickerIM (IM-specialized)
- [ ] MapNormalizer (performance)
- [ ] PhosphoScoring (complex access)

**7.2 Other tools**
- [ ] Remaining 16 tools

### Phase 8: Testing & Documentation

**8.1 Comprehensive testing**
- [ ] MSSpectrum_test.cpp - full coverage of new API
- [ ] Integration tests - full pipeline tests
- [ ] Performance benchmarks - compare to baseline from Phase 0

**8.2 Documentation**
- [ ] API documentation for new methods
- [ ] Migration guide with examples for each pattern
- [ ] Performance guide for optimal SoA usage

### Phase 9: Cleanup (Future Release)

**To be done in a subsequent major release after deprecation period:**

**9.1 Remove deprecated methods**
- [ ] Remove `getPeak()`, `setPeak()`, `toAoS()`
- [ ] Remove `MZBegin()`, `MZEnd()` (iterator versions)
- [ ] Remove `containsIMData()`, `getIMData()`
- [ ] Remove compatibility typedefs

**9.2 Final cleanup**
- [ ] Remove any remaining Peak1D dependencies
- [ ] Update all documentation to remove references to deprecated API

---

## 8. Testing Strategy

### 8.1 Unit Tests

**MSSpectrum_test.cpp updates:**

```cpp
// Test new SoA API
START_SECTION(getMZArray())
{
    MSSpectrum spec;
    spec.addPeak(100.0, 1000.0f);
    spec.addPeak(200.0, 2000.0f);

    const auto& mz = spec.getMZArray();
    TEST_EQUAL(mz.size(), 2);
    TEST_REAL_SIMILAR(mz[0], 100.0);
    TEST_REAL_SIMILAR(mz[1], 200.0);
}
END_SECTION

START_SECTION(hasIonMobilityArray())
{
    MSSpectrum spec;
    TEST_EQUAL(spec.hasIonMobilityArray(), false);

    spec.setIonMobilityArray({1.0f, 2.0f}, DriftTimeUnit::MILLISECOND);
    TEST_EQUAL(spec.hasIonMobilityArray(), true);
    TEST_EQUAL(spec.getIonMobilityUnit(), DriftTimeUnit::MILLISECOND);
}
END_SECTION

START_SECTION(sortByPosition_with_IM())
{
    MSSpectrum spec;
    spec.initializeIonMobility(DriftTimeUnit::VSSC);
    spec.addPeak(300.0, 3000.0f, 3.0f);
    spec.addPeak(100.0, 1000.0f, 1.0f);
    spec.addPeak(200.0, 2000.0f, 2.0f);

    spec.sortByPosition();

    // All arrays should be sorted together
    TEST_REAL_SIMILAR(spec.getMZ(0), 100.0);
    TEST_REAL_SIMILAR(spec.getMZ(1), 200.0);
    TEST_REAL_SIMILAR(spec.getMZ(2), 300.0);
    TEST_REAL_SIMILAR(spec.getIonMobility(0), 1.0f);
    TEST_REAL_SIMILAR(spec.getIonMobility(1), 2.0f);
    TEST_REAL_SIMILAR(spec.getIonMobility(2), 3.0f);
}
END_SECTION
```

### 8.2 Performance Benchmarks

```cpp
// Benchmark: Sequential m/z access
void benchmark_mz_access() {
    MSSpectrum spec;
    // Fill with 1M peaks
    for (Size i = 0; i < 1000000; ++i) {
        spec.addPeak(i * 0.001, i * 0.1f);
    }

    // Old way (if compatibility layer exists)
    auto start = std::chrono::high_resolution_clock::now();
    double sum_old = 0;
    for (Size i = 0; i < spec.size(); ++i) {
        sum_old += spec.getPeak(i).getMZ();  // Creates temporary
    }
    auto old_time = std::chrono::high_resolution_clock::now() - start;

    // New way
    start = std::chrono::high_resolution_clock::now();
    double sum_new = 0;
    const auto& mz = spec.getMZArray();
    for (double m : mz) {
        sum_new += m;
    }
    auto new_time = std::chrono::high_resolution_clock::now() - start;

    // Expect 2-4x improvement
    std::cout << "Old: " << old_time.count() << "ns, New: " << new_time.count() << "ns\n";
}
```

### 8.3 Integration Tests

- Roundtrip: Load mzML → modify → save → reload → compare
- Ion mobility: Load IM data → sort by IM → verify alignment
- pyOpenMS: `get_peaks()` / `set_peaks()` roundtrip

---

## 9. pyOpenMS Bindings

### 9.1 Current Implementation Analysis

**Location:** `src/pyOpenMS/addons/MSSpectrum.pyx`

Current `get_peaks()` iterates through Peak1D:
```python
def get_peaks(self):
    cdef _MSSpectrum * spec_ = self.inst.get()
    cdef size_t n = spec_.size()
    cdef np.ndarray[np.float64_t, ndim=1] mzs = np.empty(n, dtype=np.float64)
    cdef np.ndarray[np.float32_t, ndim=1] intensities = np.empty(n, dtype=np.float32)
    for i in range(n):
        mzs[i] = deref(spec_)[i].getMZ()           # Peak1D access
        intensities[i] = deref(spec_)[i].getIntensity()
    return mzs, intensities
```

### 9.2 Updated Implementation

```python
def get_peaks(self):
    """Returns (mz_array, intensity_array) as numpy arrays."""
    cdef _MSSpectrum * spec_ = self.inst.get()
    cdef const vector[double]* mz_vec = &spec_.getMZArray()
    cdef const vector[float]* int_vec = &spec_.getIntensityArray()
    cdef size_t n = mz_vec.size()

    if n == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float32)

    # Direct memory copy - no Peak1D iteration
    cdef np.ndarray[np.float64_t, ndim=1] mzs = np.empty(n, dtype=np.float64)
    cdef np.ndarray[np.float32_t, ndim=1] intensities = np.empty(n, dtype=np.float32)
    memcpy(&mzs[0], mz_vec.data(), n * sizeof(double))
    memcpy(&intensities[0], int_vec.data(), n * sizeof(float))

    return mzs, intensities

def set_peaks(self, peaks):
    """Sets peaks from (mz_array, intensity_array) tuple."""
    mzs, intensities = peaks
    cdef _MSSpectrum * spec_ = self.inst.get()

    # Convert to vectors and call bulk setter
    cdef vector[double] mz_vec
    cdef vector[float] int_vec
    mz_vec.assign(&mzs[0], &mzs[0] + len(mzs))
    int_vec.assign(&intensities[0], &intensities[0] + len(intensities))

    spec_.setData(move(mz_vec), move(int_vec))
```

### 9.3 New Ion Mobility Methods

```python
def has_ion_mobility_array(self):
    """Check if spectrum has per-peak ion mobility data."""
    return self.inst.get().hasIonMobilityArray()

def get_ion_mobility_array(self):
    """Get ion mobility values as numpy array."""
    if not self.has_ion_mobility_array():
        return None
    cdef const vector[float]* im_vec = &self.inst.get().getIonMobilityArray()
    cdef size_t n = im_vec.size()
    cdef np.ndarray[np.float32_t, ndim=1] im = np.empty(n, dtype=np.float32)
    memcpy(&im[0], im_vec.data(), n * sizeof(float))
    return im

def set_ion_mobility_array(self, np.ndarray[float, ndim=1] im, int unit):
    """Set ion mobility array with unit."""
    cdef vector[float] im_vec
    im_vec.assign(&im[0], &im[0] + len(im))
    self.inst.get().setIonMobilityArray(move(im_vec), <_DriftTimeUnit>unit)
```

### 9.4 Backward Compatibility

Keep deprecated methods with warnings:
```python
def containsIMData(self):
    """Deprecated: Use has_ion_mobility_array() instead."""
    import warnings
    warnings.warn("containsIMData() is deprecated, use has_ion_mobility_array()",
                  DeprecationWarning, stacklevel=2)
    return self.has_ion_mobility_array()
```

---

## 10. File I/O Considerations

### 10.1 MzMLHandler Priority

**Current bottleneck:** `populateSpectraWithData_()` accounts for ~50% of load time.

**Current implementation (simplified):**
```cpp
// Peak-by-peak for filtered case
for (Size n = 0; n < data.size(); n += 2) {
    Peak1D p;
    p.setMZ(mz_data[n]);
    p.setIntensity(int_data[n]);
    if (passes_filter(p)) {
        spectrum.push_back(p);  // AoS insertion
    }
}
```

**SoA implementation:**
```cpp
// Bulk operation - much faster
std::vector<double> mz_filtered;
std::vector<float> int_filtered;
mz_filtered.reserve(data.size() / 2);
int_filtered.reserve(data.size() / 2);

for (Size n = 0; n < data.size(); n += 2) {
    if (passes_filter(mz_data[n], int_data[n])) {
        mz_filtered.push_back(mz_data[n]);
        int_filtered.push_back(int_data[n]);
    }
}

spectrum.setData(std::move(mz_filtered), std::move(int_filtered));
```

### 10.2 Ion Mobility in mzML

**Current:** IM stored in FloatDataArray with CV term name
**New:** During load, detect IM array and populate dedicated member

```cpp
// In MzMLHandler
if (isIonMobilityArray(float_array)) {
    DriftTimeUnit unit = getIMUnit(float_array);
    spectrum.setIonMobilityArray(std::move(float_array), unit);
} else {
    spectrum.getFloatDataArrays().push_back(std::move(float_array));
}
```

### 10.3 Serialization Changes

When writing mzML:
```cpp
// Write IM array as FloatDataArray with proper CV term
if (spectrum.hasIonMobilityArray()) {
    FloatDataArray im_fda;
    im_fda = spectrum.getIonMobilityArray();
    im_fda.setName(getIMCVTerm(spectrum.getIonMobilityUnit()));
    writeFloatDataArray(im_fda);
}
```

---

## 11. TOPP Tools Migration

### 11.1 High-Priority Tools

**PeakPickerIM** (IM-specialized):
```cpp
// Before
const auto& im_data = spectrum.getFloatDataArrays()[im_data_index];
for (Size i = 0; i < spectrum.size(); ++i) {
    double mz = spectrum[i].getMZ();
    float im = im_data[i];
}

// After
const auto& mz = spectrum.getMZArray();
const auto& im = spectrum.getIonMobilityArray();
for (Size i = 0; i < spectrum.size(); ++i) {
    process(mz[i], im[i]);
}
```

**MapNormalizer** (performance critical):
```cpp
// Before
for (Peak1D& pk : ms) {
    pk.setIntensity(pk.getIntensity() / max);
}

// After
auto& intensity = ms.getIntensityArray();
for (float& i : intensity) {
    i /= max;
}
// Or vectorized:
std::transform(intensity.begin(), intensity.end(), intensity.begin(),
               [max](float i) { return i / max; });
```

### 11.2 Migration Checklist

For each TOPP tool:
- [ ] Identify all MSSpectrum access patterns
- [ ] Convert indexed `spec[i]` to `spec.getMZ(i)` / `spec.getIntensity(i)`
- [ ] Convert range-based loops to array iteration
- [ ] Convert `push_back(Peak1D)` to `addPeak(mz, intensity)`
- [ ] Update IM access to use new methods
- [ ] Test functionality
- [ ] Verify performance

---

## 12. Risk Assessment

### 12.1 High Risk

| Risk | Mitigation |
|------|------------|
| Breaking downstream code | Compatibility layer with deprecation warnings |
| Performance regression in edge cases | Comprehensive benchmarking |
| Ion mobility data loss | Extensive IM-specific testing |
| pyOpenMS API break | Maintain method signatures, update internals |

### 12.2 Medium Risk

| Risk | Mitigation |
|------|------------|
| Incomplete migration | Track all files, phased rollout |
| Test coverage gaps | Add tests for new API before migration |
| Documentation lag | Update docs in parallel with code |

### 12.3 Low Risk

| Risk | Mitigation |
|------|------------|
| Memory overhead (temporary dual storage) | Phase 1 only, removed in Phase 3 |
| Build time increase | Minimal header changes |

---

## 13. Success Metrics

### 13.1 Performance Targets

| Metric | Target |
|--------|--------|
| Sequential m/z iteration | 2x faster |
| TIC calculation | 4x faster |
| mzML load time | 20% faster |
| Memory per spectrum | 10% reduction |

### 13.2 Code Quality Targets

| Metric | Target |
|--------|--------|
| Test coverage | ≥95% for new API |
| Deprecated method warnings | All compatibility methods marked |
| Documentation | All new methods documented |
| Migration guide | Complete with examples |

### 13.3 Compatibility Targets

| Metric | Target |
|--------|--------|
| pyOpenMS API compatibility | 100% method availability |
| TOPP tool functionality | All tools pass tests |
| File format roundtrip | Bit-identical output |

---

## 14. Zero-Copy NumPy Integration

One of the key benefits of the SoA layout is enabling **zero-copy data transfer** between C++ and Python/NumPy. This section details the implementation strategies and safety considerations.

### 14.1 Why Zero-Copy Matters

**Current AoS Limitation:**
```
Memory Layout (AoS - Peak1D):
┌──────────┬───────────┬─────────┬──────────┬───────────┬─────────┐
│ mz₁ (8B) │ int₁ (4B) │ pad(4B) │ mz₂ (8B) │ int₂ (4B) │ pad(4B) │ ...
└──────────┴───────────┴─────────┴──────────┴───────────┴─────────┘
     ↑ m/z values are NOT contiguous - NumPy cannot view this
```

NumPy requires contiguous memory for a given dtype. With AoS, extracting m/z or intensity requires copying each value individually.

**SoA Enables Direct Views:**
```
Memory Layout (SoA):
mz_ array:        ┌──────────┬──────────┬──────────┬──────────┐
                  │ mz₁ (8B) │ mz₂ (8B) │ mz₃ (8B) │ mz₄ (8B) │ ...
                  └──────────┴──────────┴──────────┴──────────┘
                       ↑ Contiguous float64 - NumPy can view directly!

intensity_ array: ┌───────────┬───────────┬───────────┬───────────┐
                  │ int₁ (4B) │ int₂ (4B) │ int₃ (4B) │ int₄ (4B) │ ...
                  └───────────┴───────────┴───────────┴───────────┘
                       ↑ Contiguous float32 - NumPy can view directly!
```

### 14.2 Performance Impact

| Operation | Copy-based | Zero-copy | Improvement |
|-----------|------------|-----------|-------------|
| `get_peaks()` for 1M peaks | ~50 ms | ~0.001 ms | 50,000x |
| `set_peaks()` for 1M peaks | ~50 ms | ~0.001 ms | 50,000x |
| Memory for 1M peak access | 12 MB allocation | 0 bytes | ∞ |

### 14.3 Implementation Approaches

#### Approach A: Typed Memory Views (Recommended for Read)

Cython can expose C++ vector data as typed memory views:

```python
# In MSSpectrum.pyx

def get_mz_array_view(self):
    """
    Get m/z values as a NumPy array VIEW (zero-copy, read-only).

    WARNING: The returned array is a view into internal C++ memory.
    Do NOT store this array beyond the lifetime of the spectrum.
    Do NOT modify the spectrum while using this view.

    Returns:
        np.ndarray: Read-only float64 array viewing internal m/z data
    """
    cdef _MSSpectrum* spec_ = self.inst.get()
    cdef const vector[double]* mz_vec = &spec_.getMZArray()
    cdef size_t n = mz_vec.size()

    if n == 0:
        return np.empty(0, dtype=np.float64)

    # Create a NumPy array that VIEWS the C++ memory (no copy!)
    cdef np.ndarray[np.float64_t, ndim=1] arr = np.asarray(
        <np.float64_t[:n]> mz_vec.data()
    )
    arr.flags.writeable = False  # Prevent accidental modification
    return arr

def get_intensity_array_view(self):
    """Get intensity values as a NumPy array VIEW (zero-copy, read-only)."""
    cdef _MSSpectrum* spec_ = self.inst.get()
    cdef const vector[float]* int_vec = &spec_.getIntensityArray()
    cdef size_t n = int_vec.size()

    if n == 0:
        return np.empty(0, dtype=np.float32)

    cdef np.ndarray[np.float32_t, ndim=1] arr = np.asarray(
        <np.float32_t[:n]> int_vec.data()
    )
    arr.flags.writeable = False
    return arr

def get_peaks_view(self):
    """
    Get (mz, intensity) as NumPy array VIEWS (zero-copy).

    Returns:
        tuple: (mz_view, intensity_view) - both read-only
    """
    return self.get_mz_array_view(), self.get_intensity_array_view()
```

#### Approach B: Writable Views with Explicit Lock (Advanced)

For in-place modification, use an explicit locking context:

```python
# In MSSpectrum.pyx

@contextmanager
def edit_peaks(self):
    """
    Context manager for zero-copy in-place peak editing.

    Usage:
        with spectrum.edit_peaks() as (mz, intensity):
            intensity *= 2.0  # Modify in-place, no copy!
            mz += 0.001       # Shift all m/z values
        # Changes are committed when context exits

    WARNING: Do not resize arrays or store references outside the context.
    """
    cdef _MSSpectrum* spec_ = self.inst.get()
    cdef vector[double]* mz_vec = &spec_.getMZArray()
    cdef vector[float]* int_vec = &spec_.getIntensityArray()
    cdef size_t n = mz_vec.size()

    if n == 0:
        yield np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float32)
        return

    # Create writable views
    cdef np.ndarray[np.float64_t, ndim=1] mz_arr = np.asarray(
        <np.float64_t[:n]> mz_vec.data()
    )
    cdef np.ndarray[np.float32_t, ndim=1] int_arr = np.asarray(
        <np.float32_t[:n]> int_vec.data()
    )

    try:
        yield mz_arr, int_arr
    finally:
        # Update ranges after modification
        spec_.updateRanges()
```

**Usage example:**
```python
# Zero-copy intensity normalization
with spectrum.edit_peaks() as (mz, intensity):
    max_int = intensity.max()
    if max_int > 0:
        intensity /= max_int  # In-place, no copy!

# Zero-copy m/z recalibration
with spectrum.edit_peaks() as (mz, intensity):
    mz += calibration_offset  # Shift all m/z values in-place
```

#### Approach C: Buffer Protocol (Most Flexible)

Expose the C++ vectors via Python's buffer protocol for maximum interoperability:

```cpp
// In C++ header (MSSpectrum.h)
class MSSpectrum {
public:
    // Expose raw pointers for Python buffer protocol
    const double* mzData() const noexcept { return mz_.data(); }
    double* mzData() noexcept { return mz_.data(); }
    const float* intensityData() const noexcept { return intensity_.data(); }
    float* intensityData() noexcept { return intensity_.data(); }
};
```

```python
# In MSSpectrum.pyx - implement buffer protocol

def __getbuffer__(self, Py_buffer *buffer, int flags):
    """Expose m/z array via buffer protocol."""
    cdef _MSSpectrum* spec_ = self.inst.get()
    cdef vector[double]* mz_vec = &spec_.getMZArray()

    buffer.buf = <void*>mz_vec.data()
    buffer.len = mz_vec.size() * sizeof(double)
    buffer.itemsize = sizeof(double)
    buffer.format = 'd'  # double
    buffer.ndim = 1
    buffer.shape = <Py_ssize_t*>&mz_vec.size()
    buffer.strides = <Py_ssize_t*>&buffer.itemsize
    buffer.suboffsets = NULL
    buffer.readonly = 0 if (flags & PyBUF_WRITABLE) else 1
```

### 14.4 Safety Considerations

#### Lifetime Management

**Problem:** NumPy view becomes a dangling pointer if:
1. The MSSpectrum object is deleted
2. The internal vector reallocates (e.g., `addPeak()` called)

**Solution 1: Reference counting with PyCapsule**
```python
def get_mz_array_view_safe(self):
    """Safe zero-copy view with prevented deallocation."""
    import weakref

    # Get the view
    view = self._get_mz_array_view_internal()

    # Prevent spectrum deallocation while view exists
    view._spectrum_ref = self  # Strong reference

    return view
```

**Solution 2: Copy-on-resize semantics**
```cpp
// In C++ - track if views are active
class MSSpectrum {
    mutable bool views_active_ = false;

public:
    void addPeak(double mz, float intensity) {
        if (views_active_) {
            throw Exception::IllegalArgument(
                "Cannot modify spectrum while Python views are active"
            );
        }
        // ... normal implementation
    }

    void lockForViews() const { views_active_ = true; }
    void unlockViews() const { views_active_ = false; }
};
```

#### Thread Safety

Zero-copy views require careful handling in multi-threaded contexts:

```python
# Thread-safe pattern
import threading

class SpectrumViewLock:
    def __init__(self, spectrum):
        self.spectrum = spectrum
        self.lock = threading.RLock()

    @contextmanager
    def read_view(self):
        with self.lock:
            self.spectrum._lock_for_views()
            try:
                yield self.spectrum.get_peaks_view()
            finally:
                self.spectrum._unlock_views()
```

### 14.5 API Design

#### Naming Convention

| Method | Behavior | Use Case |
|--------|----------|----------|
| `get_peaks()` | Returns copies (safe) | Default, always safe |
| `get_peaks_view()` | Returns views (fast) | Performance-critical, careful use |
| `edit_peaks()` | Context manager for modification | In-place editing |
| `get_mz_array()` | Copy of m/z | Backward compatible |
| `get_mz_array_view()` | View of m/z | Zero-copy read |

#### Documentation Requirements

```python
def get_peaks_view(self):
    """
    Get peak data as zero-copy NumPy array views.

    This method returns views directly into the C++ memory, avoiding
    any data copying. This is ~50,000x faster than get_peaks() for
    large spectra.

    ⚠️  SAFETY WARNINGS:

    1. LIFETIME: Views become invalid if the spectrum is deleted or
       modified (addPeak, clear, sort, etc.). Always use views within
       a limited scope.

    2. READ-ONLY: Views are marked read-only by default. Use edit_peaks()
       context manager for modification.

    3. NO STORAGE: Do not store views in data structures. Copy if needed:
       ```
       mz_copy = spectrum.get_mz_array_view().copy()  # Safe to store
       ```

    4. THREAD SAFETY: Not thread-safe. Use external locking if sharing
       spectrum between threads.

    Returns:
        tuple: (mz_view, intensity_view)
            - mz_view: np.ndarray[float64], read-only view
            - intensity_view: np.ndarray[float32], read-only view

    Example:
        >>> # Fast read-only access
        >>> mz, intensity = spectrum.get_peaks_view()
        >>> tic = intensity.sum()  # No copy made
        >>> base_peak_mz = mz[intensity.argmax()]  # No copy made

        >>> # If you need to store, copy explicitly
        >>> mz_stored = mz.copy()

    See Also:
        - get_peaks(): Safe copy-based access (always use if unsure)
        - edit_peaks(): Context manager for zero-copy modification
    """
```

### 14.6 Benchmarks

Add benchmarks to verify zero-copy performance:

```python
# In tests/benchmark_zerocopy.py

import pyopenms
import numpy as np
import timeit

def benchmark_access():
    # Create large spectrum
    n = 1_000_000
    spec = pyopenms.MSSpectrum()
    spec.set_peaks((np.random.rand(n), np.random.rand(n).astype(np.float32)))

    # Benchmark copy-based access
    copy_time = timeit.timeit(
        lambda: spec.get_peaks(),
        number=100
    ) / 100

    # Benchmark zero-copy access
    view_time = timeit.timeit(
        lambda: spec.get_peaks_view(),
        number=100
    ) / 100

    print(f"Copy: {copy_time*1000:.3f} ms")
    print(f"View: {view_time*1000:.6f} ms")
    print(f"Speedup: {copy_time/view_time:.0f}x")

def benchmark_modification():
    n = 1_000_000
    spec = pyopenms.MSSpectrum()
    spec.set_peaks((np.random.rand(n), np.random.rand(n).astype(np.float32)))

    # Copy-modify-set pattern
    def copy_modify():
        mz, intensity = spec.get_peaks()
        intensity *= 2.0
        spec.set_peaks((mz, intensity))

    # Zero-copy modify pattern
    def zerocopy_modify():
        with spec.edit_peaks() as (mz, intensity):
            intensity *= 2.0

    copy_time = timeit.timeit(copy_modify, number=10) / 10
    view_time = timeit.timeit(zerocopy_modify, number=10) / 10

    print(f"Copy-modify-set: {copy_time*1000:.1f} ms")
    print(f"Zero-copy edit:  {view_time*1000:.3f} ms")
    print(f"Speedup: {copy_time/view_time:.0f}x")
```

### 14.7 Ion Mobility Zero-Copy

Same patterns apply to ion mobility:

```python
def get_ion_mobility_array_view(self):
    """Get ion mobility values as zero-copy view."""
    if not self.has_ion_mobility_array():
        return None

    cdef _MSSpectrum* spec_ = self.inst.get()
    cdef const vector[float]* im_vec = &spec_.getIonMobilityArray()
    cdef size_t n = im_vec.size()

    cdef np.ndarray[np.float32_t, ndim=1] arr = np.asarray(
        <np.float32_t[:n]> im_vec.data()
    )
    arr.flags.writeable = False
    return arr

@contextmanager
def edit_peaks_with_im(self):
    """Zero-copy editing including ion mobility."""
    # ... similar to edit_peaks but yields (mz, intensity, im)
```

### 14.8 Integration with DataFrame Export

Optimize DataFrame creation with zero-copy:

```python
def get_df_zerocopy(self, columns=None):
    """
    Create pandas DataFrame using zero-copy where possible.

    Note: DataFrame will hold references to spectrum data.
    Spectrum must not be modified while DataFrame is in use.
    """
    import pandas as pd

    mz, intensity = self.get_peaks_view()

    data = {
        'mz': mz,           # Zero-copy view
        'intensity': intensity,  # Zero-copy view
    }

    if self.has_ion_mobility_array():
        data['ion_mobility'] = self.get_ion_mobility_array_view()

    # Scalar columns still need to be broadcast (but that's fast)
    n = len(mz)
    data['rt'] = np.full(n, self.getRT())
    data['ms_level'] = np.full(n, self.getMSLevel(), dtype=np.uint8)

    return pd.DataFrame(data)
```

---

## Appendix A: File List

### A.1 Core Files to Modify

```
src/openms/include/OpenMS/KERNEL/MSSpectrum.h
src/openms/source/KERNEL/MSSpectrum.cpp
src/openms/include/OpenMS/KERNEL/Peak1D.h
```

### A.2 I/O Files to Modify

```
src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp
src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp
src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp
src/openms/source/FORMAT/HANDLERS/CachedMzMLHandler.cpp
src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp
```

### A.3 pyOpenMS Files to Modify

```
src/pyOpenMS/pxds/MSSpectrum.pxd
src/pyOpenMS/addons/MSSpectrum.pyx
src/pyOpenMS/tests/unittests/test_MSSpectrum.py
```

### A.4 Test Files to Modify

```
src/tests/class_tests/openms/source/MSSpectrum_test.cpp
src/tests/class_tests/openms/source/Peak1D_test.cpp
src/tests/class_tests/openms/source/MSExperiment_test.cpp
src/tests/class_tests/openms/source/PeakPickerHiRes_test.cpp
src/tests/class_tests/openms/source/PeakPickerIM_test.cpp
src/tests/class_tests/openms/source/IMDataConverter_test.cpp
```

---

## Appendix B: API Quick Reference

### New Methods (Primary API)

```cpp
// ═══════════════════════════════════════════════════════════════
// ARRAY ACCESS (reference - for C++ performance code)
// ═══════════════════════════════════════════════════════════════
const std::vector<double>& getMZArray() const noexcept;
const std::vector<float>& getIntensityArray() const noexcept;
std::vector<double>& getMZArray() noexcept;  // Triggers COW if shared
std::vector<float>& getIntensityArray() noexcept;  // Triggers COW if shared

// ═══════════════════════════════════════════════════════════════
// SHARED_PTR ACCESS (for Python zero-copy bindings)
// ═══════════════════════════════════════════════════════════════
std::shared_ptr<const std::vector<double>> getMZArrayPtr() const noexcept;
std::shared_ptr<const std::vector<float>> getIntensityArrayPtr() const noexcept;
std::shared_ptr<const std::vector<float>> getIonMobilityArrayPtr() const;

// ═══════════════════════════════════════════════════════════════
// INDEXED ACCESS (convenience)
// ═══════════════════════════════════════════════════════════════
double getMZ(Size i) const;
float getIntensity(Size i) const;
void setMZ(Size i, double mz);
void setIntensity(Size i, float intensity);

// ═══════════════════════════════════════════════════════════════
// ION MOBILITY
// ═══════════════════════════════════════════════════════════════
bool hasIonMobilityArray() const noexcept;
const std::vector<float>& getIonMobilityArray() const;
std::vector<float>& getIonMobilityArray();
DriftTimeUnit getIonMobilityUnit() const noexcept;
float getIonMobility(Size i) const;
void setIonMobility(Size i, float im);

// ═══════════════════════════════════════════════════════════════
// BULK OPERATIONS
// ═══════════════════════════════════════════════════════════════
void setData(std::vector<double> mz, std::vector<float> intensity);
void setData(std::vector<double> mz, std::vector<float> intensity,
             std::vector<float> ion_mobility, DriftTimeUnit im_unit);
void setIonMobilityArray(std::vector<float> im, DriftTimeUnit unit);
void clearIonMobility();
void initializeIonMobility(DriftTimeUnit unit);

// ═══════════════════════════════════════════════════════════════
// ADDING PEAKS
// ═══════════════════════════════════════════════════════════════
void addPeak(double mz, float intensity);
void addPeak(double mz, float intensity, float ion_mobility);

// ═══════════════════════════════════════════════════════════════
// SEARCH (index-based - replaces iterator-based MZBegin/MZEnd)
// ═══════════════════════════════════════════════════════════════
Size findFirstMZ(double mz) const;
Size findLastMZ(double mz) const;
std::pair<Size, Size> findMZRange(double mz_min, double mz_max) const;
```

### Deprecated Methods (Compatibility)

```cpp
// Peak1D access (creates temporary objects)
[[deprecated("Use getMZ(i)/getIntensity(i) instead")]]
Peak1D getPeak(Size i) const;

[[deprecated("Use setMZ(i)/setIntensity(i) instead")]]
void setPeak(Size i, const Peak1D& peak);

[[deprecated("Use getMZArray()/getIntensityArray() instead")]]
std::vector<Peak1D> toAoS() const;

// Ion mobility (old API)
[[deprecated("Use hasIonMobilityArray() instead")]]
bool containsIMData() const;

[[deprecated("Use getIonMobilityArray() and getIonMobilityUnit() instead")]]
std::pair<Size, DriftTimeUnit> getIMData() const;

// Iterator-based search (replaced by index-based)
[[deprecated("Use findFirstMZ() instead")]]
Iterator MZBegin(CoordinateType mz);

[[deprecated("Use findLastMZ() instead")]]
Iterator MZEnd(CoordinateType mz);
```

### Removed Methods (Phase 3+)

```cpp
// No longer available after Phase 3
Peak1D& operator[](Size i);           // Use getMZ(i)/getIntensity(i)
const Peak1D& operator[](Size i) const;
iterator begin();                      // No Peak1D iteration
iterator end();
Peak1D* data();                        // Use getMZArray().data()
void push_back(const Peak1D&);         // Use addPeak(mz, intensity)
void emplace_back(mz, intensity);      // Use addPeak(mz, intensity)
```

---

## Appendix C: Future Extensibility - Adding New Dimensions

This section discusses architectural considerations for extending MSSpectrum with additional per-peak dimensions beyond m/z, intensity, and ion mobility.

### C.1 Potential Future Dimensions

| Dimension | Type | Use Case |
|-----------|------|----------|
| Collision Cross Section (CCS) | `float` | Ion mobility calibrated measurements |
| Peak Width / FWHM | `float` | Peak shape characterization |
| Signal-to-Noise Ratio | `float` | Quality metrics |
| Charge State | `int8_t` | Per-peak charge annotation |
| Isotope Cluster ID | `int32_t` | Grouping isotopic peaks |
| Retention Time (per-peak) | `float` | 2D chromatographic separations |
| Peak Quality Score | `float` | Machine learning quality predictions |
| Baseline | `float` | Local baseline estimates |

### C.2 Design Options for Extensibility

#### Option A: Explicit Named Members (Current Approach for IM)

Add each new dimension as an explicit member:

```cpp
class MSSpectrum {
    std::vector<double> mz_;
    std::vector<float> intensity_;
    std::vector<float> ion_mobility_;
    std::vector<float> ccs_;              // New
    std::vector<float> peak_width_;       // New
    std::vector<int8_t> charge_state_;    // New
    // ...
};
```

**Pros:**
- Type-safe access
- Clear API: `getCCSArray()`, `getChargeStateArray()`
- Compile-time errors for invalid access
- Best performance (no indirection)

**Cons:**
- Requires code changes for each new dimension
- Class grows with each addition
- All dimensions always present (even if empty)

#### Option B: Templated Dimension Registry

Use a type-safe registry pattern:

```cpp
// Dimension tag types
struct DimensionMZ {};
struct DimensionIntensity {};
struct DimensionIonMobility {};
struct DimensionCCS {};

template<typename Dim> struct DimensionTraits;
template<> struct DimensionTraits<DimensionMZ> {
    using value_type = double;
    static constexpr const char* name = "mz";
};
template<> struct DimensionTraits<DimensionIonMobility> {
    using value_type = float;
    static constexpr const char* name = "ion_mobility";
    static constexpr bool optional = true;
};

class MSSpectrum {
    // Core storage
    std::vector<double> mz_;
    std::vector<float> intensity_;

    // Optional dimensions stored in type-erased map
    std::unordered_map<std::type_index, std::any> optional_dimensions_;

public:
    template<typename Dim>
    const std::vector<typename DimensionTraits<Dim>::value_type>&
    getDimension() const;

    template<typename Dim>
    bool hasDimension() const;

    template<typename Dim>
    void setDimension(std::vector<typename DimensionTraits<Dim>::value_type> data);
};

// Usage:
spec.setDimension<DimensionCCS>(ccs_values);
if (spec.hasDimension<DimensionCCS>()) {
    const auto& ccs = spec.getDimension<DimensionCCS>();
}
```

**Pros:**
- Type-safe at compile time
- Extensible without modifying MSSpectrum class
- Only stores dimensions that are used
- Clear, expressive API

**Cons:**
- `std::any` has overhead (type erasure, heap allocation)
- More complex implementation
- Requires C++17

#### Option C: Named Dimension Map with Type Variants

Use a string-keyed map with variant types:

```cpp
class MSSpectrum {
    std::vector<double> mz_;
    std::vector<float> intensity_;
    std::vector<float> ion_mobility_;

    // Extensible typed arrays
    using DimensionVariant = std::variant<
        std::vector<float>,
        std::vector<double>,
        std::vector<int32_t>,
        std::vector<int8_t>
    >;
    std::unordered_map<std::string, DimensionVariant> extra_dimensions_;

public:
    template<typename T>
    const std::vector<T>& getExtraDimension(const std::string& name) const {
        return std::get<std::vector<T>>(extra_dimensions_.at(name));
    }

    template<typename T>
    void setExtraDimension(const std::string& name, std::vector<T> data) {
        extra_dimensions_[name] = std::move(data);
    }

    bool hasExtraDimension(const std::string& name) const {
        return extra_dimensions_.count(name) > 0;
    }
};

// Usage:
spec.setExtraDimension<float>("ccs", ccs_values);
spec.setExtraDimension<int8_t>("charge", charge_values);
const auto& ccs = spec.getExtraDimension<float>("ccs");
```

**Pros:**
- Fully extensible at runtime
- No class modification needed for new dimensions
- Supports different data types via variant

**Cons:**
- String keys are error-prone (typos)
- Type must be specified at call site
- Some runtime overhead

#### Option D: Hybrid Approach (Recommended for Future)

Keep core dimensions as explicit members, use typed registry for extensions:

```cpp
class MSSpectrum {
    // ═══════════════════════════════════════════════════════════════
    // CORE DIMENSIONS (always present conceptually, SoA layout)
    // ═══════════════════════════════════════════════════════════════
    std::vector<double> mz_;
    std::vector<float> intensity_;

    // ═══════════════════════════════════════════════════════════════
    // WELL-KNOWN OPTIONAL DIMENSIONS (first-class support)
    // ═══════════════════════════════════════════════════════════════
    std::vector<float> ion_mobility_;
    DriftTimeUnit ion_mobility_unit_ = DriftTimeUnit::NONE;

    std::vector<float> ccs_;              // Collision cross section
    std::vector<float> peak_fwhm_;        // Peak width

    // ═══════════════════════════════════════════════════════════════
    // EXTENSIBLE DIMENSIONS (for future/custom needs)
    // ═══════════════════════════════════════════════════════════════
    struct TypedArray {
        enum class Type { Float, Double, Int32, Int8, String };
        Type type;
        std::string name;
        std::variant<
            std::vector<float>,
            std::vector<double>,
            std::vector<int32_t>,
            std::vector<int8_t>,
            std::vector<std::string>
        > data;
    };
    std::vector<TypedArray> extra_dimensions_;

public:
    // Well-known dimensions have dedicated API
    bool hasCCS() const { return !ccs_.empty(); }
    const std::vector<float>& getCCSArray() const { return ccs_; }

    // Extension dimensions use generic API
    template<typename T>
    void addDimension(const std::string& name, std::vector<T> data);

    template<typename T>
    const std::vector<T>& getDimension(const std::string& name) const;
};
```

### C.3 Synchronization Requirements

Any new dimension must participate in:

1. **Sorting:** When `sortByPosition()` is called, all dimension arrays must be permuted
2. **Selection:** `select(indices)` must filter all arrays
3. **Resize/Clear:** All arrays must be resized together
4. **Validation:** Size consistency checks across all arrays

```cpp
void MSSpectrum::applyPermutation(const std::vector<Size>& indices) {
    mz_ = permute(mz_, indices);
    intensity_ = permute(intensity_, indices);

    if (hasIonMobilityArray()) {
        ion_mobility_ = permute(ion_mobility_, indices);
    }
    if (hasCCS()) {
        ccs_ = permute(ccs_, indices);
    }

    // Permute all extension dimensions
    for (auto& dim : extra_dimensions_) {
        std::visit([&indices](auto& vec) {
            vec = permute(vec, indices);
        }, dim.data);
    }

    // Permute metadata arrays (existing)
    for (auto& arr : float_data_arrays_) {
        if (arr.size() == indices.size()) {
            arr = permute(arr, indices);
        }
    }
}
```

### C.4 File Format Considerations

New dimensions need serialization support:

| Format | Approach |
|--------|----------|
| mzML | Add as new CV-termed binary data arrays |
| mzXML | Not supported (format is legacy) |
| SQLite/sqMass | Add new columns or JSON blob |
| Cached binary | Include in serialized layout |

```xml
<!-- Example mzML encoding for CCS -->
<binaryDataArray encodedLength="...">
  <cvParam cvRef="MS" accession="MS:XXXXXXX" name="collision cross section array"/>
  <binary>...</binary>
</binaryDataArray>
```

### C.5 pyOpenMS Binding Pattern

For each new dimension, add Python methods:

```python
# In MSSpectrum.pyx
def has_ccs_array(self):
    return self.inst.get().hasCCS()

def get_ccs_array(self):
    if not self.has_ccs_array():
        return None
    cdef const vector[float]* ccs = &self.inst.get().getCCSArray()
    return np.array(<float[:ccs.size()]>ccs.data(), dtype=np.float32)

def set_ccs_array(self, np.ndarray[float, ndim=1] ccs):
    cdef vector[float] ccs_vec
    ccs_vec.assign(&ccs[0], &ccs[0] + len(ccs))
    self.inst.get().setCCSArray(move(ccs_vec))
```

### C.6 Recommendations

1. **Short-term (this refactoring):** Implement IM as first-class; keep FloatDataArrays for other metadata

2. **Medium-term:** If CCS becomes commonly used, promote it to first-class member

3. **Long-term:** If >3 optional dimensions become common:
   - Consider Option D (hybrid approach)
   - Implement compile-time dimension traits for type safety
   - Add generic iteration over all dimensions for bulk operations

4. **Design principle:**
   - Core dimensions (m/z, intensity) are always explicit members
   - Frequently-used optional dimensions get first-class API
   - Rare/experimental dimensions use the extensible mechanism
   - Never sacrifice type safety for convenience

### C.7 Migration Path for New Dimensions

When adding a new first-class dimension (e.g., CCS):

1. Add member variable and unit enum if applicable
2. Add accessor methods (`hasCCS()`, `getCCSArray()`, etc.)
3. Update `applyPermutation()` and `select()`
4. Update `clear()` and `resize()`
5. Update file I/O handlers (mzML priority)
6. Update pyOpenMS bindings
7. Add tests
8. Document CV terms and units

---

## Appendix D: shared_ptr Wrapping for Data Arrays

This section analyzes the implications of wrapping the SoA data arrays (`mz_`, `intensity_`, `ion_mobility_`) in `std::shared_ptr` for enhanced Python interoperability, particularly in the context of autowrap.

### D.1 Current Patterns in OpenMS

OpenMS already uses `shared_ptr` for data arrays in the OpenSwath interfaces:

```cpp
// From OpenMS/INTERFACES/DataStructures.h
struct BinaryDataArray {
    std::vector<double> data;
};
typedef std::shared_ptr<BinaryDataArray> BinaryDataArrayPtr;

struct Spectrum {
    std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;

    BinaryDataArrayPtr getMZArray() const {
        return binaryDataArrayPtrs[0];  // Returns shared_ptr
    }
};
```

**pyOpenMS wrapper pattern (OSSpectrum.pyx):**
```python
def get_mz_array(self):
    cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
    cdef size_t n = _r.get().data.size()
    cdef double[::1] arr = <double [:n]>_r.get().data.data()
    return np.asarray(arr.copy())  # Still copies!

def get_mz_array_mv(self):
    """Memory view - no copy but lifetime risk."""
    cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
    cdef double[::1] arr = <double [:n]>_r.get().data.data()
    return arr  # Zero-copy but shared_ptr doesn't help with lifetime
```

### D.2 Proposed Design Options

#### Option 1: Direct shared_ptr for Each Array

```cpp
class MSSpectrum {
private:
    std::shared_ptr<std::vector<double>> mz_;
    std::shared_ptr<std::vector<float>> intensity_;
    std::shared_ptr<std::vector<float>> ion_mobility_;

public:
    // Returns shared_ptr - Python can hold reference
    std::shared_ptr<std::vector<double>> getMZArrayPtr() const {
        return mz_;
    }

    // Convenience: direct access (still zero-copy)
    const std::vector<double>& getMZArray() const {
        return *mz_;
    }
};
```

**Pros:**
- Python can hold `shared_ptr`, preventing deallocation while NumPy view exists
- Enables true Copy-on-Write (COW) semantics
- Autowrap handles `shared_ptr` return types well

**Cons:**
- Pointer indirection on every access: `mz_->operator[]()` vs `mz_[]`
- Memory overhead: 16 bytes (2 pointers) per shared_ptr × 3 arrays = 48 bytes/spectrum
- Atomic reference counting overhead on copy/assignment
- Complicates internal code (must check for null, dereference)

#### Option 2: Wrapper Struct with shared_ptr (Like BinaryDataArray)

```cpp
// Reuse existing pattern
struct MZArray {
    std::vector<double> data;
};
typedef std::shared_ptr<MZArray> MZArrayPtr;

struct IntensityArray {
    std::vector<float> data;
};
typedef std::shared_ptr<IntensityArray> IntensityArrayPtr;

class MSSpectrum {
private:
    MZArrayPtr mz_;
    IntensityArrayPtr intensity_;
    // ...
};
```

**Pros:**
- Consistent with existing OpenSwath interfaces
- Clear semantic: the array is a managed resource

**Cons:**
- Same overhead as Option 1
- Extra struct indirection
- More complex type system

#### Option 3: Separate Ownership Token

```cpp
class MSSpectrum {
private:
    std::vector<double> mz_;
    std::vector<float> intensity_;
    std::vector<float> ion_mobility_;

    // Shared ownership token for Python views
    mutable std::shared_ptr<MSSpectrum*> ownership_token_;

public:
    // Returns shared_ptr to token that prevents deallocation
    std::shared_ptr<MSSpectrum*> getOwnershipToken() const {
        if (!ownership_token_) {
            ownership_token_ = std::make_shared<MSSpectrum*>(
                const_cast<MSSpectrum*>(this)
            );
        }
        return ownership_token_;
    }
};
```

**Pros:**
- No indirection for internal array access
- Single shared_ptr overhead instead of per-array
- Python holds token, spectrum stays alive

**Cons:**
- Doesn't prevent internal reallocation (addPeak could invalidate views)
- More complex lifetime semantics
- Non-standard pattern

### D.3 Autowrap Integration

Autowrap can handle `shared_ptr` return types in several ways:

**1. Direct shared_ptr return (supported):**
```cpp
// C++ header
std::shared_ptr<std::vector<double>> getMZArrayPtr() const;
```

```python
# Generated by autowrap
cdef shared_ptr[libcpp_vector[double]] _r = self.inst.get().getMZArrayPtr()
```

**2. Manual addon for zero-copy (recommended):**
```python
# In MSSpectrum.pyx addon
def get_mz_array_view(self):
    """Zero-copy view with shared ownership."""
    cdef shared_ptr[libcpp_vector[double]] arr_ptr = \
        self.inst.get().getMZArrayPtr()
    cdef size_t n = arr_ptr.get().size()
    if n == 0:
        return np.empty(0, dtype=np.float64)

    # Create numpy array viewing the data
    cdef np.ndarray[np.float64_t, ndim=1] view = np.asarray(
        <np.float64_t[:n]> arr_ptr.get().data()
    )

    # CRITICAL: Store shared_ptr in array to prevent deallocation
    # NumPy's base mechanism keeps the shared_ptr alive
    cdef object capsule = PyCapsule_New(
        <void*>new shared_ptr[libcpp_vector[double]](arr_ptr),
        NULL,
        _shared_ptr_destructor
    )
    view.base = capsule  # NumPy keeps capsule alive

    view.flags.writeable = False
    return view

# Destructor for PyCapsule
cdef void _shared_ptr_destructor(object capsule):
    cdef shared_ptr[libcpp_vector[double]]* ptr = \
        <shared_ptr[libcpp_vector[double]]*>PyCapsule_GetPointer(capsule, NULL)
    del ptr
```

### D.4 Copy-on-Write (COW) Semantics

With `shared_ptr<vector<T>>`, COW can be implemented:

```cpp
class MSSpectrum {
private:
    std::shared_ptr<std::vector<double>> mz_;

    // Ensure unique ownership before modification
    void ensureUniqueMZ() {
        if (mz_.use_count() > 1) {
            mz_ = std::make_shared<std::vector<double>>(*mz_);  // Deep copy
        }
    }

public:
    void setMZ(Size i, double value) {
        ensureUniqueMZ();  // COW: copy if shared
        (*mz_)[i] = value;
    }

    // For bulk operations, caller can check shared state
    bool isMZShared() const { return mz_.use_count() > 1; }
};
```

**Usage pattern:**
```python
# Python
view = spectrum.get_mz_array_view()  # ref_count = 2 (C++ + Python)
spectrum.setMZ(0, 500.0)             # COW triggers deep copy
# view still points to original data
# spectrum.mz_ now points to new copy
```

### D.5 Performance Impact Analysis

| Operation | Direct vector | shared_ptr<vector> | Overhead |
|-----------|---------------|-------------------|----------|
| `mz_[i]` read | 1 ptr deref | 2 ptr deref | ~1 cycle |
| `mz_.size()` | 1 ptr deref | 2 ptr deref | ~1 cycle |
| Copy spectrum | Deep copy all | Shallow copy (COW) | **Big win** |
| addPeak | Direct append | Deref + append | ~1 cycle |
| Memory/spectrum | 0 bytes | 48 bytes (3×16) | +48 bytes |
| Python view creation | N/A | O(1) | **Safe** |

**Conclusion:** The overhead is minimal for typical use. The main benefit is safe Python zero-copy views with proper lifetime management.

### D.6 Recommendation

**Recommended Approach: Option 1 with Lazy Initialization**

```cpp
class MSSpectrum {
private:
    // Use shared_ptr for Python interop, lazy-initialized
    mutable std::shared_ptr<std::vector<double>> mz_shared_;
    mutable std::shared_ptr<std::vector<float>> intensity_shared_;
    mutable std::shared_ptr<std::vector<float>> ion_mobility_shared_;

    // Internal storage (primary)
    std::vector<double> mz_;
    std::vector<float> intensity_;
    std::vector<float> ion_mobility_;

public:
    // Primary API: direct access (no overhead)
    const std::vector<double>& getMZArray() const noexcept {
        return mz_;
    }

    // Python API: returns shared_ptr for safe views
    std::shared_ptr<const std::vector<double>> getMZArrayShared() const {
        if (!mz_shared_ || mz_shared_.get() != &mz_) {
            // Create alias shared_ptr that doesn't own but prevents this from dying
            mz_shared_ = std::shared_ptr<const std::vector<double>>(
                shared_from_this(), &mz_
            );
        }
        return mz_shared_;
    }
};
```

**Alternative: Simple shared_ptr throughout (simpler, slight overhead)**

For a simpler implementation with acceptable overhead:

```cpp
class MSSpectrum {
private:
    std::shared_ptr<std::vector<double>> mz_ =
        std::make_shared<std::vector<double>>();
    std::shared_ptr<std::vector<float>> intensity_ =
        std::make_shared<std::vector<float>>();
    std::shared_ptr<std::vector<float>> ion_mobility_;  // nullptr if not used

public:
    // Fast access (1 extra deref, negligible)
    const std::vector<double>& getMZArray() const {
        return *mz_;
    }

    // For Python zero-copy
    std::shared_ptr<std::vector<double>> getMZArrayPtr() const {
        return mz_;
    }
};
```

### D.7 Implementation Considerations

**1. Null Safety:**
```cpp
// Always initialize core arrays
MSSpectrum::MSSpectrum()
    : mz_(std::make_shared<std::vector<double>>())
    , intensity_(std::make_shared<std::vector<float>>())
    , ion_mobility_(nullptr)  // Optional, null when not used
{}
```

**2. Thread Safety:**
- `shared_ptr` reference counting is atomic (thread-safe)
- Array modification is NOT thread-safe (same as current)
- For parallel algorithms, copy spectrum or use external locking

**3. Serialization:**
- No impact on mzML/mzXML I/O (serialize vector contents, not shared_ptr)

**4. Move Semantics:**
```cpp
MSSpectrum::MSSpectrum(MSSpectrum&& other) noexcept
    : mz_(std::move(other.mz_))
    , intensity_(std::move(other.intensity_))
    , ion_mobility_(std::move(other.ion_mobility_))
{
    // other's shared_ptrs are now null
}
```

### D.8 Summary Decision Matrix

| Factor | Direct vector | shared_ptr<vector> | Recommendation |
|--------|---------------|-------------------|----------------|
| Python zero-copy safety | Poor (dangling) | **Excellent** | shared_ptr |
| Internal access speed | **Optimal** | 1 cycle slower | Acceptable |
| Memory overhead | **None** | 48 bytes/spec | Acceptable |
| COW semantics | Manual | **Built-in** | shared_ptr |
| Implementation complexity | **Simple** | Moderate | Trade-off |
| Autowrap compatibility | N/A | **Native** | shared_ptr |

**Final Recommendation:** Use `shared_ptr<std::vector<T>>` for the SoA arrays. The benefits for Python interoperability and COW semantics outweigh the minimal performance overhead. The pattern is already established in OpenSwath interfaces.

---

## Document History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12-10 | Initial plan document |
| 1.1 | 2025-12-10 | Added future extensibility (Appendix C) |
| 1.2 | 2025-12-10 | Added zero-copy NumPy integration (Section 14) |
| 1.3 | 2025-12-10 | Added shared_ptr analysis (Appendix D) |
| 1.4 | 2025-12-10 | Comprehensive revision: reconciled shared_ptr decision, added mzML alignment section, updated file counts, added MZBegin/MZEnd migration, revised implementation phases, added version timeline |

---

*Document Version: 1.4*
*Date: 2025-12-10*
*Author: Claude (AI Assistant)*
