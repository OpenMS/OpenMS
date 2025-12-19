# Plan: Switching to LightTransition to Reduce In-Memory Size of TargetedExperiments

**Issue:** [#8532](https://github.com/OpenMS/OpenMS/issues/8532)
**Problem:** Loading a 10GB TSV file with ~30 million transitions requires 15-20GB+ peak memory
**Root Cause:** `ReactionMonitoringTransition` inherits from `CVTermList`, costing ~1000 bytes/transition vs ~136 bytes for `LightTransition`

---

## Executive Summary

The codebase already has a two-tier architecture with lightweight classes (`LightTransition`, `LightTargetedExperiment`), but they are underutilized. The main strategy is to:

1. **Extend lightweight loading paths** to more tools
2. **Implement lazy CV term loading** for cases requiring full metadata
3. **Optimize the conversion pipeline** to avoid intermediate heavy objects
4. **Add direct-to-light loading** bypassing `TargetedExperiment` entirely where possible

---

## Phase 1: Analysis and Baseline Measurements

### Task 1.1: Create Memory Benchmarking Infrastructure
- **File:** Create new test file `src/tests/class_tests/openms/source/TargetedExperimentMemoryTest.cpp`
- **Purpose:** Establish baseline memory measurements for various transition counts
- **Metrics to capture:**
  - Peak memory for loading 1K, 10K, 100K, 1M, 10M transitions
  - Memory per transition for both `ReactionMonitoringTransition` and `LightTransition`
  - Time to load and convert between formats

### Task 1.2: Audit Current Lightweight Path Usage
- **Files to audit:**
  - `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp` (lines 574-718)
  - `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.cpp` (lines 103-196)
  - All TOPP tools in `src/topp/` that use `TargetedExperiment`
- **Output:** Document which tools use heavy vs light paths

---

## Phase 2: Optimize TransitionTSVFile Loading

### Task 2.1: Add Direct LightTargetedExperiment Loading API
**Current state:** `TransitionTSVFile` has two `convertTSVToTargetedExperiment()` overloads - one for `TargetedExperiment` (heavy) and one for `LightTargetedExperiment` (light)

**Changes required:**
- **File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h`
  - Add new method signature: `void convertTSVToLightTargetedExperiment(const char* filename, FileTypes::Type filetype, LightTargetedExperiment& exp) const`

- **File:** `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`
  - Optimize `TSVToTargetedExperiment_()` for `LightTargetedExperiment` (lines 649-718) to:
    - Avoid creating intermediate `TSVTransition` objects where possible
    - Use move semantics for string assignments
    - Pre-allocate vectors based on estimated row count

### Task 2.2: Optimize TSVTransition Parsing
- **File:** `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`
- **Function:** `readUnstructuredTSVInput_()` (lines 192-466)
- **Changes:**
  - Add optional streaming mode to avoid loading entire file into memory
  - Implement chunked processing for very large files
  - Add parameter to skip parsing columns not needed for light loading

### Task 2.3: Add Memory-Efficient Mode to TransitionTSVFile
- **File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h`
- **Changes:**
  - Add parameter `memory_efficient_mode` (boolean)
  - When enabled, skip parsing of optional columns that only feed CV terms
  - Document which columns are skipped in this mode

---

## Phase 3: Extend LightTargetedExperiment Capabilities

### Task 3.1: Add Missing Fields to LightTransition
**File:** `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` (lines 19-117)

**Current `LightTransition` fields:**
```cpp
std::string transition_name;
std::string peptide_ref;
double library_intensity{};
double product_mz{};
double precursor_mz{};
double precursor_im{-1};
int fragment_charge{};
bool decoy{};
bool detecting_transition{};
bool quantifying_transition{};
bool identifying_transition{};
```

**Fields to add for broader tool compatibility:**
```cpp
// Retention time (critical for many workflows)
double retention_time{-1};
double retention_time_lower_bound{-1};
double retention_time_upper_bound{-1};

// Fragment annotation (for assay generation compatibility)
std::string fragment_type;  // "y", "b", etc.
int fragment_series_number{-1};

// Optional: product charge (separate from fragment_charge for clarity)
int product_charge{0};
```

### Task 3.2: Add Conversion Methods
**File:** `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.cpp`

**New methods:**
```cpp
// Convert single transition (for on-demand heavy-to-light)
static LightTransition convertTransition(const ReactionMonitoringTransition& tr);

// Convert single transition back (for tools that need to write)
static ReactionMonitoringTransition convertToFullTransition(const LightTransition& lt);

// Batch conversion with memory optimization
static void convertTransitionsInPlace(
    const std::vector<ReactionMonitoringTransition>& heavy,
    std::vector<LightTransition>& light,
    bool clear_heavy_after = true);
```

### Task 3.3: Extend LightTargetedExperiment Container
**File:** `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` (lines 184-256)

**Add methods:**
```cpp
// Get total memory footprint estimate
size_t estimateMemoryUsage() const;

// Reserve capacity for known sizes
void reserve(size_t num_transitions, size_t num_compounds, size_t num_proteins);

// Clear and release memory
void clear();
void shrink_to_fit();
```

---

## Phase 4: Update TOPP Tools to Use Lightweight Path

### Task 4.1: High-Priority Tools (OpenSWATH Workflow)
These tools process large transition libraries and would benefit most:

| Tool | File | Current Behavior | Required Change |
|------|------|------------------|-----------------|
| `OpenSwathWorkflow` | `src/topp/OpenSwathWorkflow.cpp` | Uses light path internally | Verify optimization is active |
| `OpenSwathRTNormalizer` | `src/topp/OpenSwathRTNormalizer.cpp` | Check path | Switch to light if not already |
| `OpenSwathDIAPreScoring` | `src/topp/OpenSwathDIAPreScoring.cpp` | Check path | Switch to light if not already |
| `OpenSwathChromatogramExtractor` | Part of workflow | Check path | Ensure light path used |

### Task 4.2: Medium-Priority Tools (Can Be Optimized)
These tools read transitions but don't necessarily need full CV terms:

| Tool | File | Current Behavior | Required Change |
|------|------|------------------|-----------------|
| `MRMMapper` | `src/topp/MRMMapper.cpp` | Uses TargetedExperiment | Add light loading option |
| `MRMTransitionGroupPicker` | `src/topp/MRMTransitionGroupPicker.cpp` | Uses full transitions | Evaluate if light is sufficient |
| `FeatureFinderIdentification` | `src/topp/FeatureFinderIdentification.cpp` | Uses full transitions | Add light loading option |

### Task 4.3: Tools Requiring Full Transitions (No Change)
These tools write transitions or manipulate CV terms - must keep using `TargetedExperiment`:

| Tool | Reason |
|------|--------|
| `TargetedFileConverter` | Writes TraML/TSV/PQP with full metadata |
| `OpenSwathAssayGenerator` | Adds CV terms to transitions |
| `OpenSwathDecoyGenerator` | Modifies and annotates transitions |
| `AssayGeneratorMetabo` | Creates new transitions with annotations |

---

## Phase 5: Implement Lazy CV Term Loading

### Task 5.1: Create CVTermProxy Class
**New file:** `src/openms/include/OpenMS/METADATA/CVTermProxy.h`

```cpp
/// Lightweight reference to CV terms stored elsewhere
class CVTermProxy
{
public:
    CVTermProxy() = default;

    /// Check if CV terms are loaded
    bool isLoaded() const;

    /// Load CV terms on demand
    void load();

    /// Get reference to actual CVTermList (loads if needed)
    const CVTermList& get() const;

private:
    mutable std::unique_ptr<CVTermList> cv_terms_;
    mutable bool loaded_{false};
    // Reference to source for lazy loading
    std::function<CVTermList()> loader_;
};
```

### Task 5.2: Create HybridTransition Class (Optional Advanced Feature)
**New file:** `src/openms/include/OpenMS/ANALYSIS/MRM/HybridTransition.h`

A middle-ground class that:
- Has the small footprint of `LightTransition`
- Can lazy-load CV terms when needed
- Is compatible with existing `ReactionMonitoringTransition` APIs

```cpp
class HybridTransition
{
public:
    // Essential fields (always in memory)
    String transition_name_;
    String peptide_ref_;
    double precursor_mz_{};
    double product_mz_{};
    double library_intensity_{};
    bool decoy_{};

    // Lazy-loaded fields
    mutable CVTermProxy cv_terms_;
    mutable std::unique_ptr<RetentionTime> retention_time_;

    // Flag tracking
    bool detecting_{true};
    bool quantifying_{true};
    bool identifying_{false};
};
```

---

## Phase 6: File Format Optimizations

### Task 6.1: Optimize TraML Loading
**File:** `src/openms/source/FORMAT/TraMLFile.cpp`

- Add option to load transitions without CV terms (metadata-only mode)
- Implement streaming parser for very large TraML files
- Add direct-to-light loading method

### Task 6.2: Optimize PQP Loading
**File:** `src/openms/source/ANALYSIS/OPENSWATH/TransitionPQPFile.cpp`

- Add method to load directly into `LightTargetedExperiment`
- Implement query-based loading (only load transitions for specific peptides)
- Add pagination support for very large databases

### Task 6.3: Add Binary Transition Format (Optional Future Enhancement)
Consider creating a compact binary format optimized for:
- Fast loading (memory-mapped)
- Small file size
- Direct deserialization to `LightTransition` vectors

---

## Phase 7: Algorithm Updates

### Task 7.1: Update ChromatogramExtractor
**File:** `src/openms/source/ANALYSIS/OPENSWATH/ChromatogramExtractor.cpp`

- Verify it works with `LightTargetedExperiment`
- Add overloads accepting `LightTransition` directly
- Remove any unnecessary conversions

### Task 7.2: Update SwathMapMassCorrection
**Files:**
- `src/openms/source/ANALYSIS/OPENSWATH/SwathMapMassCorrection.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/SwathMapMassCorrection_PASEF.cpp`

- Ensure both work with light transitions
- Add direct `LightTransition` overloads if missing

### Task 7.3: Update DIAScoring
**File:** `src/openms/source/ANALYSIS/OPENSWATH/DIAScoring.cpp`

- Verify compatibility with light transitions
- Optimize memory access patterns

---

## Phase 8: Testing and Validation

### Task 8.1: Unit Tests
**New test files:**
- `src/tests/class_tests/openms/source/LightTransitionTest.cpp`
- `src/tests/class_tests/openms/source/TransitionConversionTest.cpp`
- `src/tests/class_tests/openms/source/MemoryEfficientLoadingTest.cpp`

**Test coverage:**
- Verify `LightTransition` ↔ `ReactionMonitoringTransition` conversion roundtrip
- Verify all fields are preserved during conversion
- Memory usage assertions

### Task 8.2: Integration Tests
- Create test TSV files with 1M+ transitions
- Verify peak memory stays within acceptable bounds
- Verify output correctness matches heavy-path output

### Task 8.3: Performance Benchmarks
- Measure loading time for various file sizes
- Measure memory usage per transition
- Compare against baseline from Phase 1

---

## Phase 9: Documentation and Migration Guide

### Task 9.1: Update API Documentation
- Document new light loading methods
- Add memory usage guidance to class docs
- Update tool parameter documentation

### Task 9.2: Create Migration Guide
- Document which tools support light loading
- Provide examples for common workflows
- List breaking changes (if any)

### Task 9.3: Update User Documentation
- Add section on memory optimization for large libraries
- Document recommended workflows for 10M+ transition libraries

---

## Implementation Priority

### Immediate Impact (Implement First)
1. **Phase 2.1**: Direct LightTargetedExperiment loading API
2. **Phase 4.1**: Update OpenSWATH workflow tools
3. **Phase 8.1**: Basic unit tests

### Medium-Term Improvements
4. **Phase 3.1-3.3**: Extend LightTransition capabilities
5. **Phase 2.2-2.3**: Memory-efficient TSV parsing
6. **Phase 7.1-7.3**: Algorithm updates

### Long-Term Optimizations
7. **Phase 5**: Lazy CV term loading
8. **Phase 6**: File format optimizations
9. **Phase 9**: Documentation

---

## Expected Outcomes

| Metric | Current | After Phase 2-4 | After All Phases |
|--------|---------|-----------------|------------------|
| Memory per transition | ~1000 bytes | ~136 bytes | ~100 bytes |
| Peak memory for 30M transitions | 15-20 GB | 4-5 GB | 3-4 GB |
| Loading time | Baseline | -20% | -40% |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Breaking existing workflows | Medium | High | Extensive testing, backward-compatible APIs |
| Missing CV terms in light path | Low | Medium | Document which tools need full transitions |
| Performance regression | Low | Medium | Benchmark every change |
| Increased code complexity | Medium | Low | Good documentation, clear separation of concerns |

---

## Files to Modify (Summary)

### Core Files
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h`
- `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.cpp`

### TOPP Tools
- `src/topp/OpenSwathWorkflow.cpp`
- `src/topp/OpenSwathRTNormalizer.cpp`
- `src/topp/OpenSwathDIAPreScoring.cpp`
- `src/topp/MRMMapper.cpp`

### Test Files (New)
- `src/tests/class_tests/openms/source/LightTransitionTest.cpp`
- `src/tests/class_tests/openms/source/TransitionConversionTest.cpp`
