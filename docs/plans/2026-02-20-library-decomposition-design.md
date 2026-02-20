# OpenMS Library Decomposition Design

**Date:** 2026-02-20
**Status:** Approved
**Goal:** Split the monolithic OpenMS library into 3 sub-libraries for faster incremental builds and clean dependency boundaries.

## Motivation

- **Primary:** Faster incremental builds — developers only recompile the sub-library they changed
- **Secondary:** Enforce clean dependency boundaries — prevent architectural erosion

## Current State

OpenMS is a monolithic library: ~689 .cpp files compiled into a single `libOpenMS.so`. Source is logically organized into modules (KERNEL, FORMAT, ANALYSIS, etc.) with `sources.cmake` per directory, but there are no enforced API boundaries. OpenSwathAlgo is the only separate library.

## Architecture

### Three Libraries with Strict Layering

```
┌─────────────────────────────────────────────────────┐
│                    OpenMS_Algo                       │
│  ANALYSIS, COMPARISON, FEATUREFINDER, ML,            │
│  PROCESSING, QC, APPLICATIONS, INTERFACES_IMPL       │
│  + specialized format handlers (TraML, XQuest, etc.) │
│  + extern: evergreen, IsoSpec, eol-bspline           │
│  (~340 .cpp files)                                   │
├────────────────┬────────────────────────────────────┤
│                │           OpenMS_IO                 │
│                │  FORMAT (generic handlers only)      │
│                │  + FASTAContainer (from DATASTR)     │
│                │  + OnDiscMSExperiment impl           │
│                │  + extern: SQLiteCpp, nlohmann_json, │
│                │    tool_description_lib              │
│                │  (~150 .cpp files)                   │
├────────────────┴────────────────────────────────────┤
│                    OpenMS_Core                       │
│  KERNEL, CONCEPT, DATASTRUCTURES, MATH, SYSTEM,      │
│  METADATA, CHEMISTRY, IONMOBILITY                    │
│  + extern: GTE, Quadtree, SIMDe                      │
│  (~240 .cpp files)                                   │
└─────────────────────────────────────────────────────┘
```

**Dependency rule:** `Algo → IO → Core` (strict downward only, no cycles)

**Compatibility target:** `OpenMS` INTERFACE library that transitively links all three. Existing `find_package(OpenMS)` and `target_link_libraries(... OpenMS)` keep working.

### Module Assignment

| Library | Modules |
|---------|---------|
| **OpenMS_Core** | KERNEL, CONCEPT, DATASTRUCTURES, MATH, SYSTEM, METADATA, CHEMISTRY, IONMOBILITY |
| **OpenMS_IO** | FORMAT (generic handlers: MzML, MzXML, FASTA, DTA, MGF, mzIdentML, pepXML, idXML, featureXML, consensusXML, etc.) |
| **OpenMS_Algo** | ANALYSIS, COMPARISON, FEATUREFINDER, ML, PROCESSING, QC, APPLICATIONS, INTERFACES_IMPL |

### Rationale for Core Scope

METADATA and CHEMISTRY are included in Core because:
- KERNEL headers include METADATA headers 28 times across 11 files (FeatureMap, ConsensusMap, MSSpectrum, BaseFeature, etc.)
- DATASTRUCTURES classes (Adduct, Compomer, QTCluster) depend on CHEMISTRY types (EmpiricalFormula, AASequence)
- Moving METADATA out would require refactoring 28+ header includes with potential ABI breakage
- These modules represent foundational MS data types, not algorithms

### OpenSwathAlgo

Remains a separate standalone library (zero OpenMS dependencies). Linked by OpenMS_Algo.

## Boundary Fixes Required

### Files that need to move

| File | From | To | Reason |
|------|------|----|--------|
| `FASTAContainer.h/.cpp` | DATASTRUCTURES | FORMAT (IO) | Depends on FASTAFile, is format-specific |
| `PosteriorErrorProbabilityModel` | MATH | ANALYSIS (Algo) | Depends on METADATA peptide ID types |
| `QTCluster` | DATASTRUCTURES | FEATUREFINDER or ANALYSIS (Algo) | Depends on AASequence |

### Files that need refactoring

| File | Change | Reason |
|------|--------|--------|
| `OnDiscMSExperiment` | Pimpl/bridge pattern: abstract base in Core, concrete impl in IO | Currently in KERNEL but depends on FORMAT handlers |

### Specialized format handlers moving from FORMAT to Algo

These are domain-specific serialization handlers, not generic format support:
- `TraMLHandler` / `TraMLFile` (TARGETED experiment format)
- `MRMFeatureQCFile` / `MRMFeaturePickerFile` (OpenSWATH-specific)
- `XQuestResultXMLFile` (cross-linking MS)
- `FLASHDeconvFeatureFile` / `FLASHDeconvSpectrumFile` (top-down proteomics)
- `SwathFileConsumer` (SWATH-MS)
- `AbsoluteQuantitationMethodFile` (quantitation)
- `MRMFile` (MRM-specific)
- `IBSpectraFile` (isobaric labeling)
- `GNPSMGFFile` (depends on spectrum merging)

## CMake Build System

### Directory structure (source files stay in place)

```
src/openms/
├── CMakeLists.txt           # Orchestrator
├── core/
│   └── CMakeLists.txt       # openms_add_library(OpenMS_Core ...)
├── io/
│   └── CMakeLists.txt       # openms_add_library(OpenMS_IO ...)
├── algo/
│   └── CMakeLists.txt       # openms_add_library(OpenMS_Algo ...)
├── source/                  # Source files stay where they are
│   ├── KERNEL/
│   ├── FORMAT/
│   ├── ANALYSIS/
│   └── ...
└── include/OpenMS/          # Headers stay where they are
```

### Dependency declarations

```cmake
# Core
target_link_libraries(OpenMS_Core PUBLIC
    XercesC::XercesC Qt6::Core Boost::regex
    LibSVM::LibSVM OpenMP::OpenMP_CXX)
target_link_libraries(OpenMS_Core PRIVATE
    Eigen3::Eigen ZLIB::ZLIB BZip2::BZip2
    GTE Quadtree SIMDe)

# IO
target_link_libraries(OpenMS_IO PUBLIC OpenMS_Core)
target_link_libraries(OpenMS_IO PRIVATE
    Boost::iostreams Boost::date_time
    SQLiteCpp nlohmann_json tdl::tdl)

# Algo
target_link_libraries(OpenMS_Algo PUBLIC OpenMS_IO OpenMS_Core)
target_link_libraries(OpenMS_Algo PRIVATE
    OpenSwathAlgo Evergreen IsoSpec eol-bspline)

# Compatibility umbrella
add_library(OpenMS INTERFACE)
target_link_libraries(OpenMS INTERFACE OpenMS_Algo OpenMS_IO OpenMS_Core)
```

### Key points
- Headers stay in `include/OpenMS/` — no include path changes, no `#include` rewrites
- Source files stay in `source/` — only CMake aggregation changes
- Each `sources.cmake` is included by exactly one sub-library's CMakeLists.txt
- `OpenMSConfig.cmake` exports all three targets plus the umbrella

## Test System

### Test organization — split to match libraries

```
src/tests/class_tests/
├── openms_core/            # Links OpenMS_Core only
│   ├── CMakeLists.txt
│   └── source/
├── openms_io/              # Links OpenMS_IO (+ transitively Core)
│   ├── CMakeLists.txt
│   └── source/
├── openms_algo/            # Links OpenMS_Algo (+ transitively IO, Core)
│   ├── CMakeLists.txt
│   └── source/
└── openms/                 # Legacy dir during transition
```

### Key principles
- Each test links **only** against the sub-library it tests (enforces boundaries at compile time)
- Test data files remain shared
- CTest labels: `ctest -L core`, `ctest -L io`, `ctest -L algo`
- TOPP tests stay unchanged (they test executables)

## pyOpenMS Bindings

**No changes.** pyOpenMS remains a single unified Python package:
- `.pxd` files don't change
- Autowrap generation doesn't change
- `pyopenms` links the `OpenMS` umbrella target
- Users see no difference in Python

## Migration Strategy

### Phase 1 — Establish Core library

1. Create `src/openms/core/CMakeLists.txt` building `OpenMS_Core` from KERNEL, CONCEPT, DATASTRUCTURES, MATH, SYSTEM, METADATA, CHEMISTRY, IONMOBILITY sources
2. Create `OpenMS` umbrella target linking `OpenMS_Core` + temporary `OpenMS_Rest`
3. Move a handful of Core tests to `openms_core/` test dir, verify they compile linking only Core
4. Validate: full test suite passes

### Phase 2 — Extract IO library

1. Move `FASTAContainer` from DATASTRUCTURES to FORMAT
2. Refactor `OnDiscMSExperiment` (abstract base in Core, impl in IO)
3. Create `OpenMS_IO` from generic FORMAT sources
4. Move specialized format handlers to Algo tier
5. Move `PosteriorErrorProbabilityModel` from MATH to ANALYSIS
6. Move `QTCluster` from DATASTRUCTURES to FEATUREFINDER or ANALYSIS
7. Update umbrella: `OpenMS` links `Core + IO + Rest`
8. Move IO tests, validate

### Phase 3 — Remaining becomes Algo

1. Rename `OpenMS_Rest` → `OpenMS_Algo`
2. Verify full dependency chain: `Algo → IO → Core`
3. Move remaining tests
4. Update TOPP tool linking (use `OpenMS` umbrella)
5. Update `OpenMSConfig.cmake` exports

### Phase 4 — Cleanup and validation

1. Full CTest suite passes
2. pyOpenMS builds and Python tests pass
3. Remove old monolithic build path
4. Update documentation

**Key principle:** At every phase, the project builds and all tests pass. No big-bang migration.

## Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| Hidden circular dependencies discovered during split | Phase 1 catches these early; umbrella target allows fallback |
| External consumers break | Umbrella `OpenMS` target provides full backwards compatibility |
| Build system complexity increases | Each sub-CMakeLists.txt is simple; complexity is in the orchestration |
| CI pipeline needs updating | Add per-library build/test stages incrementally |
| Symbol visibility / DLL export issues on Windows | Each library gets its own export macro (OPENMS_CORE_DLLAPI, etc.) |
