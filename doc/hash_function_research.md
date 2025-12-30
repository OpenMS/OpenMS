# OpenMS Hash Function Research Report

## Executive Summary

This document presents a comprehensive analysis of hash function implementations in OpenMS and identifies classes that would benefit from hash functions. The research found **9 existing hash implementations** across 11 files, and identified **37+ core classes** that have `operator==` defined but lack hash functions.

---

## 1. Existing Hash Function Implementations

### 1.1 std::hash Template Specializations

| Class | File | Implementation |
|-------|------|----------------|
| `OpenMS::String` | `DATASTRUCTURES/String.h:491-499` | Delegates to `std::hash<std::string>` |
| `OpenMS::Index` | `ANALYSIS/ID/AhoCorasickAmbiguous.h:257-263` | Hashes underlying `uint32_t` |

### 1.2 Boost-Compatible hash_value Functions

| Type | File | Implementation |
|------|------|----------------|
| `OpenMS::String` | `DATASTRUCTURES/String.h:488` + `String.cpp` | Uses `boost::hash<std::string>` |
| `DPosition<N,T>` | `ML/CLUSTERING/HashGrid.h:501-509` | XOR of coordinate hashes |

### 1.3 Custom Hash Functor Classes

| Functor | File | Purpose |
|---------|------|---------|
| `ProteinHitAccessionHash` | `METADATA/ProteinHit.h:49-56` | Hash ProteinHit by accession |
| `ProteinHitPtrAccessionHash` | `METADATA/ProteinHit.h:58-65` | Hash ProteinHit* by accession |
| `FNV1aHasher` | `FORMAT/ControlledVocabulary.h:34-46` | Platform-independent FNV-1a for CV terms |
| `MyUIntSetHasher` | `ANALYSIS/ID/IDBoostGraph.cpp:40-47` | Hash set of vertex indices |
| `accessionHash_` | `ANALYSIS/ID/IDMergerAlgorithm.h:242-244` | Static function for ProteinHit hashing |

---

## 2. Classes Recommended for Hash Functions

### 2.1 HIGH PRIORITY (Performance-Critical)

These classes are used as keys in `std::map`/`std::set` in hot paths and would benefit significantly from `std::unordered_*` containers.

#### AASequence
- **File:** `CHEMISTRY/AASequence.h`
- **Has operator==:** Yes (line 547)
- **Current Usage:**
  - `std::set<AASequence>` in `GridFeature.h`, `QTCluster.h`
  - `std::map<AASequence, ...>` in `PeptideAndProteinQuant.h`, `ConsensusIDAlgorithm.h`, `FFIDAlgoExternalIDHandler.h`
  - Similarity cache: `std::map<std::pair<AASequence, AASequence>, double>`
- **Fields to Hash:**
  - `peptide_` (vector of Residue pointers)
  - `n_term_mod_`, `c_term_mod_` (terminal modifications)
- **Impact:** HIGH - Used in feature clustering, quantitation, consensus ID

#### EmpiricalFormula
- **File:** `CHEMISTRY/EmpiricalFormula.h`
- **Has operator==:** Yes (line 262)
- **Current Usage:**
  - `std::set<EmpiricalFormula>` in `TheoreticalSpectrumGenerator.cpp`
  - `std::map<EmpiricalFormula, String>` formula cache in spectrum generation
- **Fields to Hash:**
  - `formula_` (map of Element* to SignedSize)
  - `charge_` (Int)
- **Impact:** HIGH - Used in spectrum generation hot loop (see TODO at line 447)

#### PeptideHit
- **File:** `METADATA/PeptideHit.h`
- **Has operator==:** Yes
- **Current Usage:** Frequently compared for deduplication
- **Fields to Hash:**
  - `sequence_` (AASequence)
  - `charge_`
  - `score_`
- **Impact:** HIGH - Core identification data type

#### Peak1D / Peak2D / ChromatogramPeak
- **Files:** `KERNEL/Peak1D.h`, `KERNEL/Peak2D.h`, `KERNEL/ChromatogramPeak.h`
- **Has operator==:** Yes (default)
- **Fields to Hash:**
  - `position_` (DPosition)
  - `intensity_`
- **Impact:** MEDIUM-HIGH - Fundamental data types

### 2.2 MEDIUM PRIORITY (Frequently Used)

| Class | File | Fields to Hash |
|-------|------|----------------|
| `Element` | `CHEMISTRY/Element.h` | `atomic_number_`, `symbol_` |
| `Residue` | `CHEMISTRY/Residue.h` | `one_letter_code_`, `modification_` |
| `PeptideEvidence` | `METADATA/PeptideEvidence.h` | `accession_`, `start_`, `end_` |
| `PeptideIdentification` | `METADATA/PeptideIdentification.h` | `RT_`, `MZ_`, identifier |
| `Feature` / `BaseFeature` | `KERNEL/Feature.h` | `unique_id_`, position, charge |
| `FeatureHandle` | `KERNEL/FeatureHandle.h` | `map_index_`, `unique_id_` |
| `NASequence` | `CHEMISTRY/NASequence.h` | Nucleotide sequence vector |
| `ResidueModification` | `CHEMISTRY/ResidueModification.h` | `id_`, `full_id_` |
| `CVTerm` | `METADATA/CVTerm.h` | `accession_`, `cv_ref_` |

### 2.3 LOWER PRIORITY (Metadata/Complex)

| Class | File | Reason |
|-------|------|--------|
| `MSSpectrum` | `KERNEL/MSSpectrum.h` | Large container, rarely used as key |
| `MSChromatogram` | `KERNEL/MSChromatogram.h` | Large container |
| `ProteinIdentification` | `METADATA/ProteinIdentification.h` | Complex container type |
| `ExperimentalSettings` | `METADATA/ExperimentalSettings.h` | Metadata only |
| `DigestionEnzyme` | `CHEMISTRY/DigestionEnzyme.h` | Typically in lookup tables |

---

## 3. Identified Issues

### 3.1 Unused/Problematic Type Definition
**File:** `ANALYSIS/MAPMATCHING/QTClusterFinder.h:84-86`
```cpp
typedef std::unordered_map<
  std::pair<OpenMS::GridFeature*, OpenMS::GridFeature*>,
  double> PairDistances;
```
**Issue:** `std::pair` has no default hash in C++. This typedef appears unused but would fail to compile if instantiated.

### 3.2 Redundant Data Structures
**File:** `TheoreticalSpectrumGenerator.cpp:447`
```cpp
// TODO why do you need a separate set for the losses? Just use the keys from the formula_str_cache?
```
**Issue:** Maintaining both `std::set<EmpiricalFormula>` and `std::map<EmpiricalFormula, String>` is redundant. Adding hash to EmpiricalFormula would enable `std::unordered_map` with better performance.

---

## 4. Recommended Implementation Pattern

Based on existing OpenMS patterns, new hash functions should follow this template:

### 4.1 For Classes in std Namespace (Preferred)
```cpp
// In header file, after class definition
namespace std {
  template<> struct hash<OpenMS::ClassName> {
    std::size_t operator()(const OpenMS::ClassName& obj) const noexcept {
      std::size_t seed = 0;
      // Use boost::hash_combine or similar
      boost::hash_combine(seed, obj.getField1());
      boost::hash_combine(seed, obj.getField2());
      return seed;
    }
  };
}
```

### 4.2 For Boost Compatibility
```cpp
// In OpenMS namespace
namespace OpenMS {
  inline std::size_t hash_value(const ClassName& obj) {
    std::size_t seed = 0;
    boost::hash_combine(seed, obj.getField1());
    boost::hash_combine(seed, obj.getField2());
    return seed;
  }
}
```

### 4.3 Hash Combine Implementation
If avoiding boost dependency:
```cpp
inline void hash_combine(std::size_t& seed, std::size_t value) {
  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
```

---

## 5. Implementation Priority Roadmap

### Phase 1: Critical (Immediate Performance Benefit)
1. **AASequence** - Enables unordered containers in feature clustering, quantitation
2. **EmpiricalFormula** - Fixes spectrum generation performance issue

### Phase 2: High Value
3. **PeptideHit** - Core ID data type
4. **Peak1D**, **Peak2D**, **ChromatogramPeak** - Fundamental types

### Phase 3: Medium Value
5. **Element**, **Residue** - Chemistry building blocks
6. **PeptideEvidence** - Protein mapping
7. **FeatureHandle** - Feature referencing

### Phase 4: Completeness
8. Remaining METADATA classes
9. Remaining CHEMISTRY classes

---

## 6. Files to Modify

| Priority | Class | Header File |
|----------|-------|-------------|
| 1 | AASequence | `src/openms/include/OpenMS/CHEMISTRY/AASequence.h` |
| 1 | EmpiricalFormula | `src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h` |
| 2 | PeptideHit | `src/openms/include/OpenMS/METADATA/PeptideHit.h` |
| 2 | Peak1D | `src/openms/include/OpenMS/KERNEL/Peak1D.h` |
| 2 | Peak2D | `src/openms/include/OpenMS/KERNEL/Peak2D.h` |
| 2 | ChromatogramPeak | `src/openms/include/OpenMS/KERNEL/ChromatogramPeak.h` |
| 3 | Element | `src/openms/include/OpenMS/CHEMISTRY/Element.h` |
| 3 | Residue | `src/openms/include/OpenMS/CHEMISTRY/Residue.h` |
| 3 | PeptideEvidence | `src/openms/include/OpenMS/METADATA/PeptideEvidence.h` |
| 3 | FeatureHandle | `src/openms/include/OpenMS/KERNEL/FeatureHandle.h` |

---

## 7. Testing Recommendations

For each new hash function:
1. Verify hash consistency: `a == b` implies `hash(a) == hash(b)`
2. Test distribution: Verify reasonable spread across hash values
3. Test with unordered_set/map: Ensure proper behavior
4. Performance benchmark: Compare map vs unordered_map for typical use cases

---

## Appendix A: Current Hash Function Locations

```
src/openms/include/OpenMS/
├── ANALYSIS/ID/
│   ├── AhoCorasickAmbiguous.h      # std::hash<Index>
│   ├── ConsensusMapMergerAlgorithm.h # accessionHash_
│   └── IDMergerAlgorithm.h         # accessionHash_
├── DATASTRUCTURES/
│   └── String.h                    # std::hash<String>, hash_value
├── FORMAT/
│   └── ControlledVocabulary.h      # FNV1aHasher
├── METADATA/
│   └── ProteinHit.h               # ProteinHitAccessionHash
└── ML/CLUSTERING/
    └── HashGrid.h                  # hash_value(DPosition)

src/openms/source/
├── ANALYSIS/ID/
│   └── IDBoostGraph.cpp            # MyUIntSetHasher
└── DATASTRUCTURES/
    └── String.cpp                  # hash_value implementation
```

---

*Generated: 2025-12-30*
*Research conducted on OpenMS commit: e3c7431*
