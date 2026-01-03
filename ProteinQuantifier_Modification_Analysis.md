# ProteinQuantifier Modification Handling - Complete Analysis

## Executive Summary

This document analyzes how ProteinQuantifier currently handles peptide modifications and proposes extensions to support optional filtering, with **special attention to the distinction between fixed and variable modifications**.

**Current State:** All modification variants of a peptide are aggregated at the protein level using the unmodified sequence.

**Proposed Extension:** Add filtering capabilities that respect the biological distinction between fixed modifications (expected, should be kept) and variable modifications (may need filtering).

---

## 1. Current Behavior

### 1.1 Peptide Level Quantification
- Each modification variant stored separately
- Key: Full `AASequence` (includes all modifications)
- Example: `PEPTIDE`, `PEPTIDE(Oxidation)`, `PEPTIDE(Phospho)` are 3 separate entries

### 1.2 Protein Level Quantification
- All variants aggregated using unmodified sequence
- Location: `PeptideAndProteinQuant.cpp`, line 988
```cpp
const String peptide = pep_q.first.toUnmodifiedString();  // Modification info discarded
prot_quant_[leader_accession].peptide_abundances[peptide][...] += abundance;
```

### 1.3 Example
```
Input peptides:
  PEPTIDE           → 1000 intensity
  PEPTIDE(Oxidation) → 500 intensity
  PEPTIDE(Phospho)   → 300 intensity

Protein output:
  PEPTIDE → 1800 total (all three aggregated)
```

This is **intentional behavior** as documented in ProteinQuantifier.cpp line 71.

---

## 2. Fixed vs Variable Modifications

### 2.1 Biological Distinction

**Fixed Modifications:**
- Always expected on specific residues
- Examples: TMT/iTRAQ labels, Carbamidomethyl (C)
- Applied during sample preparation
- Should be present on all relevant residues

**Variable Modifications:**
- May or may not occur
- Examples: Oxidation (M), Phosphorylation, Acetylation
- Represent biological variation or artifacts
- Only some peptides will have these

### 2.2 Storage in OpenMS

```cpp
// ProteinIdentification::SearchParameters (line 260-261)
struct SearchParameters {
    std::vector<String> fixed_modifications;    // e.g., ["TMT6plex", "Carbamidomethyl"]
    std::vector<String> variable_modifications; // e.g., ["Oxidation", "Phospho"]
};
```

This information is typically available in:
- ProteinIdentification objects from search results
- FeatureXML/ConsensusXML/idXML files
- Can be passed via `protein_groups` parameter

### 2.3 Why This Matters for Filtering

**Problem with generic filtering:**
A simple "filter all modified peptides" approach would incorrectly exclude fixed modifications like TMT labels, breaking labeled quantification workflows.

**Example - TMT workflow:**
```
Fixed: TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)
Variable: Oxidation (M)

Input:
  PEPTIDE(TMT)                    → 1000 (expected)
  PEPTIDE(TMT)M(Oxidation)        → 500 (artifact)

Generic "filter modified" = true:
  ❌ WRONG: Excludes both → No quantification!

Filter "variable modifications" = true:
  ✅ CORRECT: Keeps TMT (fixed), excludes oxidation (variable)
  Result: 1000 intensity
```

---

## 3. Recommended Implementation

### 3.1 Phase 1: Filter Variable Modifications (4-6 hours)

**New parameter:**
```ini
filter_variable_modifications = "false"  # Exclude variable mods from protein quant
```

**Implementation approach:**

1. **Extract search parameters** (in `quantifyProteins()`)
```cpp
std::set<String> variable_mod_names_;
std::set<String> fixed_mod_names_;

if (proteins.hasSearchParameters()) {
    const auto& params = proteins.getSearchParameters();
    fixed_mod_names_.insert(params.fixed_modifications.begin(),
                           params.fixed_modifications.end());
    variable_mod_names_.insert(params.variable_modifications.begin(),
                              params.variable_modifications.end());
}
```

2. **Helper method to check for variable modifications**
```cpp
bool PeptideAndProteinQuant::hasVariableModification_(const AASequence& seq) {
    // Check residue modifications
    for (Size i = 0; i < seq.size(); ++i) {
        if (seq[i].isModified()) {
            String mod_name = seq[i].getModificationName();
            if (variable_mod_names_.find(mod_name) != variable_mod_names_.end()) {
                return true;
            }
        }
    }
    
    // Check N-terminal modification
    if (seq.hasNTerminalModification()) {
        if (variable_mod_names_.find(seq.getNTerminalModificationName()) 
            != variable_mod_names_.end()) {
            return true;
        }
    }
    
    // Check C-terminal modification
    if (seq.hasCTerminalModification()) {
        if (variable_mod_names_.find(seq.getCTerminalModificationName()) 
            != variable_mod_names_.end()) {
            return true;
        }
    }
    
    return false;
}
```

3. **Apply filter in transferPeptideDataToProteins_()** (line ~977)
```cpp
bool filter_variable = param_.getValue("filter_variable_modifications") == "true";

for (auto const& pep_q : pep_quant_) {
    if (filter_variable && hasVariableModification_(pep_q.first)) {
        OPENMS_LOG_DEBUG << "Filtering peptide with variable modification: " 
                        << pep_q.first.toString() << endl;
        continue;
    }
    
    // ... existing aggregation code
}
```

4. **Update documentation** (ProteinQuantifier.cpp line 71)
```cpp
// NEW TEXT:
// Peptides with the same sequence, but with different modifications are 
// quantified separately on the peptide level. By default, they are treated 
// as one peptide for the protein quantification (aggregated by unmodified sequence).
// The @p filter_variable_modifications parameter allows excluding variable 
// modifications (e.g., oxidation, phosphorylation) while keeping fixed 
// modifications (e.g., TMT labels, carbamidomethylation) in the protein 
// quantification.
```

### 3.2 Phase 2: Advanced Filtering (4-8 hours)

**Extended parameters:**
```ini
modification_filter_mode = "none"  
# Options:
#   none              - No filtering (default)
#   exclude_variable  - Exclude variable modifications
#   exclude_fixed     - Exclude fixed modifications (rare use case)
#   exclude_specific  - Exclude specific modifications by name
#   include_specific  - Include only specific modifications
#   unmodified_only   - Only unmodified peptides

excluded_modifications = []  # For exclude_specific mode
included_modifications = []  # For include_specific mode
```

---

## 4. Use Cases

### 4.1 TMT/iTRAQ Labeled Quantification
```bash
# Keep TMT labels (fixed), exclude oxidation artifacts (variable)
ProteinQuantifier -in tmt_data.consensusXML -out proteins.csv \
                  -filter_variable_modifications true
```

### 4.2 Label-Free Quantification
```bash
# Keep carbamidomethyl (fixed), exclude oxidation (variable)
ProteinQuantifier -in lfq_data.featureXML -out proteins.csv \
                  -filter_variable_modifications true
```

### 4.3 Phosphoproteomics Analysis
```bash
# Separate analysis for non-phosphorylated proteins
ProteinQuantifier -in phospho_data.featureXML -out non_phospho.csv \
                  -modification_filter_mode exclude_specific \
                  -excluded_modifications "Phospho (S)" "Phospho (T)" "Phospho (Y)"

# Separate analysis for phosphorylated proteins
ProteinQuantifier -in phospho_data.featureXML -out phospho.csv \
                  -modification_filter_mode include_specific \
                  -included_modifications "Phospho (S)" "Phospho (T)" "Phospho (Y)"
```

---

## 5. Edge Cases

### 5.1 Missing SearchParameters
```cpp
if (variable_mod_names_.empty() && filter_variable) {
    OPENMS_LOG_WARN << "Warning: filter_variable_modifications is enabled but "
                   << "no search parameters found. Cannot determine which "
                   << "modifications are variable. Filtering disabled." << endl;
    filter_variable = false;
}
```

### 5.2 Modification Not in SearchParameters
```cpp
// Treat unknown modifications as variable (safer to filter)
if (mod_name not in fixed_mod_names_ && mod_name not in variable_mod_names_) {
    OPENMS_LOG_DEBUG << "Unknown modification: " << mod_name << endl;
    // Treat as variable
}
```

### 5.3 Multiple Identification Runs
```cpp
// Use union of all variable modifications from all runs
for (const auto& protein_id : protein_ids) {
    if (protein_id.hasSearchParameters()) {
        const auto& params = protein_id.getSearchParameters();
        variable_mod_names_.insert(params.variable_modifications.begin(),
                                  params.variable_modifications.end());
    }
}
```

---

## 6. Testing Strategy

### 6.1 Test Case 1: TMT with Oxidation
```
Setup:
- Fixed: TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)
- Variable: Oxidation (M)

Input:
- PEPTIDE(TMT)M(TMT)              : 1000
- PEPTIDE(TMT)M(Oxidation)(TMT)   : 500

filter_variable_modifications = true

Expected:
- Protein abundance: 1000 (oxidized excluded)
- Peptide output: Both peptides listed
```

### 6.2 Test Case 2: Phosphoproteomics
```
Setup:
- Fixed: Carbamidomethyl (C)
- Variable: Phospho (S/T/Y), Oxidation (M)

Input:
- PEPTIDE(Carbamido)              : 1000
- PEPTIDE(Carbamido)S(Phospho)    : 800
- PEPTIDE(Carbamido)M(Ox)         : 200

filter_variable_modifications = true

Expected:
- Protein abundance: 1000 (only peptide with fixed mod)
```

### 6.3 Test Case 3: No SearchParameters
```
Input: FeatureXML without SearchParameters

filter_variable_modifications = true

Expected:
- Warning logged
- Filtering disabled
- All peptides used (backward compatible)
```

---

## 7. Code Locations

### 7.1 Files to Modify

1. **ProteinQuantifier.cpp** (line ~350)
   - Add parameter: `filter_variable_modifications`
   - Update documentation (line 71)

2. **PeptideAndProteinQuant.h**
   - Add member: `std::set<String> variable_mod_names_`
   - Add member: `std::set<String> fixed_mod_names_`
   - Declare: `bool hasVariableModification_(const AASequence&)`
   - Update Statistics struct for tracking

3. **PeptideAndProteinQuant.cpp**
   - Extract search parameters in `quantifyProteins()` (~line 175)
   - Implement `hasVariableModification_()` helper
   - Apply filter in `transferPeptideDataToProteins_()` (~line 977)
   - Update statistics output

### 7.2 Test Files
- Create TOPP test: `/src/tests/topp/ProteinQuantifier_<N>_*.{input,output}`
- Add test with TMT-labeled data
- Add test with phosphoproteomics data
- Add test for missing SearchParameters edge case

---

## 8. Comparison with Original Proposal

| Aspect | Original Generic | Revised Fixed/Variable |
|--------|-----------------|------------------------|
| Parameter | `filter_modified_peptides` | `filter_variable_modifications` |
| TMT workflows | ❌ Breaks (filters labels) | ✅ Works (keeps labels) |
| Oxidation filtering | ✅ Works | ✅ Works |
| Biological meaning | Unclear | Clear and standard |
| Search alignment | None | Matches search parameters |
| User understanding | Requires explanation | Familiar from search tools |
| Effort | 2-4h | 4-6h |

---

## 9. Advantages of Fixed/Variable Approach

1. **Biologically meaningful**: Aligns with proteomics experiment design
2. **Correct for labeled quant**: TMT/iTRAQ labels preserved
3. **Flexible for PTMs**: Can separate biological PTMs from artifacts
4. **User-friendly**: Concept already familiar from search tools
5. **Backward compatible**: Default keeps all modifications
6. **Search-aware**: Uses existing SearchParameters metadata
7. **Robust**: Handles mixed workflows correctly

---

## 10. Conclusion

**The distinction between fixed and variable modifications is fundamental to proteomics** and must be respected in any modification filtering feature for ProteinQuantifier.

**Recommended immediate action:** Implement Phase 1 (filter_variable_modifications parameter)

**Estimated effort:** 
- Phase 1: 4-6 hours (basic filtering)
- Phase 2: 4-8 hours (advanced options)
- Total with tests and documentation: 10-15 hours

**Benefits:**
- Solves the user pain point (filtering artifacts like oxidation)
- Maintains correct behavior for labeled quantification
- Provides foundation for advanced PTM analysis workflows
- Aligns with standard proteomics concepts and terminology

---

## References

### Key Code Locations
- ProteinQuantifier TOPP tool: `/home/runner/work/OpenMS/OpenMS/src/topp/ProteinQuantifier.cpp`
- Core quantification logic: `/home/runner/work/OpenMS/OpenMS/src/openms/source/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.cpp`
- Header definitions: `/home/runner/work/OpenMS/OpenMS/src/openms/include/OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h`
- Modification definitions: `/home/runner/work/OpenMS/OpenMS/src/openms/include/OpenMS/CHEMISTRY/ModificationDefinition.h`
- Search parameters: `/home/runner/work/OpenMS/OpenMS/src/openms/include/OpenMS/METADATA/ProteinIdentification.h` (line 260-261)

### Critical Lines
- Line 988 of `PeptideAndProteinQuant.cpp`: Where modifications are discarded (`toUnmodifiedString()`)
- Line 71 of `ProteinQuantifier.cpp`: Documentation of current behavior
- Line 260-261 of `ProteinIdentification.h`: SearchParameters with fixed/variable mod lists
