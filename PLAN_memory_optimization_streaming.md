# Memory Optimization Plan: Stream TSV Directly to LightTargetedExperiment

## Issue Reference
GitHub Issue: #8544

## Problem Statement

The current TSV transition file loading workflow stores data twice in memory:

1. **First storage**: `std::vector<TSVTransition>` - intermediate representation with **all** fields duplicated per transition row
2. **Second storage**: `LightTargetedExperiment` - efficient hierarchical structure with de-duplicated compound/protein data

### Current Memory Flow

```
TSV File (on disk)
     │
     ▼
readUnstructuredTSVInput_()
     │
     ▼
std::vector<TSVTransition>  ◄── MEMORY BOTTLENECK
     │                          - Duplicates compound data per transition
     │                          - ~30+ string fields per transition
     │                          - Full data stays in memory during conversion
     ▼
TSVToTargetedExperiment_()
     │
     ▼
LightTargetedExperiment     ◄── Efficient structure
     │                          - Compounds stored once per group_id
     │                          - Transitions reference compounds by ID
     ▼
(TSVTransition vector freed only after full conversion)
```

### Memory Duplication Analysis

For a typical large transition list with 100,000 transitions across 1,000 peptides:

| Data Field | TSVTransition (per row) | LightTargetedExperiment |
|------------|------------------------|-------------------------|
| PeptideSequence | 100,000 copies | 1,000 compounds |
| FullPeptideName | 100,000 copies | 1,000 compounds |
| ProteinName vector | 100,000 copies | 1,000 protein refs |
| GeneName | 100,000 copies | 1,000 compounds |
| CompoundName | 100,000 copies | 1,000 compounds |
| SMILES | 100,000 copies | 1,000 compounds |
| SumFormula | 100,000 copies | 1,000 compounds |
| peptide_group_label | 100,000 copies | 1,000 compounds |

**Estimated overhead**: For a 100k transition file, the TSVTransition vector can consume 500MB-1GB of memory unnecessarily.

## Proposed Solution

Stream data directly from TSV parsing into `LightTargetedExperiment`, eliminating the intermediate `TSVTransition` vector entirely.

### New Memory Flow

```
TSV File (on disk)
     │
     ▼
streamTSVToLightTargetedExperiment_()  ◄── NEW FUNCTION
     │
     │  For each line:
     │    1. Parse line into temporary single TSVTransition
     │    2. Create/lookup compound (deduplicated)
     │    3. Create LightTransition with reference
     │    4. Add to LightTargetedExperiment
     │    5. TSVTransition goes out of scope (freed)
     │
     ▼
LightTargetedExperiment (only copy of data)
```

## Implementation Plan

### Step 1: Create New Streaming Method (Header)

**File**: `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h`

Add new private method declaration:

```cpp
/** @brief Stream TSV input directly to LightTargetedExperiment (memory-efficient)
 *
 * @param[in] filename The input file
 * @param[in] filetype The type of file ("mrm" or "tsv")
 * @param[out] exp The output targeted experiment
 *
 * This method avoids creating an intermediate TSVTransition vector,
 * instead streaming each row directly into the target structure.
 */
void streamTSVToLightTargetedExperiment_(
    const char* filename,
    FileTypes::Type filetype,
    OpenSwath::LightTargetedExperiment& exp);
```

### Step 2: Refactor Line Parsing

**File**: `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`

Extract the line parsing logic into a reusable function:

```cpp
/** @brief Parse a single TSV line into a TSVTransition
 *
 * @param[in] tmp_line The parsed line fields
 * @param[in] header_dict Map of header names to column indices
 * @param[in] filetype The type of file
 * @param[out] transition The output transition
 * @param[in,out] spectrast_legacy Flag for SpectraST legacy mode
 * @return true if transition should be processed, false if it should be skipped
 */
bool parseTSVLine_(
    const std::vector<std::string>& tmp_line,
    const std::map<std::string, int>& header_dict,
    FileTypes::Type filetype,
    TSVTransition& transition,
    bool& spectrast_legacy);
```

This refactoring extracts lines 261-467 from `readUnstructuredTSVInput_()` into a standalone function.

### Step 3: Implement Streaming Conversion

**File**: `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`

Implement the new streaming method:

```cpp
void TransitionTSVFile::streamTSVToLightTargetedExperiment_(
    const char* filename,
    FileTypes::Type filetype,
    OpenSwath::LightTargetedExperiment& exp)
{
    std::ifstream data(filename);
    std::string line;
    std::map<std::string, int> header_dict;
    char delimiter = ',';

    // Setup header (same as existing code)
    if (filetype == FileTypes::MRM) {
        // ... MRM header setup
    } else {
        TextFile::getLine(data, line);
        getTSVHeader_(line, delimiter, header_dict);
    }

    // Tracking maps for deduplication (similar to TSVToTargetedExperiment_)
    std::map<String, int> compound_map;
    std::map<String, int> protein_map;

    // For mixed sequence group resolution - need two-pass or deferred resolution
    // Option A: Collect group_ids in first pass, validate, then process
    // Option B: Post-process compounds to fix labels (preferred for memory)

    bool spectrast_legacy = false;
    int line_count = 0;

    while (TextFile::getLine(data, line))
    {
        line_count++;
        std::vector<std::string> tmp_line;
        // ... parse line into tmp_line

        TSVTransition mytransition;  // Stack-allocated, freed each iteration

        if (!parseTSVLine_(tmp_line, header_dict, filetype, mytransition, spectrast_legacy))
        {
            continue;  // Skip this transition
        }

        cleanupTransitions_(mytransition);

        // Create LightTransition (same as TSVToTargetedExperiment_)
        OpenSwath::LightTransition transition;
        transition.transition_name = std::move(mytransition.transition_name);
        transition.peptide_ref = mytransition.group_id;  // Keep copy for reference
        transition.library_intensity = mytransition.library_intensity;
        transition.precursor_mz = mytransition.precursor;
        transition.product_mz = mytransition.product;
        transition.precursor_im = mytransition.drift_time;
        transition.fragment_charge = 0;
        if (!mytransition.fragment_charge.empty() && mytransition.fragment_charge != "NA")
        {
            transition.fragment_charge = mytransition.fragment_charge.toInt();
        }
        transition.decoy = mytransition.decoy;
        transition.detecting_transition = mytransition.detecting_transition;
        transition.identifying_transition = mytransition.identifying_transition;
        transition.quantifying_transition = mytransition.quantifying_transition;

        exp.transitions.push_back(std::move(transition));

        // Create compound if new (same deduplication logic)
        if (compound_map.find(mytransition.group_id) == compound_map.end())
        {
            OpenSwath::LightCompound compound;
            if (mytransition.isPeptide())
            {
                OpenMS::TargetedExperiment::Peptide tramlpeptide;
                createPeptide_(/* iterator adapter */, tramlpeptide);
                OpenSwathDataAccessHelper::convertTargetedCompound(tramlpeptide, compound);
            }
            else
            {
                OpenMS::TargetedExperiment::Compound tramlcompound;
                createCompound_(/* iterator adapter */, tramlcompound);
                OpenSwathDataAccessHelper::convertTargetedCompound(tramlcompound, compound);
            }
            exp.compounds.push_back(std::move(compound));
            compound_map[mytransition.group_id] = 0;
        }

        // Create proteins if new
        for (Size i = 0; i < mytransition.ProteinName.size(); ++i)
        {
            if (mytransition.isPeptide() &&
                protein_map.find(mytransition.ProteinName[i]) == protein_map.end())
            {
                OpenSwath::LightProtein protein;
                protein.id = std::move(mytransition.ProteinName[i]);
                protein.sequence = "";
                exp.proteins.push_back(std::move(protein));
                protein_map[mytransition.ProteinName[i]] = 0;
            }
        }

        // mytransition goes out of scope here - memory freed
    }

    // Post-process: resolve mixed sequence groups
    resolveMixedSequenceGroupsInPlace_(exp);
}
```

### Step 4: Handle Mixed Sequence Group Resolution

The current `resolveMixedSequenceGroups_()` modifies the TSVTransition vector. For streaming, we need a variant that works on `LightTargetedExperiment`:

```cpp
/** @brief Resolve mixed sequence groups in-place on LightTargetedExperiment
 *
 * Since we no longer have access to the original peptide_group_label values
 * after streaming, we need to either:
 * A) Track label->sequences during streaming and fix at the end
 * B) Store peptide_group_label in LightCompound (requires struct change)
 */
void resolveMixedSequenceGroupsInPlace_(OpenSwath::LightTargetedExperiment& exp);
```

**Note**: The `LightCompound` struct already has a `peptide_group_label` field (line 141 in TransitionExperiment.h), so this can be implemented by:
1. Building a map of `peptide_group_label -> vector<compound_ids>`
2. Checking if sequences differ within each group
3. Updating the labels accordingly

### Step 5: Adapt createPeptide_ and createCompound_ Functions

Current signatures take iterators to TSVTransition vector:
```cpp
void createPeptide_(std::vector<TSVTransition>::const_iterator tr_it, ...);
void createCompound_(std::vector<TSVTransition>::const_iterator tr_it, ...);
```

Add overloads that take a reference to a single TSVTransition:
```cpp
void createPeptide_(const TSVTransition& transition, OpenMS::TargetedExperiment::Peptide& peptide);
void createCompound_(const TSVTransition& transition, OpenMS::TargetedExperiment::Compound& compound);
```

### Step 6: Update Public API

**File**: `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp`

Modify `convertTSVToTargetedExperiment` for `LightTargetedExperiment` to use streaming:

```cpp
void TransitionTSVFile::convertTSVToTargetedExperiment(
    const char* filename,
    FileTypes::Type filetype,
    OpenSwath::LightTargetedExperiment& targeted_exp)
{
    // Use memory-efficient streaming approach
    streamTSVToLightTargetedExperiment_(filename, filetype, targeted_exp);
}
```

The `TargetedExperiment` overload can remain unchanged (it has additional metadata requirements that may need the full TSVTransition data).

### Step 7: Add Progress Reporting

Since we're streaming, we need to either:
1. Pre-count lines (fast file scan) for progress percentage
2. Use line count progress without knowing total
3. Report bytes read vs file size

Recommended: Pre-count lines with a simple loop before main processing.

## Files to Modify

| File | Changes |
|------|---------|
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h` | Add `streamTSVToLightTargetedExperiment_()`, `parseTSVLine_()` declarations |
| `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp` | Implement streaming method, refactor line parsing, add function overloads |

## Testing Strategy

1. **Unit tests**: Compare output of streaming vs. non-streaming for identical TSV inputs
2. **Memory profiling**: Use Valgrind/Massif to measure peak memory usage before/after
3. **Large file test**: Test with 100k+ transition files to verify memory reduction
4. **Regression tests**: Ensure existing OpenSWATH workflows produce identical results

## Expected Results

- **Memory reduction**: ~50% for large transition files (elimination of TSVTransition vector)
- **Performance**: Potentially faster due to reduced memory allocation/deallocation overhead
- **Functionality**: Identical output to current implementation

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Mixed sequence group resolution differs | Comprehensive unit tests comparing old vs new output |
| SpectraST file handling edge cases | Preserve exact parsing logic in `parseTSVLine_()` |
| Progress reporting less accurate | Pre-count lines or use byte-based progress |
| API compatibility | Keep existing signatures, only change internal implementation |

## Future Optimizations (Out of Scope)

1. **String interning**: Use a string pool for repeated values (protein names, gene names)
2. **Memory-mapped I/O**: For very large files, use mmap for faster reading
3. **Parallel parsing**: Parse file in chunks with multiple threads (complex due to line dependencies)
4. **Reserve vector capacity**: Pre-allocate vectors based on estimated counts from file size

## Implementation Order

1. Extract `parseTSVLine_()` function (refactoring only, no behavior change)
2. Add `createPeptide_`/`createCompound_` overloads for single TSVTransition reference
3. Implement `streamTSVToLightTargetedExperiment_()`
4. Implement `resolveMixedSequenceGroupsInPlace_()`
5. Update `convertTSVToTargetedExperiment()` to use streaming
6. Add unit tests
7. Perform memory profiling validation
