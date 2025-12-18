# Plan: Streaming-Based TSV Loading for TargetedExperiment

## Problem Statement

Loading a 10GB TSV file (~30M transitions) into `TargetedExperiment` requires 15-20GB+ peak memory due to:
1. `ReactionMonitoringTransition` inherits `CVTermList` (~1000 bytes/transition vs 136 bytes for `LightTransition`)
2. All transitions loaded simultaneously before processing
3. Multiple processing passes each create new transition vectors

## Goal

Reduce peak memory from ~15-20GB to ~500MB-1GB by processing transitions in streaming fashion, grouped by peptide.

## Architecture Overview

```
┌─────────────┐    ┌──────────────────┐    ┌─────────────────┐    ┌────────────┐
│  TSV File   │───▶│ Streaming Reader │───▶│ Peptide Grouper │───▶│  Processor │
└─────────────┘    └──────────────────┘    └─────────────────┘    └─────┬──────┘
                                                                        │
                   ┌──────────────────┐    ┌─────────────────┐          │
                   │   Output File    │◀───│ Streaming Writer│◀─────────┘
                   └──────────────────┘    └─────────────────┘

                   Persistent State (small):
                   ┌─────────────────────────────────────────┐
                   │ • peptide_map (deduplication)           │
                   │ • protein_map (deduplication)           │
                   │ • compound_map (deduplication)          │
                   │ • Peptide/Protein/Compound vectors      │
                   └─────────────────────────────────────────┘
```

## Key Insight

Transitions are processed **per peptide group** in `MRMAssay::detectingTransitions()`. We can:
1. Stream-read transitions
2. Buffer only transitions for the **current peptide group**
3. Process the complete group
4. Write results and discard buffer
5. Move to next group

## Implementation Steps

### Step 1: Create StreamingTSVReader Class

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVReader.h`

```cpp
class OPENMS_DLLAPI StreamingTransitionTSVReader
{
public:
  /// Callback signature for processing a complete peptide group
  using PeptideGroupCallback = std::function<void(
    const String& group_id,
    std::vector<TSVTransition>& transitions,
    const Peptide& peptide,
    const std::vector<Protein>& proteins
  )>;

  /// Stream-read TSV file, calling callback for each complete peptide group
  void readGrouped(
    const String& filename,
    PeptideGroupCallback callback,
    bool require_sorted = false  // If false, will sort internally
  );

  /// Get deduplicated peptides/proteins/compounds after streaming
  const std::vector<Peptide>& getPeptides() const;
  const std::vector<Protein>& getProteins() const;
  const std::vector<Compound>& getCompounds() const;

private:
  std::map<String, Size> peptide_map_;    // group_id -> index
  std::map<String, Size> protein_map_;    // protein_name -> index
  std::map<String, Size> compound_map_;   // compound_id -> index

  std::vector<Peptide> peptides_;
  std::vector<Protein> proteins_;
  std::vector<Compound> compounds_;

  void flushGroup_(const String& group_id,
                   std::vector<TSVTransition>& buffer,
                   PeptideGroupCallback& callback);
};
```

**Key behaviors:**
- Reads TSV line-by-line (memory: ~1 line at a time)
- Buffers transitions until `group_id` changes
- Creates/deduplicates Peptide/Protein on first occurrence
- Calls callback with complete peptide group
- Caller processes and discards; reader continues

### Step 2: Create StreamingMRMAssay Class

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingMRMAssay.h`

Refactor `MRMAssay` operations to work on single peptide groups:

```cpp
class OPENMS_DLLAPI StreamingMRMAssay
{
public:
  /// Process a single peptide group through the full pipeline
  void processGroup(
    const String& group_id,
    std::vector<ReactionMonitoringTransition>& transitions,
    const TargetedExperiment::Peptide& peptide,
    // Output: filtered transitions ready for writing
    std::vector<ReactionMonitoringTransition>& output_transitions
  );

  /// Configuration (same as current MRMAssay parameters)
  void setParameters(
    double precursor_mz_threshold,
    double product_mz_threshold,
    int max_transitions,
    // ... other params
  );

private:
  // Reannotate transitions for this group
  void reannotateGroup_(std::vector<ReactionMonitoringTransition>& transitions,
                        const TargetedExperiment::Peptide& peptide);

  // Filter by m/z bounds
  void restrictGroup_(std::vector<ReactionMonitoringTransition>& transitions);

  // Select top N by intensity
  void selectTopTransitions_(std::vector<ReactionMonitoringTransition>& transitions,
                             int max_transitions);
};
```

### Step 3: Create StreamingTransitionTSVWriter Class

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.h`

```cpp
class OPENMS_DLLAPI StreamingTransitionTSVWriter
{
public:
  /// Open output file and write header
  void open(const String& filename);

  /// Write a batch of transitions (called per peptide group)
  void writeTransitions(const std::vector<ReactionMonitoringTransition>& transitions,
                        const Peptide& peptide,
                        const std::vector<Protein>& proteins);

  /// Finalize and close file
  void close();

private:
  std::ofstream ofs_;
  bool header_written_ = false;
};
```

### Step 4: Refactor OpenSwathAssayGenerator Main Flow

**File:** `src/topp/OpenSwathAssayGenerator.cpp`

```cpp
// NEW: Streaming mode for large files
if (use_streaming_mode)  // Could auto-detect based on file size
{
  StreamingTransitionTSVReader reader;
  StreamingTransitionTSVWriter writer;
  StreamingMRMAssay assay;

  assay.setParameters(precursor_mz_threshold, product_mz_threshold, ...);
  writer.open(out_tsv);

  reader.readGrouped(in_tsv, [&](const String& group_id,
                                  std::vector<TSVTransition>& tsv_transitions,
                                  const Peptide& peptide,
                                  const std::vector<Protein>& proteins)
  {
    // Convert TSVTransition -> ReactionMonitoringTransition (light version)
    std::vector<ReactionMonitoringTransition> transitions;
    convertTransitions_(tsv_transitions, transitions);

    // Process this group through the pipeline
    std::vector<ReactionMonitoringTransition> output;
    assay.processGroup(group_id, transitions, peptide, output);

    // Write results immediately
    writer.writeTransitions(output, peptide, proteins);

    // Memory freed when lambda exits
  });

  writer.close();
}
```

### Step 5: Handle Unsorted Input Files

TSV files may not be sorted by `group_id`. Two approaches:

**Option A: External Sort (Recommended for very large files)**
```cpp
// Pre-sort using external merge sort
void externalSortByGroupId(const String& input, const String& sorted_output);
```

**Option B: Two-Pass Approach**
```cpp
// Pass 1: Build index of group_id -> file offsets
std::map<String, std::vector<std::streampos>> group_offsets;

// Pass 2: For each group, seek and read all its transitions
for (auto& [group_id, offsets] : group_offsets) {
  std::vector<TSVTransition> group_transitions;
  for (auto offset : offsets) {
    file.seekg(offset);
    // read transition
  }
  callback(group_id, group_transitions, ...);
}
```

**Option C: Memory-bounded buffer with spill-to-disk**
```cpp
// Buffer up to N transitions in memory
// If group incomplete and buffer full, spill to temp file
// Process complete groups immediately
```

### Step 6: Create LightReactionMonitoringTransition

For streaming, we don't need full `CVTermList` inheritance. Create a lightweight variant:

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/LightReactionMonitoringTransition.h`

```cpp
struct OPENMS_DLLAPI LightReactionMonitoringTransition
{
  String name;
  String peptide_ref;
  double precursor_mz;
  double product_mz;
  int product_charge;
  double library_intensity;
  double rt_calibrated;
  String annotation;
  bool decoy;
  bool detecting;
  bool identifying;
  bool quantifying;

  // ~200 bytes vs ~1000 bytes for full RMT
};
```

## Memory Budget Analysis

| Component | Count | Size | Total |
|-----------|-------|------|-------|
| Peptide dedup map | 1M entries | 80 bytes | 80 MB |
| Protein dedup map | 100K entries | 80 bytes | 8 MB |
| Peptide objects | 1M | 500 bytes | 500 MB |
| Protein objects | 100K | 200 bytes | 20 MB |
| **Current group buffer** | ~1000 transitions | 200 bytes | **200 KB** |
| File I/O buffers | - | - | ~10 MB |
| **Total Peak Memory** | | | **~620 MB** |

Compared to current: **~15-20 GB → ~620 MB** (25-30x reduction)

## Files to Create/Modify

### New Files
1. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVReader.h`
2. `src/openms/source/ANALYSIS/OPENSWATH/StreamingTransitionTSVReader.cpp`
3. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.h`
4. `src/openms/source/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.cpp`
5. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingMRMAssay.h`
6. `src/openms/source/ANALYSIS/OPENSWATH/StreamingMRMAssay.cpp`
7. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/LightReactionMonitoringTransition.h`

### Modified Files
1. `src/topp/OpenSwathAssayGenerator.cpp` - Add streaming mode option
2. `src/openms/source/ANALYSIS/OPENSWATH/MRMAssay.cpp` - Extract per-group logic
3. `CMakeLists.txt` files - Add new sources

## Testing Strategy

1. **Unit tests** for each new class
2. **Integration test**: Compare streaming vs non-streaming output for small file (must be identical)
3. **Memory profiling**: Verify peak memory stays under 1GB for 10GB input
4. **Performance benchmark**: Streaming may be slower due to I/O; measure overhead

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Unsorted input files | High | Implement external sort or two-pass approach |
| IPF/UIS requires cross-group data | Medium | Buffer required groups or disable in streaming mode |
| PQP output needs all data | Medium | Stream to temp TSV, then convert to PQP |
| Different output between modes | High | Comprehensive diff testing |

## Implementation Order

1. **Phase 1**: `StreamingTransitionTSVReader` + basic callback
2. **Phase 2**: `StreamingTransitionTSVWriter` + end-to-end TSV→TSV
3. **Phase 3**: `StreamingMRMAssay` with reannotate/restrict/detect
4. **Phase 4**: Integration into `OpenSwathAssayGenerator`
5. **Phase 5**: Handle unsorted input (external sort)
6. **Phase 6**: Support PQP output format
7. **Phase 7**: Performance optimization and testing
