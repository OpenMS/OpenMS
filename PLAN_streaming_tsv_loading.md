# Plan: Streaming-Based TSV Loading for TargetedExperiment

## Problem Statement

Loading a 10GB TSV file (~30M transitions) into `TargetedExperiment` requires 15-20GB+ peak memory due to:
1. `ReactionMonitoringTransition` inherits `CVTermList` (~1000 bytes/transition vs 136 bytes for `LightTransition`)
2. All transitions loaded simultaneously before processing
3. Multiple processing passes each create new transition vectors

## Goal

Reduce peak memory from ~15-20GB to ~500MB-1GB by processing transitions in streaming fashion, grouped by peptide.

## Architecture Overview: Two-Pass Design

Decoy generation requires **global knowledge** of all target peptides to check for collisions. This necessitates a two-pass architecture:

```
                                    PASS 1: Index Building
┌─────────────┐    ┌──────────────────┐    ┌─────────────────────────────────────┐
│  TSV File   │───▶│ Sequential Scan  │───▶│ Build:                              │
└─────────────┘    └──────────────────┘    │ • Collision Index (seq+charge→id)   │
                                           │ • Group Order Index (group→offsets) │
                                           │ • Peptide/Protein dedup maps        │
                                           │ • Validate group consecutiveness    │
                                           └─────────────────────────────────────┘
                                                            │
                                                            ▼
                                    PASS 2: Streaming Processing
┌─────────────┐    ┌──────────────────┐    ┌─────────────────┐    ┌────────────┐
│  TSV File   │───▶│ Grouped Reader   │───▶│ Process Group:  │───▶│  Writer    │
└─────────────┘    └──────────────────┘    │ • Reannotate    │    └────────────┘
                          │                │ • Restrict      │
                          │                │ • Detect        │
                          ▼                │ • Gen Decoys    │
                   ┌──────────────────┐    └─────────────────┘
                   │ Collision Index  │◀──── Check & Update
                   │ (from Pass 1)    │
                   └──────────────────┘
```

## Key Insight: Why Two Passes?

### Decoy Collision Checking (MRMDecoy.cpp:419-497)

```cpp
// Collision index built from ALL targets before any decoy generation
std::unordered_map<std::string, std::string> allPeptideSequences;
for (const auto& peptide : all_peptides) {
  allPeptideSequences[getModifiedPeptideSequence_(peptide) +
                      String(peptide.getChargeState())] = peptide.id;
}

// When generating decoy, check against ALL targets + previously generated decoys
if (allPeptideSequences.find(decoy_key) != allPeptideSequences.end()) {
  // COLLISION - skip this decoy
}
```

**Collision key**: `modified_sequence + charge_state` (e.g., `"M(UniMod:4)PEPTIDEK+2"`)

### Transition-Level Collision (MRMAssay.cpp:686-694 for IPF/UIS)

```cpp
// IonMapT: SWATH_window → unmodified_sequence → [(fragment_mz, modified_sequence), ...]
IonMapT TargetIonMap;

// Check if decoy fragment m/z overlaps with any target fragment m/z
if (targetIonMap[swath][sequence].contains_mz_within_threshold(decoy_mz)) {
  // Skip overlapping decoy transition
}
```

## Implementation Steps

### Step 1: First Pass - Index Building

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVIndexer.h`

```cpp
class OPENMS_DLLAPI TransitionTSVIndexer
{
public:
  struct IndexResult {
    // Collision detection index: (modified_seq + charge) → peptide_id
    std::unordered_map<std::string, std::string> collision_index;

    // Group ordering: group_id → [file_offset, ...]
    // Empty if groups are consecutive (no reordering needed)
    std::map<String, std::vector<std::streampos>> group_offsets;

    // Deduplicated entities (lightweight versions for Pass 1)
    std::vector<LightPeptide> peptides;
    std::vector<LightProtein> proteins;
    std::vector<LightCompound> compounds;

    // Deduplication maps: id → index
    std::unordered_map<String, Size> peptide_map;
    std::unordered_map<String, Size> protein_map;
    std::unordered_map<String, Size> compound_map;

    // Statistics
    Size total_transitions = 0;
    Size total_groups = 0;
    bool groups_are_consecutive = true;  // False if reordering needed

    // For IPF/UIS: target ion map (SWATH → seq → [(mz, mod_seq), ...])
    IonMapT target_ion_map;
  };

  /// First pass: build all indices without storing transitions
  IndexResult buildIndex(const String& filename);

private:
  // Lightweight peptide for collision index (no CVTermList)
  struct LightPeptide {
    String id;
    String sequence;
    int charge;
    std::vector<std::pair<int, int>> modifications;  // (position, unimod_id)
  };

  String buildCollisionKey_(const LightPeptide& pep);
  void validateGroupConsecutiveness_(/* ... */);
};
```

**Pass 1 Operations:**

1. **Sequential scan** - Read each line once, O(n)
2. **Extract group_id** - Track current vs previous group
3. **Validate consecutiveness** - If group_id seen before after different group, mark `groups_are_consecutive = false`
4. **Build collision index** - For each unique peptide: `collision_index[seq+charge] = peptide_id`
5. **Record group offsets** - If not consecutive: `group_offsets[group_id].push_back(file_position)`
6. **Deduplicate** - Create lightweight Peptide/Protein/Compound objects
7. **Build target ion map** (if IPF enabled) - Store fragment m/z per SWATH window

**Memory for Pass 1** (10GB file, ~1M unique peptides):

| Component | Size |
|-----------|------|
| collision_index | ~100 MB (1M entries × 100 bytes) |
| group_offsets (if needed) | ~240 MB (30M offsets × 8 bytes) |
| peptide_map | ~50 MB |
| protein_map | ~5 MB |
| LightPeptide vector | ~100 MB |
| target_ion_map (IPF) | ~200 MB |
| **Total Pass 1** | **~500-700 MB** |

### Step 2: Handle Non-Consecutive Groups

If `groups_are_consecutive == false`, we have three options:

**Option A: External Sort (Recommended for very large files)**

```cpp
void externalSortByGroupId(const String& input, const String& sorted_output) {
  // 1. Split file into chunks that fit in memory
  // 2. Sort each chunk by group_id
  // 3. Merge sorted chunks using k-way merge
  // Result: sorted file where all transitions of a group are consecutive
}
```

**Option B: Use Group Offset Index**

```cpp
// In Pass 2, instead of sequential read:
for (auto& [group_id, offsets] : index.group_offsets) {
  std::vector<TSVTransition> group_transitions;
  for (auto offset : offsets) {
    file.seekg(offset);
    group_transitions.push_back(readLine(file));
  }
  processGroup(group_id, group_transitions);
}
```

**Tradeoff**: Random I/O is slower than sequential, but avoids creating temp file.

**Option C: Abort and Require Sorted Input**

```cpp
if (!index.groups_are_consecutive) {
  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "TSV file must have consecutive transitions per group for streaming mode. "
    "Use 'sort -t$'\\t' -k<group_col>' or disable streaming mode.", filename);
}
```

### Step 3: Second Pass - Streaming Processing

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionProcessor.h`

```cpp
class OPENMS_DLLAPI StreamingTransitionProcessor
{
public:
  using GroupCallback = std::function<void(
    const String& group_id,
    std::vector<LightTransition>& transitions,
    std::vector<LightTransition>& decoy_transitions  // Generated decoys
  )>;

  /// Configure with index from Pass 1
  void setIndex(TransitionTSVIndexer::IndexResult&& index);

  /// Configure processing parameters
  void setParameters(const Param& param);

  /// Second pass: stream-process with decoy generation
  void processFile(
    const String& filename,
    GroupCallback callback
  );

private:
  TransitionTSVIndexer::IndexResult index_;

  // Decoy generation state (updated during processing)
  std::unordered_map<std::string, std::string> decoy_collision_index_;

  // Process single group
  void processGroup_(
    const String& group_id,
    std::vector<LightTransition>& transitions,
    std::vector<LightTransition>& output_transitions,
    std::vector<LightTransition>& decoy_transitions
  );

  // Decoy generation with collision checking
  void generateDecoys_(
    const LightPeptide& peptide,
    const std::vector<LightTransition>& target_transitions,
    std::vector<LightTransition>& decoy_transitions
  );

  // Check collision against targets AND previously generated decoys
  bool checkCollision_(const String& collision_key) {
    return index_.collision_index.count(collision_key) > 0 ||
           decoy_collision_index_.count(collision_key) > 0;
  }
};
```

**Pass 2 Operations per Group:**

1. Read all transitions for current group
2. Convert to `LightTransition` objects
3. Reannotate (match to theoretical ion series)
4. Restrict (filter by m/z bounds)
5. Select top N (by intensity)
6. **Generate decoys** using collision index from Pass 1
7. Check each decoy against `collision_index` (targets) AND `decoy_collision_index` (previous decoys)
8. Write targets + decoys to output
9. Free memory, move to next group

### Step 4: Streaming Writer

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.h`

```cpp
class OPENMS_DLLAPI StreamingTransitionTSVWriter
{
public:
  void open(const String& filename);

  /// Write target and decoy transitions for one group
  void writeGroup(
    const String& group_id,
    const std::vector<LightTransition>& targets,
    const std::vector<LightTransition>& decoys,
    const LightPeptide& peptide,
    const std::vector<LightProtein>& proteins
  );

  void close();

private:
  std::ofstream ofs_;
  bool header_written_ = false;
};
```

### Step 5: Integration - OpenSwathAssayGenerator

**File:** `src/topp/OpenSwathAssayGenerator.cpp`

```cpp
// Streaming mode for large files
if (use_streaming_mode || file_size > streaming_threshold)
{
  OPENMS_LOG_INFO << "Using streaming mode for large file..." << std::endl;

  // === PASS 1: Build indices ===
  OPENMS_LOG_INFO << "Pass 1: Building collision index and validating group order..." << std::endl;
  TransitionTSVIndexer indexer;
  auto index = indexer.buildIndex(in_tsv);

  OPENMS_LOG_INFO << "Found " << index.total_transitions << " transitions in "
                  << index.total_groups << " groups" << std::endl;
  OPENMS_LOG_INFO << "Collision index size: " << index.collision_index.size() << std::endl;

  // Handle non-consecutive groups
  if (!index.groups_are_consecutive) {
    if (allow_external_sort) {
      OPENMS_LOG_INFO << "Groups not consecutive, performing external sort..." << std::endl;
      String sorted_file = File::getTempDirectory() + "/" + File::getUniqueName() + ".tsv";
      externalSortByGroupId(in_tsv, sorted_file);
      in_tsv = sorted_file;
      // Rebuild index for sorted file (groups now consecutive)
      index = indexer.buildIndex(in_tsv);
    } else {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Groups not consecutive. Enable -allow_external_sort or pre-sort input.", in_tsv);
    }
  }

  // === PASS 2: Stream processing ===
  OPENMS_LOG_INFO << "Pass 2: Processing transitions with decoy generation..." << std::endl;
  StreamingTransitionProcessor processor;
  processor.setIndex(std::move(index));
  processor.setParameters(getParam_());

  StreamingTransitionTSVWriter writer;
  writer.open(out_tsv);

  Size groups_processed = 0;
  processor.processFile(in_tsv, [&](const String& group_id,
                                     std::vector<LightTransition>& targets,
                                     std::vector<LightTransition>& decoys)
  {
    writer.writeGroup(group_id, targets, decoys,
                      processor.getPeptide(group_id),
                      processor.getProteins(group_id));

    if (++groups_processed % 10000 == 0) {
      OPENMS_LOG_INFO << "Processed " << groups_processed << " groups..." << std::endl;
    }
  });

  writer.close();
  OPENMS_LOG_INFO << "Streaming complete. Processed " << groups_processed << " groups." << std::endl;
}
```

### Step 6: Light Data Structures

**File:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/LightTargetedTypes.h`

```cpp
/// Lightweight transition (~136 bytes vs ~1000 bytes for ReactionMonitoringTransition)
struct OPENMS_DLLAPI LightTransition
{
  String name;                    // 32 bytes
  String peptide_ref;             // 32 bytes
  double precursor_mz;            // 8 bytes
  double product_mz;              // 8 bytes
  double library_intensity;       // 8 bytes
  double rt_calibrated;           // 8 bytes
  int precursor_charge;           // 4 bytes
  int product_charge;             // 4 bytes
  String annotation;              // 32 bytes
  uint8_t flags;                  // 1 byte (decoy, detecting, identifying, quantifying)

  bool isDecoy() const { return flags & 0x01; }
  bool isDetecting() const { return flags & 0x02; }
  bool isIdentifying() const { return flags & 0x04; }
  bool isQuantifying() const { return flags & 0x08; }
};

/// Lightweight peptide for collision checking (~100 bytes vs ~500 bytes)
struct OPENMS_DLLAPI LightPeptide
{
  String id;
  String sequence;
  String protein_refs;            // Comma-separated for compactness
  int charge;
  double rt;
  std::vector<std::pair<int16_t, int32_t>> mods;  // (position, unimod_id)

  String getCollisionKey() const;
};

/// Lightweight protein (~50 bytes vs ~200 bytes)
struct OPENMS_DLLAPI LightProtein
{
  String id;
  String accession;
};
```

## Memory Budget Analysis (Two-Pass)

### Pass 1 Peak Memory

| Component | Count | Size | Total |
|-----------|-------|------|-------|
| Collision index | 1M entries | 100 bytes | 100 MB |
| Group offsets (worst case) | 30M entries | 8 bytes | 240 MB |
| LightPeptide vector | 1M | 100 bytes | 100 MB |
| LightProtein vector | 100K | 50 bytes | 5 MB |
| Dedup maps | 1.1M entries | 50 bytes | 55 MB |
| Target ion map (IPF) | 1M × 5 ions | 40 bytes | 200 MB |
| **Pass 1 Total** | | | **~700 MB** |

### Pass 2 Peak Memory

| Component | Count | Size | Total |
|-----------|-------|------|-------|
| Index (from Pass 1) | - | - | ~500 MB |
| Decoy collision index | ~1M entries | 100 bytes | 100 MB |
| Current group buffer | ~1000 trans | 136 bytes | 136 KB |
| File I/O buffers | - | - | ~10 MB |
| **Pass 2 Total** | | | **~620 MB** |

### Comparison

| Mode | Peak Memory |
|------|-------------|
| Current (load all) | 15-20 GB |
| **Two-Pass Streaming** | **~700 MB** |

**20-30x reduction** in peak memory usage.

## Files to Create/Modify

### New Files
1. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVIndexer.h`
2. `src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVIndexer.cpp`
3. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionProcessor.h`
4. `src/openms/source/ANALYSIS/OPENSWATH/StreamingTransitionProcessor.cpp`
5. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.h`
6. `src/openms/source/ANALYSIS/OPENSWATH/StreamingTransitionTSVWriter.cpp`
7. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/LightTargetedTypes.h`
8. `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ExternalSort.h` (optional)

### Modified Files
1. `src/topp/OpenSwathAssayGenerator.cpp` - Add streaming mode option
2. `src/openms/source/ANALYSIS/OPENSWATH/MRMDecoy.cpp` - Extract collision key generation
3. `CMakeLists.txt` files - Add new sources

## Testing Strategy

1. **Unit tests** for each new class
2. **Collision correctness**: Verify decoy collision detection matches non-streaming mode
3. **Integration test**: Compare streaming vs non-streaming output (must be identical)
4. **Memory profiling**: Verify peak memory stays under 1GB for 10GB input
5. **Non-consecutive input**: Test with shuffled group order
6. **Performance benchmark**: Measure two-pass overhead vs single-pass

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Non-consecutive groups | High | External sort or offset-based random access |
| IPF/UIS cross-group collision | High | Build target_ion_map in Pass 1 |
| Decoy order dependence | Medium | Process groups in deterministic order |
| PQP output needs final data | Medium | Stream to temp TSV, then convert |
| Two-pass I/O overhead | Low | Pass 1 is lightweight; sequential I/O is fast |

## Implementation Order

1. **Phase 1**: `LightTargetedTypes.h` - Lightweight data structures
2. **Phase 2**: `TransitionTSVIndexer` - Pass 1 index building
3. **Phase 3**: Group consecutiveness validation + external sort
4. **Phase 4**: `StreamingTransitionProcessor` - Pass 2 processing
5. **Phase 5**: Decoy generation with collision checking
6. **Phase 6**: `StreamingTransitionTSVWriter` - Output streaming
7. **Phase 7**: Integration into `OpenSwathAssayGenerator`
8. **Phase 8**: IPF/UIS support with target_ion_map
9. **Phase 9**: PQP output format support
10. **Phase 10**: Performance optimization and comprehensive testing

## Command-Line Interface

```bash
# Auto-detect streaming mode based on file size
OpenSwathAssayGenerator -in large_library.tsv -out filtered.tsv

# Force streaming mode
OpenSwathAssayGenerator -in library.tsv -out filtered.tsv -streaming_mode

# With external sort for unsorted input
OpenSwathAssayGenerator -in unsorted.tsv -out filtered.tsv -streaming_mode -allow_external_sort

# Set streaming threshold (default: 1GB)
OpenSwathAssayGenerator -in library.tsv -out filtered.tsv -streaming_threshold 500000000
```
