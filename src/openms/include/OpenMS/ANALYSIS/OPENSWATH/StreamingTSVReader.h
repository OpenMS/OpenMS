// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <functional>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{

  /**
    @brief Memory-efficient streaming reader for large transition TSV files.

    This class implements a two-pass streaming approach to load transition libraries
    with minimal memory footprint. Instead of loading all transitions into memory,
    it processes them group-by-group.

    Pass 1 (buildIndex): Scans the file to:
    - Build collision detection index for decoy generation
    - Validate that transitions are grouped consecutively by peptide
    - Deduplicate peptides, proteins, and compounds

    Pass 2 (processGroups): Streams through the file calling a callback for each
    complete peptide group, allowing immediate processing and output.

    @note Requires transitions of the same group to appear consecutively in the file.
          If groups are split, an exception is thrown with instructions to pre-sort.

    @ingroup Analysis
  */
  class OPENMS_DLLAPI StreamingTSVReader :
    public ProgressLogger
  {
  public:

    /// Result of the first pass index building
    struct IndexResult
    {
      /// Collision detection index: (modified_seq + charge) -> peptide_id
      std::unordered_map<std::string, std::string> collision_index;

      /// Deduplicated compounds (peptides/metabolites)
      std::vector<OpenSwath::LightCompound> compounds;

      /// Deduplicated proteins
      std::vector<OpenSwath::LightProtein> proteins;

      /// Deduplication maps: id -> index in vectors
      std::unordered_map<std::string, size_t> compound_map;
      std::unordered_map<std::string, size_t> protein_map;

      /// Statistics
      size_t total_transitions = 0;
      size_t total_groups = 0;
    };

    /// Callback signature for processing a complete peptide group
    using GroupCallback = std::function<void(
      const std::string& group_id,
      std::vector<OpenSwath::LightTransition>& transitions,
      const OpenSwath::LightCompound& compound,
      const std::vector<std::string>& protein_ids
    )>;

    /// Constructor
    StreamingTSVReader();

    /// Destructor
    ~StreamingTSVReader() override;

    /**
      @brief First pass: build indices without storing transitions.

      Scans the entire TSV file to:
      1. Build collision index for decoy generation
      2. Validate consecutive group ordering
      3. Deduplicate compounds and proteins

      @param filename Path to the TSV file
      @return IndexResult containing all indices and deduplicated entities
      @throws Exception::InvalidValue if groups are not consecutive
    */
    IndexResult buildIndex(const String& filename);

    /**
      @brief Second pass: stream-process groups with callback.

      Reads the file group-by-group, calling the callback for each complete
      peptide group. Memory usage is bounded by the size of a single group.

      @param filename Path to the TSV file
      @param index The index built by buildIndex()
      @param callback Function to call for each complete group
    */
    void processGroups(
      const String& filename,
      const IndexResult& index,
      GroupCallback callback
    );

    /**
      @brief Build collision key for a compound (for decoy checking).

      Format: modified_sequence + charge_state
      Example: "M(UniMod:4)PEPTIDEK+2"

      @param compound The compound to build the key for
      @return Collision key string
    */
    static std::string buildCollisionKey(const OpenSwath::LightCompound& compound);

  private:

    /// Parse TSV header and determine column positions
    void parseHeader_(const std::string& line, char& delimiter,
                      std::map<std::string, int>& header_dict) const;

    /// Parse a single TSV line into transition data
    void parseLine_(const std::string& line, char delimiter,
                    const std::map<std::string, int>& header_dict,
                    OpenSwath::LightTransition& transition,
                    std::string& group_id,
                    std::string& sequence,
                    std::string& full_peptide_name,
                    int& charge,
                    std::vector<std::string>& protein_names,
                    std::string& compound_name,
                    double& rt) const;

    /// Build modified sequence string for collision key
    static std::string buildModifiedSequence_(const OpenSwath::LightCompound& compound);

    /// Known header names for TSV columns
    static const std::vector<std::string> header_names_;
  };

} // namespace OpenMS
