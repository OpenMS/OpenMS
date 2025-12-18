// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/StreamingTSVReader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/StreamingTSVWriter.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <unordered_map>
#include <random>

namespace OpenMS
{

  /**
    @brief Streaming processor for transition libraries.

    Processes large transition TSV files with minimal memory usage by:
    1. Building an index (collision detection, group validation) in Pass 1
    2. Processing each peptide group in Pass 2:
       - Filtering transitions (m/z bounds, intensity)
       - Selecting top N transitions per group
       - Generating decoys with collision checking
       - Writing output immediately

    @ingroup Analysis
  */
  class OPENMS_DLLAPI StreamingTransitionProcessor :
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:

    /// Constructor
    StreamingTransitionProcessor();

    /// Destructor
    ~StreamingTransitionProcessor() override;

    /**
      @brief Process a TSV file in streaming mode.

      Two-pass processing:
      - Pass 1: Build collision index, validate group ordering
      - Pass 2: Stream-process groups with filtering and decoy generation

      @param input_file Path to input TSV file
      @param output_file Path to output TSV file
    */
    void process(const String& input_file, const String& output_file);

    /**
      @brief Get statistics from the last processing run.
    */
    struct Statistics
    {
      size_t input_transitions = 0;
      size_t input_groups = 0;
      size_t output_transitions = 0;
      size_t output_groups = 0;
      size_t decoys_generated = 0;
      size_t decoys_collided = 0;
    };

    const Statistics& getStatistics() const { return stats_; }

  protected:

    /// Update parameters
    void updateMembers_() override;

  private:

    /**
      @brief Process a single peptide group.

      Applies filtering and generates decoys for one group.

      @param group_id Peptide group identifier
      @param transitions Input transitions for this group
      @param compound The compound (peptide/metabolite)
      @param protein_ids Associated protein identifiers
      @param collision_index Index for checking decoy collisions
      @param output_transitions Output: filtered target transitions
      @param output_decoys Output: generated decoy transitions
    */
    void processGroup_(
      const std::string& group_id,
      std::vector<OpenSwath::LightTransition>& transitions,
      const OpenSwath::LightCompound& compound,
      const std::vector<std::string>& protein_ids,
      const std::unordered_map<std::string, std::string>& collision_index,
      std::vector<OpenSwath::LightTransition>& output_transitions,
      std::vector<OpenSwath::LightTransition>& output_decoys
    );

    /**
      @brief Filter transitions by m/z and intensity bounds.
    */
    void filterTransitions_(
      std::vector<OpenSwath::LightTransition>& transitions
    );

    /**
      @brief Select top N transitions by intensity.
    */
    void selectTopTransitions_(
      std::vector<OpenSwath::LightTransition>& transitions,
      size_t max_transitions
    );

    /**
      @brief Generate decoy transitions for a group.

      Uses pseudo-reverse method by default: reverses sequence while
      keeping C-terminal amino acid fixed (for tryptic peptides).

      @param compound Target compound
      @param target_transitions Target transitions to generate decoys from
      @param collision_index Index to check for collisions
      @param decoy_transitions Output: generated decoy transitions
      @return Number of collisions detected
    */
    size_t generateDecoys_(
      const OpenSwath::LightCompound& compound,
      const std::vector<OpenSwath::LightTransition>& target_transitions,
      const std::unordered_map<std::string, std::string>& collision_index,
      std::vector<OpenSwath::LightTransition>& decoy_transitions
    );

    /**
      @brief Create pseudo-reverse decoy sequence.

      Reverses the sequence except for the C-terminal amino acid
      (to maintain tryptic cleavage properties).
    */
    std::string pseudoReverseSequence_(const std::string& sequence);

    /**
      @brief Check if a decoy collides with targets or existing decoys.
    */
    bool checkCollision_(
      const std::string& collision_key,
      const std::unordered_map<std::string, std::string>& target_index
    );

    // Parameters
    double product_mz_threshold_ = 0.05;
    double precursor_mz_threshold_ = 0.05;
    double min_product_mz_ = 350.0;
    double max_product_mz_ = 2000.0;
    int max_transitions_ = 6;
    bool generate_decoys_ = true;
    std::string decoy_tag_ = "DECOY_";
    double precursor_mz_shift_ = 0.0;
    double product_mz_shift_ = 20.0;

    // State
    Statistics stats_;
    std::unordered_map<std::string, std::string> decoy_collision_index_;
    std::mt19937 rng_;
  };

} // namespace OpenMS
