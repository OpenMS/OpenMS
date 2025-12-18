// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>
#include <vector>

namespace OpenMS
{

  /**
    @brief Streaming writer for transition TSV files.

    Writes transitions incrementally, group-by-group, without requiring
    all data to be in memory at once.

    @ingroup Analysis
  */
  class OPENMS_DLLAPI StreamingTSVWriter
  {
  public:

    /// Constructor
    StreamingTSVWriter();

    /// Destructor
    ~StreamingTSVWriter();

    /**
      @brief Open output file and write header.

      @param filename Path to output TSV file
    */
    void open(const String& filename);

    /**
      @brief Write transitions for one peptide group.

      @param group_id The peptide group identifier
      @param transitions Target transitions for this group
      @param decoys Decoy transitions for this group
      @param compound The compound (peptide/metabolite) for this group
      @param proteins Protein identifiers for this group
    */
    void writeGroup(
      const std::string& group_id,
      const std::vector<OpenSwath::LightTransition>& transitions,
      const std::vector<OpenSwath::LightTransition>& decoys,
      const OpenSwath::LightCompound& compound,
      const std::vector<std::string>& protein_ids
    );

    /**
      @brief Finalize and close the output file.
    */
    void close();

    /// Get number of transitions written
    size_t getTransitionsWritten() const { return transitions_written_; }

    /// Get number of groups written
    size_t getGroupsWritten() const { return groups_written_; }

  private:

    /// Write the TSV header
    void writeHeader_();

    /// Write a single transition line
    void writeTransition_(
      const OpenSwath::LightTransition& transition,
      const OpenSwath::LightCompound& compound,
      const std::vector<std::string>& protein_ids
    );

    std::ofstream ofs_;
    bool header_written_ = false;
    size_t transitions_written_ = 0;
    size_t groups_written_ = 0;
  };

} // namespace OpenMS
