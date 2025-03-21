// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>

namespace OpenMS
{

  /**
     @brief Class for reading Percolator tab-delimited output files.

     For PSM-level output, the file extension should be ".psms".
  */
  class OPENMS_DLLAPI PercolatorOutfile
  {

  public:

    /// Types of Percolator scores
    enum ScoreType { QVALUE, POSTERRPROB, SCORE, SIZE_OF_SCORETYPE };

    /// Names of Percolator scores (to match ScoreType)
    static const std::string score_type_names[SIZE_OF_SCORETYPE];

    /// Return a score type given its name
    static enum ScoreType getScoreType(String score_type_name);

    /// Constructor
    PercolatorOutfile();

    /// Loads a Percolator output file
    void load(const String& filename, ProteinIdentification& proteins,
              std::vector<PeptideIdentification>& peptides,
              SpectrumMetaDataLookup& lookup,
              enum ScoreType output_score = QVALUE);

  private:
    /// Converts the peptide string to an 'AASequence' instance
    void getPeptideSequence_(String peptide, AASequence& seq) const;

    /// Resolve cases where N-terminal modifications may be misassigned to the first residue (for X! Tandem results)
    void resolveMisassignedNTermMods_(String& peptide) const;
  };

} // namespace OpenMS

