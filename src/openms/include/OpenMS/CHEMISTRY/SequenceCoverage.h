// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Development Team $
// $Authors: Srikanth K N $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief Compute sequence coverage of a protein by peptide sequences.

    This utility computes the fraction of amino acids in a protein sequence
    that are covered by a given set of peptide sequences.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI SequenceCoverage
  {
  public:
    /**
      @brief Compute sequence coverage.

      @param[in] protein Protein sequence.
      @param[in] peptides Peptide sequences.
      @return Coverage percentage in the range [0, 100].
    */
    static double getCoverage(
      const AASequence& protein,
      const std::vector<AASequence>& peptides
    );
  };

} // namespace OpenMS
