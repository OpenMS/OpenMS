// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <map>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

namespace OpenMS
{
  /**
   * @brief Map a protein position to all observed modifications and associated statistics 
   * 
   * For example, allows storing that position 10 in the protein carries Oxidation (M) 
   * and was observed in 123 PSMs.
   * Note: Implementation uses a std::map, thus accessing a location not present in the map 
   * with operator[] will value construct an empty map at that location. Like with std::maps
   * check first if the key exists.
   */
  struct OPENMS_DLLAPI ProteinModificationSummary
  {
    /// basic modification statistic
    struct OPENMS_DLLAPI Statistics
    {
      bool operator==(const Statistics& rhs) const;
      size_t count = 0;  ///< total number of PSMs supporting the modification at this position
      double frequency = 0.0; ///< PSMs with modification / total number of PSMs
      double FLR = 0.0; ///< false localization rate
      double probability = 0.0; ///< (localization) probability
    };

    /// comparison operator
    bool operator==(const ProteinModificationSummary& rhs) const;

    using ModificationsToStatistics = std::map<const ResidueModification*, Statistics>;
    using AALevelModificationSummary = std::map<size_t, ModificationsToStatistics>;

    /// position -> modification -> statistic (counts, etc.)
    AALevelModificationSummary AALevelSummary;
  };
}

