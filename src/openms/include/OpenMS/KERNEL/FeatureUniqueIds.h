// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class ConsensusMap;
  class FeatureMap;

  /**
    @brief Checks that every feature of a map has a valid unique id that no other feature of the map shares.

    The intended contract is that every Feature and ConsensusFeature in a map has a valid (non-zero) unique id
    that is unique within the map, subordinates included. Not every producer establishes it yet, and file writers
    do not enforce it; producers that do (e.g. FeatureFinderAlgorithmPicked, SeedListGenerator) verify their
    output with check() in their class tests. Feature handles are not checked.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI FeatureUniqueIds
  {
  public:
    /**
      @brief Throws Exception::InvalidValue if a feature or subordinate of @p map has an invalid id, or if two of them share an id

      The message names the number of invalid ids or the first duplicated id.
    */
    static void check(const FeatureMap& map);

    /// Same as above, for the consensus features of @p map (feature handles are not checked)
    static void check(const ConsensusMap& map);
  };
} // namespace OpenMS
