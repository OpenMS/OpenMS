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

    Contract: every Feature and ConsensusFeature in a map that a library algorithm returns has a valid
    (non-zero) unique id, and the ids are unique within the map, subordinates included. Algorithms that
    produce feature maps establish it before returning; class tests verify it with check().

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
