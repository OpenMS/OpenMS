// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <numeric>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Implementation detail of the sort(Predicate) members of MSSpectrum, MSChromatogram and Mobilogram.

      Returns the indices [0, @p n) stably sorted by @p lambda (<tt>[](Size index_1, Size index_2) --> bool</tt>).
      The result is a permutation (in range and duplicate-free by construction), suitable for selectUnchecked().
      Not part of the public API.
    */
    template<class Predicate>
    std::vector<Size> stableSortedIndices(Size n, const Predicate& lambda)
    {
      std::vector<Size> indices(n);
      std::iota(indices.begin(), indices.end(), 0);
      std::stable_sort(indices.begin(), indices.end(), lambda);
      return indices;
    }
  } // namespace Internal
} // namespace OpenMS
