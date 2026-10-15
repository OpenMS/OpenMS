// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/IMSExperimentRandomAccess.h>

namespace OpenMS
{
  IMSExperimentRandomAccess::~IMSExperimentRandomAccess() = default;

  MSSpectrum IMSExperimentRandomAccess::getSpectrumPeaks(Size index)
  {
    return getSpectrum(index);
  }

  Size IMSExperimentRandomAccess::findFirstRT(double rt) const
  {
    const Size n = getNrSpectra();
    if (isSortedByRT())
    {
      Size lo = 0, hi = n;
      while (lo < hi)
      {
        const Size mid = lo + (hi - lo) / 2;
        if (getSpectrumRT(mid) < rt) { lo = mid + 1; }
        else { hi = mid; }
      }
      return lo;
    }
    for (Size i = 0; i < n; ++i)
    {
      if (getSpectrumRT(i) >= rt) { return i; }
    }
    return n;
  }
} // namespace OpenMS
