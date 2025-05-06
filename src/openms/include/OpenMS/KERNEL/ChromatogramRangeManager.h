// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Administrator $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/RangeManager.h>

namespace OpenMS
{
  /**
    @brief Range manager for chromatograms
    
    This class manages retention time, mz, and intensity ranges for multiple chromatograms.
    
    @ingroup Kernel
  */
  class OPENMS_DLLAPI ChromatogramRangeManager : public RangeManager<RangeRT, RangeIntensity, RangeMZ>
  {
  public:
    /// Base type
    using BaseType = RangeManager<RangeRT, RangeIntensity, RangeMZ>;

    /// Default constructor
    ChromatogramRangeManager() = default;
    
    /// Copy constructor
    ChromatogramRangeManager(const ChromatogramRangeManager& source) = default;
    
    /// Move constructor
    ChromatogramRangeManager(ChromatogramRangeManager&& source) = default;
    
    /// Assignment operator
    ChromatogramRangeManager& operator=(const ChromatogramRangeManager& source) = default;
    
    /// Move assignment operator
    ChromatogramRangeManager& operator=(ChromatogramRangeManager&& source) = default;
    
    /// Destructor
    ~ChromatogramRangeManager() = default;
  };

} // namespace OpenMS