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
    
    This class manages retention time, m/z, and intensity ranges for multiple chromatograms.
    It extends the basic RangeManager to provide specialized functionality for chromatogram data.
    
    The ChromatogramRangeManager is used in conjunction with the SpectrumRangeManager in MSExperiment
    to provide separate range tracking for chromatograms and spectra. This separation allows for
    more efficient and targeted range operations on specific data types.
    
    @see RangeManager
    @see SpectrumRangeManager
    @see MSExperiment
    @ingroup Kernel
  */
  class OPENMS_DLLAPI ChromatogramRangeManager : public RangeManager<RangeRT, RangeIntensity, RangeMZ>
  {
  public:
    /// Base type
    using BaseType = RangeManager<RangeRT, RangeIntensity, RangeMZ>;

  };

} // namespace OpenMS