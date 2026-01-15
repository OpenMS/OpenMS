// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure storing a single satellite peak
   *
   * The satellite peak is part of a centroided MSExperiment.
   * Hence indices rt_idx_ and mz_idx_ are sufficient to specify RT, m/z and intensity.
   * 
   * @see MultiplexFilteredPeak, MultiplexSatelliteProfile
   */
  class OPENMS_DLLAPI MultiplexSatelliteCentroided
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexSatelliteCentroided(size_t rt_idx, size_t mz_idx);
    
    /**
     * @brief returns the m/z index of the satellite peak
     */
    size_t getMZidx() const;
     
    /**
     * @brief returns the RT index of the satellite peak
     */
    size_t getRTidx() const;
    
    private:
     
    /**
     * @brief indices of the satellite peak position in the centroided experiment
     * 
     * Spectral index and peak index within the spectrum of the satellite peak.
     */
    size_t rt_idx_;
    size_t mz_idx_;
    
  };
  
}

