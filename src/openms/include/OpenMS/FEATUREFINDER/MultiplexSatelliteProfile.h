// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  /**
   * @brief data structure storing a single satellite data point
   *
   * The satellite data point is a spline-interpolated point of profile MSExperiment.
   * The triplet of RT, m/z and intensity is therefore stored explicitly.
   * 
   * @see MultiplexFilteredPeak, MultiplexSatelliteCentroided
   */
  class OPENMS_DLLAPI MultiplexSatelliteProfile
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexSatelliteProfile(float rt, double mz, float intensity);
    
    /**
     * @brief returns the RT of the satellite data point
     */
    float getRT() const;
    
    /**
     * @brief returns the m/z of the satellite data point
     */
    double getMZ() const;
    
    /**
     * @brief returns the intensity of the satellite data point
     */
    float getIntensity() const;
    
    private:
     
    /**
     * @brief position and intensity of the data point within the spline-interpolated experiment
     */
    float rt_;
    double mz_;
    float intensity_;
    
  };
  
}

