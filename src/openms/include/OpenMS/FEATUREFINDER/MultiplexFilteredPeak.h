// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteProfile.h>

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure storing a single peak that passed all filters
   * 
   * Each filter result corresponds to a successful search for a particular
   * peak pattern in the centroided data. The actual m/z shifts seen in the filter
   * result might differ from the theoretical shifts listed in the peak pattern.
   * 
   * Each MultiplexFilteredPeak consists of a primary peak and a set of satellite peaks.
   * The primary peak is a peak in the mono-isotopic masstrace of the lightest peptide
   * in the multiplet. The satellite peaks are peaks that form the m/z shift pattern
   * relative to the primary peak within a retention time range rt_band_. They are the
   * evidence on which grounds a peak may pass the filters.
   *
   * Note that in both centroid and profile mode, centroided data are filtered. (One of
   * the first steps in the profile mode algorithm is the peak picking of the profile
   * data.) Consequently in both modes, centroided peaks make up a final filtered peak.
   * @see size(). In profile mode, we additional store the profile data points that make
   * up these peak. @see sizeProfile().
   * 
   * @see MultiplexPeakPattern
   */
  class OPENMS_DLLAPI MultiplexFilteredPeak
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexFilteredPeak(double mz, float rt, size_t mz_idx, size_t rt_idx);

    /**
     * @brief returns m/z of the peak
     */
    double getMZ() const;
     
    /**
     * @brief returns RT of the peak
     */
    float getRT() const;
    
    /**
     * @brief returns the index of the peak in the spectrum
     */
    size_t getMZidx() const;
     
    /**
     * @brief returns the index of the corresponding spectrum in the MS experiment
     */
    size_t getRTidx() const;
    
    /**
     * @brief add a satellite peak
     */
    void addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx);
    
    void addSatellite(const MultiplexSatelliteCentroided& satellite, size_t pattern_idx);
    
    /**
     * @brief add a satellite data point
     */
    void addSatelliteProfile(float rt, double mz, float intensity, size_t pattern_idx);
    
    void addSatelliteProfile(const MultiplexSatelliteProfile& satellite, size_t pattern_idx);
    
    /**
     * @brief check if the peak (rt_idx, mz_idx) is already in the set of satellite peaks
     */
    bool checkSatellite(size_t rt_idx, size_t mz_idx) const;
    
    /**
     * @brief return all satellite peaks
     *
     * @see also satellites_
     */
    const std::multimap<size_t, MultiplexSatelliteCentroided >& getSatellites() const;
    
    /**
     * @brief return all satellite data points
     *
     * @see also satellites_profile_
     */
    const std::multimap<size_t, MultiplexSatelliteProfile >& getSatellitesProfile() const;
    
    /**
     * @brief return number of satellite peaks
     */
    size_t size() const;
    
    /**
     * @brief return number of satellite data points
     */
    size_t sizeProfile() const;
    
    private:
    /**
     * @brief position of the primary peak
     * 
     * Position of the primary peak in the m/z-RT plane in [Th, sec].
     * It is the input for the subsequent clustering step. 
     */
    double mz_;
    float rt_;
    
    /**
     * @brief indices of the primary peak position in the centroided experiment
     * 
     * Spectral index and peak index within the spectrum of the primary peak.
     * The indices are used to check the blacklist.
     */
    size_t mz_idx_;
    size_t rt_idx_;
    
    /**
     * @brief set of satellites
     * 
     * Mapping from a pattern index i.e. a specific mass trace to all peaks forming
     * the pattern. The primary peak is part of the satellite peak set.
     * 
     * pattern_idx -> (rt_idx, mz_idx)
     *
     * Typically peaks of the same mass trace show up in neighbouring spectra. The algorithm
     * considers spectra in the RT range @p rt_band. Consequently, the same @p pattern_idx key
     * will have multiple associated satellites, and a multimap is required.
     * 
     * Note that we store only indices, not iterators or pointers. We filter
     * 'white' experiments, but all indices refer to the original experiment.
     * White experiments are temporary (for each pattern), but the original
     * @p exp_picked_ experiment is permanent.
     */
    std::multimap<size_t, MultiplexSatelliteCentroided > satellites_;
    
    /**
     * @brief set of profile satellites (used on profile data only)
     *
     * Mapping from a pattern index i.e. a specific mass trace to all spline-interpolated
     * data points forming the pattern. Basically, when profile data are available as input,
     * we scan over the profile of each satellite peak (see MultiplexSatelliteCentroided above)
     * and decide if it passes the filters or not.
     *
     * pattern_idx -> (rt, mz, intensity)
     *
     * Typically peaks of the same mass trace show up in neighbouring spectra. The algorithm
     * considers spectra in the RT range @p rt_band. Consequently, the same @p pattern_idx key
     * will have multiple associated satellites, and a multimap is required.
     */
    std::multimap<size_t, MultiplexSatelliteProfile > satellites_profile_;
 
  };
  
}
