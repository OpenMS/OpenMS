// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief filters centroided and profile data for peak patterns
   *
   * The algorithm searches for patterns of multiple peptides in the data.
   * The peptides appear as characteristic patterns of isotopic peaks in
   * MS1 spectra. We first search the centroided data, and optionally in
   * a second step the spline interpolated profile data. For each
   * peak pattern the algorithm generates a filter result.
   *
   * @see MultiplexIsotopicPeakPattern
   * @see MultiplexFilterResult
   * @see MultiplexFiltering
   */
  class OPENMS_DLLAPI MultiplexFilteringProfile :
    public MultiplexFiltering
  {
public:
    /**
     * @brief constructor
     *
     * @param exp_profile    experimental data in profile mode
     * @param exp_centroided    experimental data in centroid mode
     * @param boundaries    peak boundaries for exp_centroided
     * @param patterns    patterns of isotopic peaks to be searched for
     * @param isotopes_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param isotopes_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param intensity_cutoff    intensity cutoff
     * @param rt_band    RT range used for filtering
     * @param mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param peptide_similarity    similarity score for two peptides in the same multiplet
     * @param averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p).
     * @param averagine_type    The averagine model to use, current options are RNA DNA or peptide.
     *
     * @throw Exception::IllegalArgument if profile and centroided data do not contain same number of spectra
     * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
     */
    MultiplexFilteringProfile(MSExperiment& exp_profile, const MSExperiment& exp_centroided, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries,
                              const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band,
                              double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type="peptide");

    /**
     * @brief filter for patterns
     * (generates a filter result for each of the patterns)
     *
     * @throw Exception::IllegalArgument if number of peaks and number of peak boundaries differ
     *
     * @see MultiplexIsotopicPeakPattern
     * @see MultiplexFilteredMSExperiment
     */
    std::vector<MultiplexFilteredMSExperiment> filter();

    /**
     * @brief returns the intensity-filtered peak boundaries
     */
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& getPeakBoundaries();

private:
    /**
     * @brief averagine filter for profile mode
     *
     * @param pattern    m/z pattern to search for
     * @param peak    peak to be filtered
     * @param satellites_profile    spline-interpolated satellites of the peak. If they pass, they will be added to the peak.
     *
     * @return if this filter was passed i.e. the correlation coefficient is greater than averagine_similarity_
     */
    bool filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const;

    /**
     * @brief peptide correlation filter for profile mode
     *
     * @param pattern    m/z pattern to search for
     * @param satellites_profile    spline-interpolated satellites of the peak. If they pass, they will be added to the peak.
     *
     * @return if this filter was passed i.e. the correlation coefficient is greater than averagine_similarity_
     */
    bool filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const;
    
    /**
     * @brief spline interpolated profile data and peak boundaries
     */
    std::vector<SplineInterpolatedPeaks> exp_spline_profile_;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_;

  };

}

