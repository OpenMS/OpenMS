// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineInterpolatedPeaks.h>

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
     * @param satellites    spline-interpolated satellites of the peak. If they pass, they will be added to the peak.
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than <averagine_similarity_>
     */
    bool filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const;

    /**
     * @brief peptide correlation filter for profile mode
     *
     * @param pattern    m/z pattern to search for
     * @param satellites    spline-interpolated satellites of the peak. If they pass, they will be added to the peak.
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than <averagine_similarity_>
     */
    bool filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const;
    
    /**
     * @brief spline interpolated profile data and peak boundaries
     */
    std::vector<SplineInterpolatedPeaks> exp_spline_profile_;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_;

  };

}

