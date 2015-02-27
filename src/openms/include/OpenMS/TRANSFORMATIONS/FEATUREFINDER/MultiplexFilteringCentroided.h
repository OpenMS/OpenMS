// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERINGCENTROIDED_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERINGCENTROIDED_H

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief filters centroided data for peak patterns
   *
   * The algorithm searches for patterns of multiple peptides in the data.
   * The peptides appear as characteristic patterns of isotopic peaks in
   * MS1 spectra. We search the centroided data for such patterns.
   * For each peak pattern the algorithm generates a filter result.
   *
   * @see MultiplexPeakPattern
   * @see MultiplexFilterResult
   * @see MultiplexFiltering
   */
  class OPENMS_DLLAPI MultiplexFilteringCentroided :
    public MultiplexFiltering
  {
public:
    /**
     * @brief constructor
     *
     * @param exp_picked    experimental data in centroid mode
     * @param patterns    patterns of isotopic peaks to be searched for
     * @param peaks_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param peaks_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param missing_peaks    flag for missing peaks
     * @param intensity_cutoff    intensity cutoff
     * @param mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param peptide_similarity    similarity score for two peptides in the same multiplet
     * @param averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p). 
     */
    MultiplexFilteringCentroided(const MSExperiment<Peak1D>& exp_picked, const std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling);

    /**
     * @brief filter for patterns
     * (generates a filter result for each of the patterns)
     *
     * @see MultiplexPeakPattern, MultiplexFilterResult
     */
    std::vector<MultiplexFilterResult> filter();

private:
    /**
     * @brief non-local intensity filter
     *
     * Checks if the intensities at the pattern positions are above the intensity cutoff.
     * We check not only at m/z but at all pattern positions i.e. non-locally.
     * (In filter 1 we checked that peaks do exist at these positions.
     *  In filter 2 we checked that the mono-isotopic peak intensities are above the threshold.)
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param spectrum_index    index of the spectrum in exp_picked_ and boundaries_
     * @param mz_shifts_actual_indices    indices of peaks corresponding to the pattern
     * @param intensities_actual    output for the spline-interpolated intensities at the actual m/z shift positions
     * @param peaks_found_in_all_peptides    number of isotopic peaks seen for each peptide (peaks)
     *
     * @return number of isotopic peaks seen for each peptide (profile)
     */
    int nonLocalIntensityFilter(MultiplexPeakPattern pattern, int spectrum_index, const std::vector<int>& mz_shifts_actual_indices, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides) const;

    /**
     * @brief returns the index of a peak at m/z
     * (for initialisation of peak registry)
     *
     * @param spectrum_index    index of the spectrum in exp_picked_ and boundaries_
     * @param mz    m/z position of the peak
     * @param scaling    rescaling of the peak boundaries
     *
     * @return index of the peak in spectrum
     */
    int getPeakIndex(int spectrum_index, double mz, double scaling) const;

  };

}

#endif /* MULTIPLEXFILTERINGCENTROIDED_H_ */
