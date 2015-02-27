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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERING_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERING_H

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief base class for filtering centroided and profile data for peak patterns
   *
   * The algorithm searches for patterns of multiple peptides in the data.
   * The peptides appear as characteristic patterns of isotopic peaks in
   * MS1 spectra. We first search the centroided data, and optionally in
   * a second step the spline interpolated profile data. For each
   * peak pattern the algorithm generates a filter result.
   *
   * The algorithm differs slightly for centroided and profile input data.
   * This base class comprises code common to both. The two child classes
   * MultiplexFilteringCentroided and MultiplexFilteringProfile contain
   * specific functions and the primary filter() method.
   *
   * @see MultiplexPeakPattern
   * @see MultiplexFilterResult
   * @see MultiplexFilteringCentroided
   * @see MultiplexFilteringProfile
   */
  class OPENMS_DLLAPI MultiplexFiltering :
    public ProgressLogger
  {
public:
    /**
     * @brief structure for peak position in neighbouring spectra
     *
     * Each peak can be annotated with a peak reference. The reference
     * contains the indices of peaks with same m/z in the previous and
     * next spectrum. If none is present -1.
     *
     * @see blacklistPeaks
     */
    struct PeakReference
    {
      int index_in_previous_spectrum;
      int index_in_next_spectrum;
    };

    /**
     * @brief structure for peak blacklisting
     *
     * Each peak can be blacklisted. It will not pass the blacklist filter
     * unless it is of the correct charge and position in the pattern.
     */
    struct BlackListEntry
    {
      bool black;
      int black_exception_mass_shift_index;
      int black_exception_charge;
      int black_exception_mz_position;
    };

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
    MultiplexFiltering(const MSExperiment<Peak1D>& exp_picked, const std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling);

protected:
    /**
     * @brief position and blacklist filter
     *
     * Checks if there are peaks at positions corresponding to the pattern
     * and that these peaks are not blacklisted.
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param spectrum    index of the spectrum in exp_picked_ and boundaries_
     * @param peak_position    m/z positions of the peaks in spectrum
     * @param peak    index of the peak in peak_position
     * @param mz_shifts_actual    output for actual m/z shifts seen in the spectrum (will differ slightly from expected m/z shifts in pattern)
     * @param mz_shifts_actual_indices    output for indices of peaks corresponding to the pattern
     *
     * @return number of isotopic peaks seen for each peptide
     */
    int positionsAndBlacklistFilter(MultiplexPeakPattern pattern, int spectrum, std::vector<double> peak_position, int peak, std::vector<double>& mz_shifts_actual, std::vector<int>& mz_shifts_actual_indices) const;

    /**
     * @brief mono-isotopic peak intensity filter
     *
     * Quick check if the intensities of the mono-isotopic peaks are
     * above the intensity cutoff.
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param spectrum_index    index of the spectrum in exp_picked_ and boundaries_
     * @param mz_shifts_actual_indices    indices of peaks corresponding to the pattern
     *
     * @return true if all intensities above threshold
     */
    bool monoIsotopicPeakIntensityFilter(MultiplexPeakPattern pattern, int spectrum_index, const std::vector<int>& mz_shifts_actual_indices) const;

    /**
     * @brief zeroth peak filter
     *
     * The mono-isotopic peak is the first peak of each peptide. A peak one m/z shift to the left (e.g. 0.5Th for 2+)
     * is called zeroth peak. High-intensity zeroth peaks indicate incorrect pattern matches. A different pattern is
     * likely to be a better fit.
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param intensities_actual    spline-interpolated intensities at the actual m/z shift positions
     *
     * @return true if there are high-intensity zeroth peaks
     */
    bool zerothPeakFilter(MultiplexPeakPattern pattern, const std::vector<double>& intensities_actual) const;

    /**
     * @brief peptide similarity filter
     *
     * The algorithm takes only MS1 spectra into account i.e. we have no knowledge of the peptide sequences.
     * But we do know that peptides in a pair should have the same sequence and hence the same isotopic distributions.
     * The filter checks the similarity of the lightest peptide with all of the other peptides of the pattern.
     * (In high-complexity samples two peptides can have the correct mass shift by chance. Such accidental pairs
     * show different isotopic distributions and are therefore filtered out.)
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param intensities_actual    spline-interpolated intensities at the actual m/z shift positions
     * @param peaks_found_in_all_peptides_spline    number of isotopic peaks seen for each peptide (profile)
     *
     * @return true if peptide isotope patterns are similar
     */
    bool peptideSimilarityFilter(MultiplexPeakPattern pattern, const std::vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline) const;

    /**
     * @brief averagine similarity filter
     *
     * Checks similarity of the isotopic distribution with the expected averagine distribution.
     * Does the isotope distribution look like a peptide?
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param intensities_actual    spline-interpolated intensities at the actual m/z shift positions
     * @param peaks_found_in_all_peptides_spline    number of isotopic peaks seen for each peptide (profile)
     * @param mz    m/z at which the averagine distribution is calculated
     *
     * @return true if isotope distribution looks like an average peptide
     */
    bool averagineSimilarityFilter(MultiplexPeakPattern pattern, const std::vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, double mz) const;

    /**
     * @brief blacklist peaks
     *
     * If a datapoint passes all filters, the corresponding peak in this and the two neighbouring spectra is blacklisted.
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param spectrum    index of the spectrum in exp_picked_ and boundaries_
     * @param peaks_found_in_all_peptides_spline    number of isotopic peaks seen for each peptide (profile)
     */
    void blacklistPeaks(MultiplexPeakPattern pattern, int spectrum, const std::vector<int>& mz_shifts_actual_indices, int peaks_found_in_all_peptides_spline);

    /**
     * @brief returns the index of a peak at m/z
     * (finds not only a valid peak, i.e. within certain m/z deviation, but the best of the valid peaks)
     *
     * @param peak_position    m/z position of the peaks
     * @param start    index in peak_position for starting the search
     * @param mz    m/z position of the peak
     * @param scaling    rescaling of limits
     *
     * @return index of the peak in spectrum
     */
    int getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling) const;

    /**
     * @brief returns similarity of two isotope patterns
     * (simple Pearson correlation coefficient)
     *
     * @param pattern1   isotope pattern 1
     * @param pattern2   isotope pattern 2
     *
     * @return similarity (+1 best, -1 worst)
     */
    double getPatternSimilarity(std::vector<double> pattern1, std::vector<double> pattern2) const;

    /**
     * @brief returns similarity of an isotope pattern and an averagine pattern at mass m
     *
     * @param pattern   isotope pattern
     * @param m    mass at which the averagine distribution is calculated
     *
     * @return similarity (+1 best, -1 worst)
     */
    double getAveragineSimilarity(std::vector<double> pattern, double m) const;

    /**
    * @brief centroided experimental data
    */
    MSExperiment<Peak1D> exp_picked_;

    /**
    * @brief auxiliary structs for navigation and blacklisting
    */
    std::vector<std::vector<PeakReference> > registry_;
    std::vector<std::vector<BlackListEntry> > blacklist_;

    /**
     * @brief list of peak patterns
     */
    std::vector<MultiplexPeakPattern> patterns_;

    /**
     * @brief minimum number of isotopic peaks per peptide
     */
    int peaks_per_peptide_min_;

    /**
     * @brief maximum number of isotopic peaks per peptide
     */
    int peaks_per_peptide_max_;

    /**
     * @brief flag for missing peaks
     */
    bool missing_peaks_;

    /**
     * @brief intensity cutoff
     */
    double intensity_cutoff_;

    /**
     * @brief m/z shift tolerance
     */
    double mz_tolerance_;

    /**
     * @brief unit for m/z shift tolerance (ppm - true, Da - false)
     */
    bool mz_tolerance_unit_;

    /**
      * @brief peptide similarity
      */
    double peptide_similarity_;

    /**
     * @brief averagine similarity
     */
    double averagine_similarity_;

    /**
     * @brief averagine similarity scaling
     */
    double averagine_similarity_scaling_;

  };

}

#endif /* MULTIPLEXFILTERING_H_ */
