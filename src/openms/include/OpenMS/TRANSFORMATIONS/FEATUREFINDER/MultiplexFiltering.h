// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
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
   * @see MultiplexIsotopicPeakPattern
   * @see MultiplexFilterResult
   * @see MultiplexFilteringCentroided
   * @see MultiplexFilteringProfile
   */
  class OPENMS_DLLAPI MultiplexFiltering :
    public ProgressLogger
  {
public:
    /**
     * @brief type for peak blacklisting
     * 
     * white    white in this and subsequent patterns
     * grey     white in this pattern and black in subsequent patterns
     * black    black in this and in subsequent patterns
     * 
     * We assume that one peak cannot belong to two or more patterns
     * i.e. peptides at the same time.
     */
    enum BlacklistEntry
    {
      white,
      grey,
      black
    };
    
    /**
     * @brief index mapping from a 'white' experiment to its original experiment
     * 
     * An MSExperiment contains a set of spectra each containing a number of peaks.
     * In the course of the filtering, some peaks are blacklisted since they are
     * identified to belong to a certain pattern i.e. peptide. An experiment in
     * which blacklisted peaks are removed is called 'white'. White spectra
     * contain fewer peaks than their corresponding primary spectra. Consequently,
     * their indices are shifted. The type maps a peak index in a 'white'
     * spectrum back to its original spectrum.
     */
    typedef std::vector<std::map<int, int> > White2Original;

    /**
     * @brief constructor
     *
     * @param exp_picked    experimental data in centroid mode
     * @param patterns    patterns of isotopic peaks to be searched for
     * @param isotopes_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param isotopes_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param missing_peaks    flag for missing peaks
     * @param intensity_cutoff    intensity cutoff
     * @param rt_band    RT range used for filtering
     * @param rt_band_fraction
     * @param mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param peptide_similarity    similarity score for two peptides in the same multiplet
     * @param averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p). 
     */
    MultiplexFiltering(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, bool missing_peaks, double intensity_cutoff, double rt_band, double rt_band_fraction, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type="peptide");

protected:
    /**
     * @brief construct an MS experiment from exp_picked_ containing
     * peaks which have not previously blacklisted in blacklist_
     * 
     * @param mapping    index mapping of 'white' peak positions to their position in the corresponding, original spectrum 
     */
    MSExperiment getWhiteMSExperiment_(White2Original& mapping); 

    /**
     * @brief check for significant peak
     *
     * @param mz    position where the peak is expected
     * @param mz_tolerance    m/z tolerance within the peak may lie
     * @param it_rt    pointer to the spectrum
     * @param intensity_first_peak    intensity to compare to
     *
     * @return boolean if there is a significant peak
     */
    bool checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak) const;

    /**
     * @brief check if there are enough peaks in the RT band to form the pattern
     *
     * Checks if there are peaks at m/z positions corresponding to the pattern
     * and that the primary peak position is not blacklisted.
     *
     * @param it_mz    m/z iterator of the primary
     * @param it_rt_begin    RT iterator of the very first spectrum of the experiment (needed to determine indices)
     * @param it_rt_band_begin    RT iterator of the first spectrum in the RT band
     * @param it_rt_band_end    RT iterator of the last spectrum in the RT band
     * @param pattern    m/z pattern to search for
     * @param peak    filter result output
     *
     * @return boolean if this filter was passed i.e. there are <isotopes_per_peptide_min_> or more mass traces which form the pattern.
     */
    bool filterPeakPositions_(const MSSpectrum<Peak1D>::ConstIterator& it_mz, White2Original& index_mapping, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak) const;

    /**
     * @brief blacklist this peak
     * 
     * Blacklist all satellites associated with this peak.
     * 
     * @param peak    peak to be blacklisted
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak);
    
    /**
     * @brief blacklist this peak
     * 
     * Each of the satellites is associated with a specific mass trace. We blacklist
     * all peaks in these mass traces (even if they are not a satellite) extending them
     * by a margin <rt_band_>.  
     * 
     * @param peak    peak to be blacklisted
     * @param pattern_idx    index of the pattern in <patterns_>
     */
    void blacklistPeak2_(const MultiplexFilteredPeak& peak, unsigned pattern_idx);

    /**
     * @brief turn grey blacklist_ entries into black ones
     * 
     * Grey entries function as white in the current pattern but black in subsequent patterns,
     * i.e. at the end of a pattern these entries need to be turned black.
     */
    void ungreyBlacklist_();

    /**
     * @brief check if the satellite peaks conform with the averagine model
     *
     * Check if the intensities of the satellite peaks correlate with the peak intensities
     * of the averagine model. We check both Pearson and Spearman rank correlation.
     *
     * @param pattern    m/z pattern to search for
     * @param peak    peak with set of satellite peaks
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than <averagine_similarity_>
     */
    bool filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const;
    
    /**
     * @brief check if corresponding satellite peaks of different peptides show a good correlation
     *
     * Different peptides in the same multiplet have the same amino acid sequence and should therefore exhibit very similar
     * isotope distributions. The filter checks if satellite peaks corresponding to different isotopes in different peptide
     * features show a strong correlation. The filter is of course ignored for singlet feature detection.
     *
     * @param pattern    m/z pattern to search for
     * @param peak    peak with set of satellite peaks
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than <peptide_similarity_>
     */
    bool filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    bool monoIsotopicPeakIntensityFilter_(const MultiplexIsotopicPeakPattern& pattern, int spectrum_index, const std::vector<int>& mz_shifts_actual_indices) const;

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
    bool zerothPeakFilter_(const MultiplexIsotopicPeakPattern& pattern, const std::vector<double>& intensities_actual) const;

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
    bool peptideSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const std::vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline) const;

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
    bool averagineSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const std::vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, double mz) const;

    /**
     * @brief blacklist peaks
     *
     * If a datapoint passes all filters, the corresponding peak in this and the two neighbouring spectra is blacklisted.
     *
     * @param pattern    pattern of isotopic peaks to be searched for
     * @param spectrum    index of the spectrum in exp_picked_ and boundaries_
     * @param peaks_found_in_all_peptides_spline    number of isotopic peaks seen for each peptide (profile)
     */
    //void blacklistPeaks_(const MultiplexIsotopicPeakPattern& pattern, int spectrum, const std::vector<int>& mz_shifts_actual_indices, int peaks_found_in_all_peptides_spline);

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
    int getPeakIndex_(const std::vector<double>& peak_position, int start, double mz, double scaling) const;

    /**
     * @brief returns similarity of two isotope patterns
     * (simple Pearson correlation coefficient)
     *
     * @param pattern1   isotope pattern 1
     * @param pattern2   isotope pattern 2
     *
     * @return similarity (+1 best, -1 worst)
     */
    double getPatternSimilarity_(const std::vector<double>& pattern1, const std::vector<double>& pattern2) const;

    /**
     * @brief returns similarity of an isotope pattern and an averagine pattern at mass m
     *
     * @param pattern   isotope pattern
     * @param m    mass at which the averagine distribution is calculated
     *
     * @return similarity (+1 best, -1 worst)
     */

    double getAveragineSimilarity_(const std::vector<double>& pattern, double m) const;

    /**
    * @brief centroided experimental data
    */
    MSExperiment exp_picked_;

    /**
    * @brief auxiliary structs for blacklisting
    */
    std::vector<std::vector<BlacklistEntry> > blacklist_;

    /**
     * @brief list of peak patterns
     */
    std::vector<MultiplexIsotopicPeakPattern> patterns_;

    /**
     * @brief minimum number of isotopic peaks per peptide
     */
    size_t isotopes_per_peptide_min_;

    /**
     * @brief maximum number of isotopic peaks per peptide
     */
    size_t isotopes_per_peptide_max_;

    /**
     * @brief flag for missing peaks
     */
    bool missing_peaks_;

    /**
     * @brief intensity cutoff
     */
    double intensity_cutoff_;

    /**
     * @brief RT range used for filtering
     */
    double rt_band_;
    
    /**
     * @brief
     */
    double rt_band_fraction_;
    
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

    /**
     * @brief type of averagine to use
     */
    String averagine_type_;

  };

}

#endif /* MULTIPLEXFILTERING_H_ */
