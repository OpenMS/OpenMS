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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <functional>
#include <optional>
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
   * @see MultiplexFilteredMSExperiment
   * @see MultiplexFilteringCentroided
   * @see MultiplexFilteringProfile
   */
  class OPENMS_DLLAPI MultiplexFiltering :
    public ProgressLogger
  {
public:
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
     * @param[in] exp_centroided    experimental data in centroid mode
     * @param[in] patterns    patterns of isotopic peaks to be searched for
     * @param[in] isotopes_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param[in] isotopes_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param[in] intensity_cutoff    intensity cutoff
     * @param[in] rt_band    RT range used for filtering
     * @param[in] mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param[in] mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param[in] peptide_similarity    similarity score for two peptides in the same multiplet
     * @param[in] averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param[in] averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p).
     * @param[in] averagine_type Averagine model to use: 'peptide', 'RNA', 'DNA'
     */
    MultiplexFiltering(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min,
                       int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit,
                       double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, std::string averagine_type="peptide");
    /**
     * @brief returns the intensity-filtered, centroided spectral data
     */
    MSExperiment& getCentroidedExperiment();
    
    /**
     * @brief returns the blacklisted, centroided peaks
     */
    MSExperiment getBlacklist();

protected:
    /**
     * @brief an entry of @em blacklist_ (spectrum index and peak index in @em exp_centroided_) and a value read from or written to it
     */
    struct BlacklistEntry
    {
      size_t rt_idx;
      size_t mz_idx;
      int value;
    };

    /**
     * @brief construct an MS experiment from exp_centroided_ containing
     * peaks which have not been previously blacklisted in blacklist_
     * 
     * In addition, construct an index mapping of 'white' peak positions
     * to their position in the corresponding, original spectrum.
     */
    void updateWhiteMSExperiment_();

    /**
     * @brief check for significant peak
     *
     * @param[in] mz    position where the peak is expected
     * @param[in] mz_tolerance    m/z tolerance within the peak may lie
     * @param[in] it_rt    pointer to the spectrum
     * @param[in] intensity_first_peak    intensity to compare to
     *
     * @return -1 (if there is no significant peak), or peak index mz_idx (if there is a significant peak)
     */
    int checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak) const;

    /**
     * @brief check if there are enough peaks in the RT band to form the pattern
     *
     * Checks if there are peaks at m/z positions corresponding to the pattern
     * and that the primary peak position is not blacklisted.
     *
     * @param[in] mz    m/z of the primary peak
     * @param[in] it_rt_begin    RT iterator of the very first spectrum of the experiment (needed to determine indices)
     * @param[in] it_rt_band_begin    RT iterator of the first spectrum in the RT band
     * @param[in] it_rt_band_end    RT iterator of the spectrum after the last spectrum in the RT band
     * @param[in] pattern    m/z pattern to search for
     * @param[out] peak    filter result output
     * @param[out] blacklist_reads    if not null, every entry of @em blacklist_ read is appended to it, with the value read
     *
     * @return boolean if this filter was passed i.e. there are @em isotopes_per_peptide_min_ or more mass traces which form the pattern.
     */
    bool filterPeakPositions_(double mz, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak, std::vector<BlacklistEntry>* blacklist_reads = nullptr) const;

    /**
     * @brief blacklist this peak
     *
     * Blacklist all satellites associated with this peak.
     *
     * @param[in] peak    peak to be blacklisted
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak);

    /**
     * @brief blacklist this peak
     *
     * Each of the satellites is associated with a specific mass trace. We blacklist
     * all peaks in these mass traces (even if they are not a satellite) extending them
     * by a margin @em rt_band_.
     *
     * @param[in] peak    peak to be blacklisted
     * @param[in] pattern_idx    index of the pattern in @em patterns_
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak, unsigned pattern_idx);

    /**
     * @brief the entries of @em blacklist_ that blacklistPeak_() sets, with their new values, in the order it sets them
     *
     * They do not depend on the blacklist.
     *
     * @param[in] peak    peak to be blacklisted
     * @param[in] pattern_idx    index of the pattern in @em patterns_
     * @param[out] entries    the entries are appended to it
     */
    void blacklistEntries_(const MultiplexFilteredPeak& peak, unsigned pattern_idx, std::vector<BlacklistEntry>& entries) const;

    /**
     * @brief decides on one peak of a 'white' spectrum
     *
     * Called with the peak index in the white spectrum and, if not null, the list to which every entry of
     * @em blacklist_ read is to be appended (see filterPeakPositions_()). Returns the peak with its satellites
     * if it passes all filters.
     */
    typedef std::function<std::optional<MultiplexFilteredPeak>(SignedSize, std::vector<BlacklistEntry>*)> PeakFilter;

    /**
     * @brief filter all peaks of one 'white' spectrum, add the ones that pass to @p result and blacklist them
     *
     * Whether a peak passes depends on the blacklist, which every peak that passes extends. The result is that of
     * filtering the peaks one after another in order of m/z (as a single thread does), whatever the number of threads:
     * The peaks are filtered in parallel against the blacklist as it is when the spectrum is reached; every decision
     * records the blacklist entries it read, and every peak that passes the entries it would blacklist. Then the peaks
     * are accepted in order of m/z. A decision stands unless a peak accepted before it in this spectrum has changed one
     * of the entries it read; then the peak is filtered again against the current blacklist. With one thread, the
     * peaks are simply filtered in order of m/z.
     *
     * @param[in] idx_rt    index of the spectrum in @em exp_centroided_white_
     * @param[in] pattern_idx    index of the pattern in @em patterns_
     * @param[in] filter_peak    decides on one peak of the spectrum
     * @param[out] result    filter result the peaks that pass are added to
     */
    void filterSpectrum_(size_t idx_rt, unsigned pattern_idx, const PeakFilter& filter_peak, MultiplexFilteredMSExperiment& result);

    /**
     * @brief check if the satellite peaks conform with the averagine model
     *
     * Check if the intensities of the satellite peaks correlate with the peak intensities
     * of the averagine model. We check both Pearson and Spearman rank correlation.
     *
     * @param[in] pattern    m/z pattern to search for
     * @param[in] peak    peak with set of satellite peaks
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than @em averagine_similarity_
     */
    bool filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const;
    
    /**
     * @brief check if corresponding satellite peaks of different peptides show a good correlation
     *
     * Different peptides in the same multiplet have the same amino acid sequence and should therefore exhibit very similar
     * isotope distributions. The filter checks if satellite peaks corresponding to different isotopes in different peptide
     * features show a strong correlation. The filter is of course ignored for singlet feature detection.
     *
     * @param[in] pattern    m/z pattern to search for
     * @param[in] peak    peak with set of satellite peaks
     *
     * @return boolean if this filter was passed i.e. the correlation coefficient is greater than @em peptide_similarity_
     */
    bool filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const;

    /**
     * @brief centroided experimental data
     */
    MSExperiment exp_centroided_;

    /**
     * @brief auxiliary structs for blacklisting
     */
    std::vector<std::vector<int> > blacklist_;
    
    /**
     * @brief for each peak of the spectrum being filtered, the entries of @em blacklist_ its decision read (see filterSpectrum_())
     *
     * Kept from spectrum to spectrum, so that the lists keep their memory.
     */
    std::vector<std::vector<BlacklistEntry>> blacklist_reads_;

    /**
     * @brief for each peak of the spectrum being filtered that passes, the entries of @em blacklist_ it sets (see filterSpectrum_())
     */
    std::vector<std::vector<BlacklistEntry>> blacklist_writes_;

    /**
     * @brief "white" centroided experimental data
     *
     * subset of all peaks of @em exp_centroided_ which are not blacklisted in @em blacklist_
     */
    MSExperiment exp_centroided_white_;
    
    /**
     * @brief mapping of peak indices from a 'white' experiment @em exp_centroided_white_ to its original experiment @em exp_centroided_
     */
    White2Original exp_centroided_mapping_;

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
     * @brief intensity cutoff
     */
    double intensity_cutoff_;

    /**
     * @brief RT range used for filtering
     */
    double rt_band_;
    
    /**
     * @brief m/z shift tolerance
     */
    double mz_tolerance_;

    /**
     * @brief unit for m/z shift tolerance (ppm - true, Da - false)
     */
    bool mz_tolerance_unit_in_ppm_;

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

    std::string averagine_type_;

  };

}

