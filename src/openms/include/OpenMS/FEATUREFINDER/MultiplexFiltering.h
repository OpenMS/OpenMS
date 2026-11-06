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
#include <map>
#include <optional>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>


namespace OpenMS
{
  class IsotopeDistribution;

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
     * spectrum back to its original spectrum: element [s][i] is the index in spectrum s
     * of the original experiment of peak i of white spectrum s.
     */
    typedef std::vector<std::vector<int> > White2Original;

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
     *
     * One spectrum per spectrum of the centroided experiment, with the blacklisted peaks in order of m/z. Call it
     * before swapping the centroided experiment out (getCentroidedExperiment()).
     */
    MSExperiment getBlacklist();

    /**
     * @brief a filtered peak that blacklisted a peak, and the mass trace of its pattern it reserved the peak for
     */
    struct BlacklistWriter
    {
      /// index of the pattern, i.e. of its filter result in the list filter() returns
      unsigned pattern_idx;
      /// index of the filtered peak in the filter result of the pattern (MultiplexFilteredMSExperiment::getPeak()); it fits
      /// in an int, since the filtered peaks of all patterns together are limited to INT_MAX (see applyBlacklistEntries_())
      unsigned peak_idx;
      /// mass trace of the pattern the peak was reserved for (peptide * isotopes_per_peptide_max + isotope)
      int mass_trace;
    };

    /**
     * @brief who blacklisted a peak: the first and the last filtered peak that reserved it (the same if only one did)
     */
    struct BlacklistOwner
    {
      /// the first filtered peak that reserved the peak
      BlacklistWriter first;
      /// the last one, whose mass trace the blacklist holds
      BlacklistWriter last;
    };

    /**
     * @brief the owners of the blacklisted peaks, aligned with getBlacklist(): one list per spectrum, one entry per peak
     *
     * A pattern reserves the peaks of its mass traces over their RT range extended by the RT band, also peaks that
     * another filtered peak reserved before, so a peak can be reserved several times. The owners are recorded when a
     * peak passes, in the order in which peaks pass, so they do not depend on the number of threads.
     *
     * Only the first and the last writer of a peak are kept. A filtered peak that reserved it in between is lost: if its
     * multiplet is reported and the clusters of the first and the last writer are dropped, the peak counts as reserved by
     * no reported feature (FeatureFinderMultiplexAlgorithm::getBlacklist()). That takes at least three writers of one peak,
     * which in practice requires an RT band > 0 (@em rt_band_, FeatureFinderMultiplex's algorithm:rt_band).
     *
     * This copy holds the owners of all blacklisted peaks at once; getBlacklistOwners(Size) returns those of one spectrum.
     */
    std::vector<std::vector<BlacklistOwner>> getBlacklistOwners() const;

    /**
     * @brief the owners of the blacklisted peaks of one spectrum, aligned with that spectrum of getBlacklist()
     *
     * The same entries as getBlacklistOwners()[@p spectrum_index], without a copy of the owners of the other spectra.
     *
     * @param[in] spectrum_index    index of the spectrum in the centroided experiment (and in getBlacklist())
     *
     * @throws Exception::IndexOverflow if @p spectrum_index is not below the number of spectra
     */
    std::vector<BlacklistOwner> getBlacklistOwners(Size spectrum_index) const;

    /**
     * @brief take the lower isotopologue of deuterium labels into account where a peak one isotope below a peptide is checked
     *
     * A label with n deuterium (2H) atoms of atom fraction p < 1 leaves molecules that lack one 2H. Their monoisotopic peak
     * lies 1.0063/z below the peptide's, 0.0029/z from its M-1, so within the m/z tolerance in which the filter looks for a
     * significant peak one isotope below the peptide (which would make the peptide's peaks the isotopes of another envelope).
     * For one Lys4 at p = 0.96 it is about 15 % of the monoisotopic peak, for two about 28 % (at 2 kDa), against the 30 % above
     * which such a peak is significant. With this setting, a peak there is significant only above 30 % plus that expected part
     * (deuteriumLowerFraction() of the averagine distribution at the peptide's mass). Peptides whose labels have no 2H keep 30 %.
     *
     * Call it after the construction, before filter().
     *
     * @param[in] deuterium_per_label    number of 2H atoms per label of the patterns' label sets (short names, e.g. Lys4 -> 4,
     *                                   see MultiplexDeltaMassesGenerator::getDeuteriumCount()); labels that are not listed have none
     * @param[in] enrichment    2H atom fraction of the labels, in (0, 1]; 1 (the default) ignores the lower isotopologue
     *
     * @throw Exception::InvalidValue if @p enrichment is not in (0, 1]
     */
    void setDeuteriumLabels(const std::map<std::string, unsigned>& deuterium_per_label, double enrichment);

    /**
     * @brief the peak one isotope below a deuterium-labelled peptide that its own label isotopologues produce, relative to its monoisotopic peak
     *
     * The molecules that lack j of the @p deuterium 2H atoms of the label (binomial in the 2H atom fraction @p enrichment,
     * q_j = C(n, j) ((1 - p) / p)^j relative to the fully labelled ones) have the peptide's isotope distribution I shifted
     * by -j * 1.0063 Da: their isotope j - 1 coincides with the peptide's M-1 within 0.0029 Da per 2H, their isotope j
     * with its M+0. The result is sum_{j>=1} q_j I_{j-1} / sum_{j>=0} q_j I_j, about r / (1 + r I_1/I_0) for r = n (1 - p) / p.
     *
     * @param[in] deuterium    number of 2H atoms of the peptide's labels
     * @param[in] enrichment    2H atom fraction of the labels, in (0, 1]
     * @param[in] isotopes    the peptide's isotope distribution without the labels' 2H (e.g. averagine)
     *
     * @return 0 if @p deuterium is 0, @p enrichment is not below 1 or @p isotopes has fewer than two isotopes
     */
    static double deuteriumLowerFraction(unsigned deuterium, double enrichment, const IsotopeDistribution& isotopes);

    /**
     * @brief the 2H enrichment at which deuteriumLowerFraction() gives @p fraction (its inverse in the enrichment)
     *
     * The fraction decreases monotonically with the enrichment, from its value at @p min_enrichment to 0 at 1. The result is
     * found by bisection and clamped to [@p min_enrichment, 1]: 1 for a fraction of at most 0, @p min_enrichment for a
     * fraction at least as large as at @p min_enrichment.
     *
     * @param[in] deuterium    number of 2H atoms of the peptide's labels (> 0)
     * @param[in] fraction    the peak one isotope below the peptide relative to its monoisotopic peak, e.g. measured
     * @param[in] isotopes    the peptide's isotope distribution without the labels' 2H (at least two isotopes)
     * @param[in] min_enrichment    lower bound of the enrichment, in (0, 1)
     *
     * @return the enrichment, or NaN if @p deuterium is 0 or @p isotopes has fewer than two isotopes
     */
    static double deuteriumEnrichmentFromLowerFraction(unsigned deuterium, double fraction, const IsotopeDistribution& isotopes,
                                                       double min_enrichment = 0.5);

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
     * @brief the owners of a peak of @em blacklist_ that more than one filtered peak reserved (indices in @em filtered_peaks_)
     */
    struct BlacklistOwnerRecord
    {
      /// the first filtered peak that reserved the peak
      int first;
      /// the value @em first left in @em blacklist_ (read when another filtered peak reserved the peak for the first time)
      int first_value;
      /// the last filtered peak that reserved the peak, whose value @em blacklist_ holds
      int last;
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
     * @param[in] fraction    a peak is significant if its intensity exceeds this fraction of @p intensity_first_peak
     *
     * @return -1 (if there is no significant peak), or the index of the significant peak in the spectrum @p it_rt points to
     * (in filterPeakPositions_() a spectrum of @em exp_centroided_white_, see @em exp_centroided_mapping_)
     */
    int checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak,
                                 double fraction = 0.3) const;

    /**
     * @brief the part of the peak one isotope below peptide @p peptide of @p pattern (monoisotopic m/z @p mz) that its own
     * deuterium labels produce, relative to its monoisotopic peak (see setDeuteriumLabels()); 0 without 2H
     */
    double deuteriumLowerFraction_(const MultiplexIsotopicPeakPattern& pattern, size_t peptide, double mz) const;

    /**
     * @brief check if there are enough peaks in the RT band to form the pattern
     *
     * Checks if there are peaks at m/z positions corresponding to the pattern
     * and that the primary peak position is not blacklisted.
     *
     * @param[in] mz    m/z of the primary peak
     * @param[in] it_rt_begin    RT iterator of the very first spectrum of @em exp_centroided_white_ (needed to determine indices;
     * peak indices in its spectra are mapped to @em exp_centroided_ with @em exp_centroided_mapping_)
     * @param[in] it_rt_band_begin    RT iterator of the first spectrum in the RT band (of @em exp_centroided_white_)
     * @param[in] it_rt_band_end    RT iterator of the spectrum after the last spectrum in the RT band (of @em exp_centroided_white_)
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
     * Each of the satellites is associated with a specific mass trace. We blacklist
     * all peaks in these mass traces (even if they are not a satellite) extending them
     * by a margin @em rt_band_.
     *
     * @param[in] peak    peak to be blacklisted
     * @param[in] pattern_idx    index of the pattern in @em patterns_
     * @param[in] peak_idx    index of @p peak in the filter result of the pattern, recorded as owner of the peaks
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak, unsigned pattern_idx, size_t peak_idx);

    /**
     * @brief write @p entries to @em blacklist_, with the filtered peak (@p pattern_idx, @p peak_idx) as their owner
     *
     * It becomes the last owner of each peak. As long as a single filtered peak has reserved a peak, @em blacklist_owner_ holds
     * that owner. When a second filtered peak reserves it, the peak gets a record in @em blacklist_owner_records_ with the first
     * owner, the value the first owner left in @em blacklist_, and the last owner. Later owners only replace the last owner of
     * the record, so any owner in between is lost (see getBlacklistOwners()).
     *
     * Owners and records are indexed by int, so the filtered peaks of all patterns together, and the peaks reserved by more
     * than one of them, are limited to INT_MAX each.
     *
     * @throws Exception::InvalidSize if this filtered peak, or the records its entries may add, would exceed that limit
     */
    void applyBlacklistEntries_(const std::vector<BlacklistEntry>& entries, unsigned pattern_idx, size_t peak_idx);

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
     * @brief averagine isotope distribution (of the type @em averagine_type_, with @em isotopes_per_peptide_max_ isotopes)
     * of a molecule of @p mass
     *
     * The returned reference is valid until the next call on the same thread. Peptide distributions are computed once
     * per element composition of the averagine model.
     *
     * @throw Exception::InvalidParameter for an unknown averagine type
     */
    const IsotopeDistribution& averagineDistribution_(double mass) const;
    
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
     * @brief for every peak of @em blacklist_, who reserved it
     *
     * -1: nobody (the peak is white); >= 0: the only filtered peak that reserved it so far (index in @em filtered_peaks_);
     * <= -2: more than one filtered peak reserved it, and -2 - value is the index of its record in @em blacklist_owner_records_.
     */
    std::vector<std::vector<int>> blacklist_owner_;

    /**
     * @brief the first and the last owner of every peak reserved by more than one filtered peak (see @em blacklist_owner_)
     */
    std::vector<BlacklistOwnerRecord> blacklist_owner_records_;

    /**
     * @brief every filtered peak that passed, in the order they passed: index of its pattern and its index in the pattern's filter result
     */
    std::vector<std::pair<unsigned, unsigned>> filtered_peaks_;
    
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
     * @brief the decision of filterSpectrum_() on each peak of the spectrum (reused, as blacklist_reads_)
     */
    std::vector<std::optional<MultiplexFilteredPeak>> passed_;

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

    /**
     * @brief number of 2H atoms of the labels of each peptide of each mass shift pattern (index: MultiplexIsotopicPeakPattern::getMassShiftIndex()),
     * empty for a pattern without 2H (setDeuteriumLabels())
     */
    std::vector<std::vector<unsigned>> deuterium_;

    /**
     * @brief 2H atom fraction of the deuterium labels (setDeuteriumLabels()); 1: their lower isotopologue is ignored
     */
    double deuterium_enrichment_ = 1.0;

  };

}

