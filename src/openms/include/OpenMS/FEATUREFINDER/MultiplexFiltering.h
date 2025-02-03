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
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/serialization/strong_typedef.hpp>

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
     * @param exp_centroided    experimental data in centroid mode
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
     * @param averagine_type Averagine model to use: 'peptide', 'RNA', 'DNA'
     */
    MultiplexFiltering(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min,
                       int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit,
                       double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type="peptide");
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
     * @param mz    position where the peak is expected
     * @param mz_tolerance    m/z tolerance within the peak may lie
     * @param it_rt    pointer to the spectrum
     * @param intensity_first_peak    intensity to compare to
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
     * @param mz    m/z of the primary peak
     * @param it_rt_begin    RT iterator of the very first spectrum of the experiment (needed to determine indices)
     * @param it_rt_band_begin    RT iterator of the first spectrum in the RT band
     * @param it_rt_band_end    RT iterator of the spectrum after the last spectrum in the RT band
     * @param pattern    m/z pattern to search for
     * @param peak    filter result output
     *
     * @return boolean if this filter was passed i.e. there are @em isotopes_per_peptide_min_ or more mass traces which form the pattern.
     */
    bool filterPeakPositions_(double mz, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak) const;

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
     * by a margin @em rt_band_.
     * 
     * @param peak    peak to be blacklisted
     * @param pattern_idx    index of the pattern in @em patterns_
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak, unsigned pattern_idx);
    
    /**
     * @brief check if the satellite peaks conform with the averagine model
     *
     * Check if the intensities of the satellite peaks correlate with the peak intensities
     * of the averagine model. We check both Pearson and Spearman rank correlation.
     *
     * @param pattern    m/z pattern to search for
     * @param peak    peak with set of satellite peaks
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
     * @param pattern    m/z pattern to search for
     * @param peak    peak with set of satellite peaks
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

    String averagine_type_;

  };

}

