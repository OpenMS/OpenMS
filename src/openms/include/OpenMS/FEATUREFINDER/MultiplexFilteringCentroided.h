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
   * @see MultiplexIsotopicPeakPattern
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
     * @param[in] exp_centroided    experimental data in centroid mode
     * @param[in] patterns    patterns of isotopic peaks to be searched for
     * @param[in] isotopes_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param[in] isotopes_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param[in] intensity_cutoff    intensity cutoff
     * @param[in] rt_band    RT range for filtering
     * @param[in] mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param[in] mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param[in] peptide_similarity    similarity score for two peptides in the same multiplet
     * @param[in] averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param[in] averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p).
     * @param[in] averagine_type    The averagine model to use, current options are RNA DNA or peptide.
     */
    MultiplexFilteringCentroided(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type="peptide");

    /**
     * @brief filter for patterns
     * (generates a filter result for each of the patterns)
     *
     * @see MultiplexIsotopicPeakPattern, MultiplexFilterResult
     */
    std::vector<MultiplexFilteredMSExperiment> filter();

  };

}

