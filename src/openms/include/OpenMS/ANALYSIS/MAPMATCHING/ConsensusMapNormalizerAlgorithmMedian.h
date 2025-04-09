// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>

namespace OpenMS
{

  /**
   * @brief Algorithms of ConsensusMapNormalizer
   *
   */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmMedian
  {
    /// our friend can use our passesWhitelist_() method in order to avoid code duplication (without overdesigning this)
    friend class ConsensusMapNormalizerAlgorithmThreshold;

private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian(const ConsensusMapNormalizerAlgorithmMedian & copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian & operator=(const ConsensusMapNormalizerAlgorithmMedian & rhs);

public:
    /// default constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian();

    /// destructor is not implemented -> private
    virtual ~ConsensusMapNormalizerAlgorithmMedian();

    /**
     * @brief The NormalizationMethod enum
     * Whether to scale to same median using division/multiplication or shift using subtraction/addition
     */
    enum NormalizationMethod { NM_SCALE, NM_SHIFT };

    /**
     * @brief normalizes the maps of the consensusMap
     * @param map ConsensusMap
     * @param method whether to use scaling or shifting to same median
     * @param acc_filter string describing the regular expression for filtering accessions
     * @param desc_filter string describing the regular expression for filtering descriptions
     */
    static void normalizeMaps(ConsensusMap & map, NormalizationMethod method, const String& acc_filter, const String& desc_filter);

    /**
     * @brief computes medians of all maps and returns index of map with most features
     * @param map ConsensusMap
     * @param medians vector of medians to be filled
     * @param acc_filter string describing the regular expression for filtering accessions
     * @param desc_filter string describing the regular expression for filtering descriptions
     * @return index of map with largest number of features
     */
    static Size computeMedians(const ConsensusMap & map, std::vector<double> & medians, const String& acc_filter, const String& desc_filter);

    /**
     * @brief returns whether consensus feature passes filters
     * returns whether consensus feature @p cf_it in @p map passes accession
     * regexp @p acc_filter and description regexp @p desc_filter
     * @param cf_it consensus feature
     * @param map consensus map
     * @param acc_filter string describing the regular expression for filtering accessions
     * @param desc_filter string describing the regular expression for filtering descriptions
     */
    static bool passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const String& acc_filter, const String& desc_filter);
  };

} // namespace OpenMS

