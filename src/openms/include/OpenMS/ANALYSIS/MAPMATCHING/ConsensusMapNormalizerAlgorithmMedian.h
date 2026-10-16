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

#include <map>

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
     * @param[in] map ConsensusMap
     * @param[in] method whether to use scaling or shifting to same median
     * @param[in] acc_filter string describing the regular expression for filtering accessions
     * @param[in] desc_filter string describing the regular expression for filtering descriptions
     */
    static void normalizeMaps(ConsensusMap & map, NormalizationMethod method, const std::string& acc_filter, const std::string& desc_filter);

    /**
     * @brief computes medians of all maps and returns index of map with most features
     *
     * The maps are numbered by their position among the column headers (see mapPositions_()); for maps
     * numbered 0 ... N-1 this is the map index.
     *
     * @param[in] map ConsensusMap
     * @param[out] medians vector of medians to be filled, one per map (by position)
     * @param[in] acc_filter string describing the regular expression for filtering accessions
     * @param[in] desc_filter string describing the regular expression for filtering descriptions
     * @return position of the map with the largest number of features (ColumnHeader::size)
     * @exception Exception::ElementNotFound if a feature handle's map index has no column header
     */
    static Size computeMedians(const ConsensusMap & map, std::vector<double> & medians, const std::string& acc_filter, const std::string& desc_filter);

    /**
     * @brief returns whether consensus feature passes filters
     * returns whether consensus feature @p cf_it in @p map passes accession
     * regexp @p acc_filter and description regexp @p desc_filter
     * @param[in] cf_it consensus feature
     * @param[in] map consensus map
     * @param[in] acc_filter string describing the regular expression for filtering accessions
     * @param[in] desc_filter string describing the regular expression for filtering descriptions
     */
    static bool passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const std::string& acc_filter, const std::string& desc_filter);

    /**
     * @brief Position of every map index among the column headers of @p map (in map index order)
     *
     * Map indices need not be consecutive: FileFilter keeps them when it extracts maps (-consensus:map).
     * The per-map vectors of the normalization algorithms are therefore indexed by these positions.
     * For maps numbered 0 ... N-1 the position is the map index.
     *
     * @param[in] map consensus map
     * @return map index -> position
     */
    static std::map<UInt64, Size> mapPositions_(const ConsensusMap& map);

    /**
     * @brief Position of @p map_index in @p positions (from mapPositions_())
     * @exception Exception::ElementNotFound if @p map_index has no column header
     */
    static Size mapPosition_(const std::map<UInt64, Size>& positions, UInt64 map_index);
  };

} // namespace OpenMS

