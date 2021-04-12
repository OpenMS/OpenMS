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

