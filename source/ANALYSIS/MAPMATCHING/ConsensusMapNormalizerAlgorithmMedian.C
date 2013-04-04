// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace std;

namespace OpenMS
{
  ConsensusMapNormalizerAlgorithmMedian::ConsensusMapNormalizerAlgorithmMedian()
  {
  }

  ConsensusMapNormalizerAlgorithmMedian::~ConsensusMapNormalizerAlgorithmMedian()
  {

  }

  vector<double> ConsensusMapNormalizerAlgorithmMedian::computeNormalizationFactors(const ConsensusMap & map)
  {
    Size number_of_maps = map.getFileDescriptions().size();
    vector<vector<double> > feature_int(number_of_maps);
    //get map with most features, reserve space for feature_int (unequal vector lengths, 0-features omitted)
    UInt map_with_most_features = 0;
    for (UInt i = 0; i < number_of_maps; i++)
    {
      feature_int[i].reserve(map.getFileDescriptions()[i].size);
      if (map.getFileDescriptions()[i].size > map.getFileDescriptions()[map_with_most_features].size)
      {
        map_with_most_features = i;
      }
    }
    //fill feature_int with intensities
    ConsensusMap::ConstIterator cf_it;
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        feature_int[f_it->getMapIndex()].push_back(f_it->getIntensity());
      }
    }
    //compute medians factors
    vector<double> medians(number_of_maps);
    for (UInt j = 0; j < number_of_maps; j++)
    {
      vector<double> & ints_j = feature_int[j];
      gsl_sort(&ints_j.front(), 1, ints_j.size());
      medians[j] = gsl_stats_median_from_sorted_data(&ints_j.front(), 1, ints_j.size());
    }
    //compute normalization factors
    vector<double> normalization_factors(number_of_maps);
    for (UInt j = 0; j < number_of_maps; ++j)
    {
      normalization_factors[j] = medians[map_with_most_features] / medians[j];
    }

    return normalization_factors;
  }

  void ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(ConsensusMap & map)
  {
    ConsensusMap::Iterator cf_it;
    ProgressLogger progresslogger;
    progresslogger.setLogType(ProgressLogger::CMD);
    progresslogger.startProgress(0, map.size(), "normalizing maps");
    vector<double> factors = computeNormalizationFactors(map);
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      progresslogger.setProgress(cf_it - map.begin());
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        f_it->asMutable().setIntensity(f_it->getIntensity() * factors[f_it->getMapIndex()]);
      }
    }
    progresslogger.endProgress();
  }

}
