// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/StatisticFunctions.h>

using namespace std;

namespace OpenMS
{
  ConsensusMapNormalizerAlgorithmThreshold::ConsensusMapNormalizerAlgorithmThreshold() = default;

  ConsensusMapNormalizerAlgorithmThreshold::~ConsensusMapNormalizerAlgorithmThreshold() = default;

  vector<double> ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(const ConsensusMap& map, const double& ratio_threshold, const String& acc_filter, const String& desc_filter)
  {
    Size number_of_features = map.size();
    Size number_of_maps = map.getColumnHeaders().size();
    vector<vector<double> > feature_int(number_of_maps);

    //get map with most features, resize feature_int
    UInt map_with_most_features_idx = 0;
    ConsensusMap::ColumnHeaders::const_iterator map_with_most_features = map.getColumnHeaders().find(0);
    for (UInt i = 0; i < number_of_maps; i++)
    {
      feature_int[i].resize(number_of_features);
      ConsensusMap::ColumnHeaders::const_iterator it = map.getColumnHeaders().find(i);
      if (it->second.size > map_with_most_features->second.size)
      {
        map_with_most_features = it;
        map_with_most_features_idx = i;
      }
    }

    //fill feature_int with intensities
    Size pass_counter = 0;
    ConsensusMap::ConstIterator cf_it;
    UInt idx = 0;
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it, ++idx)
    {
      if (!ConsensusMapNormalizerAlgorithmMedian::passesFilters_(cf_it, map, acc_filter, desc_filter))
      {
        continue;
      }
      ++pass_counter;

      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        feature_int[f_it->getMapIndex()][idx] = f_it->getIntensity();
      }
    }

    OPENMS_LOG_INFO << endl << "Using " << pass_counter << "/" << map.size() <<  " consensus features for computing normalization coefficients" << endl << endl;

    //determine ratio
    vector<double> ratio_vector(number_of_maps);
    for (UInt j = 0; j < number_of_maps; j++)
    {
      vector<double> ratios;
      for (UInt k = 0; k < number_of_features; ++k)
      {
        if (feature_int[map_with_most_features_idx][k] != 0.0 && feature_int[j][k] != 0.0)
        {
          double ratio = feature_int[map_with_most_features_idx][k] / feature_int[j][k];
          if (ratio > ratio_threshold && ratio < 1 / ratio_threshold)
          {
            ratios.push_back(ratio);
          }
        }
      }
      if (ratios.empty())
      {
        OPENMS_LOG_WARN << endl << "Not enough features passing filters. Cannot compute normalization coefficients for all maps. Result will be unnormalized." << endl << endl;
        return vector<double>(number_of_maps, 1.0);
      }
      ratio_vector[j] = Math::mean(ratios.begin(), ratios.end());
    }
    return ratio_vector;
  }

  void ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(ConsensusMap& map, const vector<double>& ratios)
  {
    ConsensusMap::Iterator cf_it;
    ProgressLogger progresslogger;
    progresslogger.setLogType(ProgressLogger::CMD);
    progresslogger.startProgress(0, map.size(), "normalizing maps");
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      progresslogger.setProgress(cf_it - map.begin());
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        f_it->asMutable().setIntensity(f_it->getIntensity() * ratios[f_it->getMapIndex()]);
      }
    }
    progresslogger.endProgress();
  }

}
