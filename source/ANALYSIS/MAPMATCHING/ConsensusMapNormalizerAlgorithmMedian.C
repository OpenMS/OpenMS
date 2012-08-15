// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

  vector<double> ConsensusMapNormalizerAlgorithmMedian::computeNormalizationFactors(const ConsensusMap& map)
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
      vector<double>& ints_j = feature_int[j];
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

  void ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(ConsensusMap& map)
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
