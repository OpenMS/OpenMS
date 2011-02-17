// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithm.h>
#include <gsl/gsl_statistics.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace std;

namespace OpenMS
{
	vector<double> ConsensusMapNormalizerAlgorithm::computeCorrelation(const ConsensusMap& map, const double& ratio_threshold)
	{
		UInt number_of_features = map.size();
		UInt number_of_maps = map.getFileDescriptions().size();
		vector<vector<double> > feature_int(number_of_maps);
		//get map with most features, resize feature_int
		UInt map_with_most_features = 0;
		for (UInt i = 0; i < number_of_maps; i++)
		{
			feature_int[i].resize(number_of_features);
			if (map.getFileDescriptions()[i].size > map.getFileDescriptions()[map_with_most_features].size)
			{
				map_with_most_features = i;
			}
		}
		//fill feature_int with intensities
		ConsensusMap::ConstIterator cf_it;
		UInt idx = 0;
		for (cf_it = map.begin(); cf_it != map.end(); ++cf_it, ++idx)
		{
			ConsensusFeature::HandleSetType::const_iterator f_it;
			for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
			{
				feature_int[f_it->getMapIndex()][idx] = f_it->getIntensity();
			}
		}
		//determine ratio
		vector<double> ratio_vector(number_of_maps);
		for (UInt j = 0; j < number_of_maps; j++)
		{
			vector<double> ratios;
			for (UInt k = 0; k < number_of_features; ++k)
			{
				if (feature_int[map_with_most_features][k] != 0.0 && feature_int[j][k] != 0.0)
				{	
					double ratio = feature_int[map_with_most_features][k] / feature_int[j][k];
					if (ratio > ratio_threshold && ratio < 1/ratio_threshold)
					{
						ratios.push_back(ratio);
					}
				}
			}
			ratio_vector[j] = gsl_stats_mean(&ratios.front(), 1, ratios.size());
		}
		return ratio_vector;
	}

	void ConsensusMapNormalizerAlgorithm::normalizeMaps(ConsensusMap& map, const vector<double>& ratios)
	{
		ConsensusMap::Iterator cf_it;
		ProgressLogger progresslogger;
		progresslogger.setLogType(ProgressLogger::CMD);
		progresslogger.startProgress(0, map.size(), "normalizing maps");
		for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
		{
			progresslogger.setProgress(cf_it - map.begin());
			ConsensusFeature::HandleSetType::iterator f_it;
			for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
			{	
				f_it->asMutable().setIntensity(f_it->getIntensity() * ratios[f_it->getMapIndex()]);
			}
		}
		progresslogger.endProgress();
	}

} 
