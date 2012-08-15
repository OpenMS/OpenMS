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

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <cmath>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace std;

namespace OpenMS
{
  ConsensusMapNormalizerAlgorithmQuantile::ConsensusMapNormalizerAlgorithmQuantile()
  {
  }

  ConsensusMapNormalizerAlgorithmQuantile::~ConsensusMapNormalizerAlgorithmQuantile()
  {

  }

  void ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(ConsensusMap& map)
	{
    //extract feature intensities
    vector<vector<double> > feature_ints;
    extractIntensityVectors(map, feature_ints);
    Size number_of_maps = feature_ints.size();

    //determine largest number of features in any map
    UInt largest_number_of_features = 0;
    for (UInt i = 0; i < number_of_maps; ++i)
    {
      if (feature_ints[i].size() > largest_number_of_features)
      {
        largest_number_of_features = feature_ints[i].size();
      }
    }

    //resample n data points from each sorted intensity distribution (from the different maps), n = maximum number of features in any map
    vector<vector<double> > resampled_sorted_data;
    for (UInt i = 0; i < number_of_maps; ++i)
    {
      vector<double> sorted = feature_ints[i];
      gsl_sort(&sorted.front(), 1, sorted.size());
      vector<double> resampled(largest_number_of_features);
      resample(sorted, resampled, largest_number_of_features);
      resampled_sorted_data.push_back(resampled);
    }

    //compute reference distribution from all resampled distributions
    vector<double> reference_distribution(largest_number_of_features);
    for (UInt i = 0; i < number_of_maps; ++i)
    {
      for (UInt j = 0; j < largest_number_of_features; ++j)
      {
        reference_distribution[j] += (resampled_sorted_data[i][j] / (double)number_of_maps);
      }
    }

    //for each map: resample from the reference distribution down to the respective original size again
    vector<vector<double> > normalized_sorted_ints(number_of_maps);
    for (UInt i = 0; i < number_of_maps; ++i)
    {
      vector<double> ints;
      resample(reference_distribution, ints, feature_ints[i].size());
      normalized_sorted_ints[i] = ints;
    }

    //set the intensities of feature_ints to the normalized intensities
    for (UInt i = 0; i < number_of_maps; ++i)
    {
      vector<Size> sort_indices(feature_ints[i].size());
      gsl_sort_index(&sort_indices.front(), &feature_ints[i].front(), 1, feature_ints[i].size());
      Size k = 0;
      for (Size j = 0; j < sort_indices.size(); ++j)
      {
        Size idx = sort_indices[j];
        feature_ints[i][idx] = normalized_sorted_ints[i][k++];
      }
    }

    //write new feature intensities to the consensus map
    setNormalizedIntensityValues(feature_ints, map);
	}

  void ConsensusMapNormalizerAlgorithmQuantile::resample(const vector<double>& data_in, vector<double>& data_out, UInt n_resampling_points)
  {
    data_out.clear();
    data_out.resize(n_resampling_points);

    if (n_resampling_points == 0)
    {
      return;
    }

    data_out[0] = data_in.front();
    data_out[n_resampling_points - 1] = data_in.back();
    double delta = (double)(data_in.size() - 1) / (double)(n_resampling_points - 1);
    for (UInt i = 1; i < n_resampling_points - 1; ++i)
    {
      double pseudo_index = (double)i * delta;
      double left_index = (UInt)floor(pseudo_index);
      double right_index = (UInt)ceil(pseudo_index);
      if (left_index == right_index)
      {
        data_out[i] = data_in[left_index];
      }
      else
      {
        double weight_left = 1.0 - (pseudo_index - (double)left_index);
        double weight_right = 1.0 - ((double)right_index - pseudo_index);
        data_out[i] = weight_left * data_in[left_index] + weight_right * data_in[right_index];
      }
    }
  }

  void ConsensusMapNormalizerAlgorithmQuantile::extractIntensityVectors(const ConsensusMap& map, vector<vector<double> >& out_intensities)
  {
    //reserve space for out_intensities (unequal vector lengths, 0-features omitted)
    Size number_of_maps = map.getFileDescriptions().size();
    out_intensities.clear();
    out_intensities.resize(number_of_maps);
    for (UInt i = 0; i < number_of_maps; i++)
    {
      out_intensities[i].reserve(map.getFileDescriptions()[i].size);
    }
    //fill out_intensities
    ConsensusMap::ConstIterator cf_it;
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        out_intensities[f_it->getMapIndex()].push_back(f_it->getIntensity());
      }
    }
  }

  void ConsensusMapNormalizerAlgorithmQuantile::setNormalizedIntensityValues(const vector<vector<double> >& feature_ints, ConsensusMap& map)
  {
    //assumes the input map and feature_ints are in the same order as in the beginning,
    //although feature_ints has normalized values now (but the same ranks as before)
    Size number_of_maps = map.getFileDescriptions().size();
    ConsensusMap::ConstIterator cf_it;
    vector<Size> progress_indices(number_of_maps);
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        Size map_idx = f_it->getMapIndex();
        double intensity = feature_ints[map_idx][progress_indices[map_idx]++];
        f_it->asMutable().setIntensity(intensity);
      }
    }
  }
} 
