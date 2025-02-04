// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace std;

namespace OpenMS
{
  ConsensusMapNormalizerAlgorithmQuantile::ConsensusMapNormalizerAlgorithmQuantile() = default;

  ConsensusMapNormalizerAlgorithmQuantile::~ConsensusMapNormalizerAlgorithmQuantile() = default;

  void ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(ConsensusMap& map)
  {
    //extract feature intensities
    vector<vector<double> > feature_ints;
    extractIntensityVectors(map, feature_ints);
    Size number_of_maps = feature_ints.size();

    //determine largest number of features in any map
    Size largest_number_of_features = 0;
    for (Size i = 0; i < number_of_maps; ++i)
    {
      if (feature_ints[i].size() > largest_number_of_features)
      {
        largest_number_of_features = feature_ints[i].size();
      }
    }

    //resample n data points from each sorted intensity distribution (from the different maps), n = maximum number of features in any map
    vector<vector<double> > resampled_sorted_data;
    for (Size i = 0; i < number_of_maps; ++i)
    {
      vector<double> sorted = feature_ints[i];
      std::sort(sorted.begin(), sorted.end());
      vector<double> resampled(largest_number_of_features);
      resample(sorted, resampled, static_cast<UInt>(largest_number_of_features));
      resampled_sorted_data.push_back(resampled);
    }

    //compute reference distribution from all resampled distributions
    vector<double> reference_distribution(largest_number_of_features);
    for (Size i = 0; i < number_of_maps; ++i)
    {
      for (Size j = 0; j < largest_number_of_features; ++j)
      {
        reference_distribution[j] += (resampled_sorted_data[i][j] / (double)number_of_maps);
      }
    }

    //for each map: resample from the reference distribution down to the respective original size again
    vector<vector<double> > normalized_sorted_ints(number_of_maps);
    for (Size i = 0; i < number_of_maps; ++i)
    {
      vector<double> ints;
      resample(reference_distribution, ints, static_cast<UInt>(feature_ints[i].size()));
      normalized_sorted_ints[i] = ints;
    }

    //set the intensities of feature_ints to the normalized intensities
    for (Size i = 0; i < number_of_maps; ++i)
    {
      // We do not want to change the order in feature_ints[i] but normalized_sorted_ints
      // comes sorted, so we transfer the values in feature_ints[i] into pairs that store
      // the value and the index in feature_ints[i]. Then we sort the vector of pair and as
      // a result store the indexes of feature_ints[i] in a sorted order in sort_indices.
      std::vector<std::pair<double, UInt> > sort_pairs;
      sort_pairs.reserve(feature_ints[i].size());
      for (Size j = 0; j < feature_ints[i].size(); ++j)
      {
        sort_pairs.emplace_back(feature_ints[i][j], j);
      }
      std::sort(sort_pairs.begin(), sort_pairs.end());
      vector<Size> sort_indices;
      sort_indices.reserve(sort_pairs.size());
      for (Size j = 0; j < sort_pairs.size(); ++j)
      {
        sort_indices.push_back(sort_pairs.at(j).second);
      }
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
    Size number_of_maps = map.getColumnHeaders().size();
    out_intensities.clear();
    out_intensities.resize(number_of_maps);
    for (UInt i = 0; i < number_of_maps; i++)
    {
      ConsensusMap::ColumnHeaders::const_iterator it = map.getColumnHeaders().find(i);
      if (it == map.getColumnHeaders().end()) throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(i));
      out_intensities[i].reserve(it->second.size);
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
    Size number_of_maps = map.getColumnHeaders().size();
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
