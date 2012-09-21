// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
      gsl_sort(&sorted.front(), 1, sorted.size());
      vector<double> resampled(largest_number_of_features);
      resample(sorted, resampled, largest_number_of_features);
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
      resample(reference_distribution, ints, feature_ints[i].size());
      normalized_sorted_ints[i] = ints;
    }

    //set the intensities of feature_ints to the normalized intensities
    for (Size i = 0; i < number_of_maps; ++i)
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
