// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/MATH/MathFunctions.h>

using namespace std;

namespace OpenMS
{

void KDTreeFeatureMaps::addFeature(Size mt_map_index, const BaseFeature* feature)
{
  map_index_.push_back(mt_map_index);
  features_.push_back(feature);
  rt_.push_back(feature->getRT());

  KDTreeFeatureNode mt_node(this, size() - 1);
  kd_tree_.insert(mt_node);
}

const BaseFeature* KDTreeFeatureMaps::feature(Size i) const
{
  return features_[i];
}

double KDTreeFeatureMaps::rt(Size i) const
{
  return rt_[i];
}

double KDTreeFeatureMaps::mz(Size i) const
{
  return features_[i]->getMZ();
}

float KDTreeFeatureMaps::intensity(Size i) const
{
  return features_[i]->getIntensity();
}

Int KDTreeFeatureMaps::charge(Size i) const
{
  return features_[i]->getCharge();
}

Size KDTreeFeatureMaps::mapIndex(Size i) const
{
  return map_index_[i];
}

Size KDTreeFeatureMaps::size() const
{
  return features_.size();
}

Size KDTreeFeatureMaps::treeSize() const
{
  return kd_tree_.size();
}

Size KDTreeFeatureMaps::numMaps() const
{
  return num_maps_;
}

void KDTreeFeatureMaps::clear()
{
  features_.clear();
  map_index_.clear();
  kd_tree_.clear();
}

void KDTreeFeatureMaps::optimizeTree()
{
  kd_tree_.optimize();
}

void KDTreeFeatureMaps::getNeighborhood(Size index, vector<Size>& result_indices, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map, double max_pairwise_log_fc) const
{
  pair<double, double> rt_win = Math::getTolWindow(rt(index), rt_tol, false);
  pair<double, double> mz_win = Math::getTolWindow(mz(index), mz_tol, mz_ppm);

  vector<Size> tmp_result;
  Size ignored_map_index = include_features_from_same_map ? numeric_limits<Size>::max() : map_index_[index];
  queryRegion(rt_win.first, rt_win.second, mz_win.first, mz_win.second, tmp_result, ignored_map_index);

  if (max_pairwise_log_fc < 0.0)
  {
    result_indices.insert(result_indices.end(), tmp_result.begin(), tmp_result.end());
  }
  else // max log fold change check enabled
  {
    double int_1 = features_[index]->getIntensity();

    for (vector<Size>::const_iterator it = tmp_result.begin(); it != tmp_result.end(); ++it)
    {
      double int_2 = features_[*it]->getIntensity();
      double abs_log_fc = fabs(log10(int_2 / int_1));

      // abs_log_fc could assume +nan or +inf if negative
      // or zero intensity features were present, but
      // this shouldn't cause a problem. they just wouldn't
      // be used.
      if (abs_log_fc <= max_pairwise_log_fc)
      {
        result_indices.push_back(*it);
      }
    }
  }
}

void KDTreeFeatureMaps::queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, vector<Size>& result_indices, Size ignored_map_index) const
{
  // set up tolerance window as region for the 2D tree
  FeatureKDTree::_Region_ region;
  region._M_low_bounds[0] = rt_low;
  region._M_high_bounds[0] = rt_high;
  region._M_low_bounds[1] = mz_low;
  region._M_high_bounds[1] = mz_high;

  // range-query tolerance window
  vector<KDTreeFeatureNode> tmp_result;
  kd_tree_.find_within_range(region, back_insert_iterator<vector<KDTreeFeatureNode> >(tmp_result));

  // add indices to result
  result_indices.clear();
  for (vector<KDTreeFeatureNode>::const_iterator it = tmp_result.begin(); it != tmp_result.end(); ++it)
  {
    Size found_index = it->getIndex();
    if (ignored_map_index == numeric_limits<Size>::max() || map_index_[found_index] != ignored_map_index)
    {
      result_indices.push_back(found_index);
    }
  }
}

void KDTreeFeatureMaps::applyTransformations(const vector<TransformationModelLowess*>& trafos)
{
  for (Size i = 0; i < size(); ++i)
  {
    rt_[i] = trafos[map_index_[i]]->evaluate(features_[i]->getRT());
  }
}

void KDTreeFeatureMaps::updateMembers_()
{
}

}
