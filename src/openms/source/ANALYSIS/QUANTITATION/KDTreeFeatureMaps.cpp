// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <boost/container/vector.hpp>
#include <queue>

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

double KDTreeFeatureMaps::rtTolerance() const
{
  return rt_tol_secs_;
}

double KDTreeFeatureMaps::mzTolerance() const
{
  return mz_tol_;
}

bool KDTreeFeatureMaps::mzPPM() const
{
  return mz_ppm_;
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

void KDTreeFeatureMaps::getNeighborhood(Size index, vector<Size>& result_indices, bool ignore_map_index) const
{
  pair<double, double> rt_win = Math::getTolWindow(rt(index), rt_tol_secs_, false);
  pair<double, double> mz_win = Math::getTolWindow(mz(index), mz_tol_, mz_ppm_);

  // set up tolerance window as region for the 2D tree
  FeatureKDTree::_Region_ region;
  region._M_low_bounds[0] = rt_win.first;
  region._M_high_bounds[0] = rt_win.second;
  region._M_low_bounds[1] = mz_win.first;
  region._M_high_bounds[1] = mz_win.second;

  // range-query tolerance window
  vector<KDTreeFeatureNode> tmp_result;
  kd_tree_.find_within_range(region, back_insert_iterator<vector<KDTreeFeatureNode> >(tmp_result));

  // unless ignore_map_index: add only compatible MTs from *other* maps to final result
  result_indices.clear();
  for (vector<KDTreeFeatureNode>::const_iterator it = tmp_result.begin(); it != tmp_result.end(); ++it)
  {
    Size found_index = it->getIndex();
    if (ignore_map_index || map_index_[found_index] != map_index_[index])
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

}
