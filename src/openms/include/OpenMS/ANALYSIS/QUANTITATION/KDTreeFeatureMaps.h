// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_ANALYSIS_QUANTITATION_KDTREEFEATUREMAPS_H
#define OPENMS_ANALYSIS_QUANTITATION_KDTREEFEATUREMAPS_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>

namespace OpenMS
{

/// Stores a set of features, together with a 2D tree for fast search
class OPENMS_DLLAPI KDTreeFeatureMaps : public DefaultParamHandler
{

public:

  /// 2D tree on features
  typedef KDTree::KDTree<2,KDTreeFeatureNode> FeatureKDTree;

  /// Default constructor
  KDTreeFeatureMaps() :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
  }

  /// Constructor
  template <typename MapType>
  KDTreeFeatureMaps(const std::vector<MapType>& maps, const Param& param) :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
    setParameters(param);
    addMaps(maps);
  }

  /// Destructor
  ~KDTreeFeatureMaps() override
  {
  }

  /// Add @p maps and balance kd-tree
  template <typename MapType>
  void addMaps(const std::vector<MapType>& maps)
  {
    num_maps_ = maps.size();

    for (Size i = 0; i < num_maps_; ++i)
    {
      const MapType& m = maps[i];
      for (typename MapType::const_iterator it = m.begin(); it != m.end(); ++it)
      {
        addFeature(i, &(*it));
      }
    }
    optimizeTree();
  }

  /// Add feature
  void addFeature(Size mt_map_index, const BaseFeature* feature);

  /// Return pointer to feature i
  const BaseFeature* feature(Size i) const;

  /// RT
  double rt(Size i) const;

  /// m/z
  double mz(Size i) const;

  /// Intensity
  float intensity(Size i) const;

  /// Charge
  Int charge(Size i) const;

  /// Map index
  Size mapIndex(Size i) const;

  /// Number of features stored
  Size size() const;

  /// Number of points in the tree
  Size treeSize() const;

  /// Number of maps
  Size numMaps() const;

  /// Clear all data
  void clear();

  /// Optimize the kD tree
  void optimizeTree();

  /// Fill @p result with indices of all features compatible (wrt. RT, m/z, map index) to the feature with @p index
  void getNeighborhood(Size index, std::vector<Size>& result_indices, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map = false, double max_pairwise_log_fc = -1.0) const;

  /// Fill @p result with indices of all features within the specified boundaries
  void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, std::vector<Size>& result_indices, Size ignored_map_index = std::numeric_limits<Size>::max()) const;

  /// Apply RT transformations
  void applyTransformations(const std::vector<TransformationModelLowess*>& trafos);

protected:

  void updateMembers_() override;

  /// Feature data
  std::vector<const BaseFeature*> features_;

  /// Map indices
  std::vector<Size> map_index_;

  /// (Potentially transformed) retention times
  std::vector<double> rt_;

  /// Number of maps
  Size num_maps_;

  /// 2D tree on features from all input maps.
  FeatureKDTree kd_tree_;

};
}

#endif // OPENMS_ANALYSIS_QUANTITATION_KDTREEFEATUREMAPS_H
