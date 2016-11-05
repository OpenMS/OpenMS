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

#ifndef OPENMS_ANALYSIS_QUANTITATION_KDTREEDATA_H
#define OPENMS_ANALYSIS_QUANTITATION_KDTREEDATA_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeNode.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

namespace OpenMS
{

/// 2D tree on features
typedef KDTree::KDTree<2,KDTreeNode> FeatureKDTree;

/// Stores a set of features, together with a 2D tree for fast search
class OPENMS_DLLAPI KDTreeData
{

public:

  /// Constructor
  template <typename MapType>
  KDTreeData(const std::vector<MapType>& maps, const Param& param) :
    rt_tol_secs_((double)(param.getValue("rt_tol"))),
    mz_tol_((double)(param.getValue("mz_tol"))),
    mz_ppm_(param.getValue("mz_unit").toString() == "ppm"),
    num_maps_(maps.size())
  {
    for (Size i = 0; i < maps.size(); ++i)
    {
      const MapType& m = maps[i];
      for (typename MapType::const_iterator it = m.begin(); it != m.end(); ++it)
      {
        addFeature(i, &(*it));
      }
    }
    optimizeTree();
  }

  /// Destructor
  ~KDTreeData()
  {
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

  /// RT tolerance
  double rtTolerance() const;

  /// mz tolerance
  double mzTolerance() const;

  /// mz tolerance ppm?
  bool mzPPM() const;

  /// Clear all data
  void clear();

  /// Optimize the kD tree
  void optimizeTree();

  /// Fill @p result with indices of all features compatible (wrt. RT, m/z, map index) to the feature with @p index
  void getNeighborhood(Size index, std::vector<Size>& result_indices, bool ignore_map_index = false) const;

  /// Return window around @p val given tolerance @p tol
  std::pair<double, double> getTolWindow(double val, double tol, bool ppm) const;

  /// Apply RT transformations
  void applyTransformations(const std::vector<TransformationModelLowess*>& trafos);

protected:

  /// Feature data
  std::vector<const BaseFeature*> features_;

  /// Map indices
  std::vector<Size> map_index_;

  /// (Potentially transformed) retention times
  std::vector<double> rt_;

  /// Number of maps
  Size num_maps_;

  /// RT tolerance in seconds
  double rt_tol_secs_;

  /// m/z tolerance in ppm or Da
  double mz_tol_;

  /// m/z tolerance unit ppm?
  bool mz_ppm_;

  /// 2D tree on features from all input maps.
  FeatureKDTree kd_tree_;

};
}

#endif // OPENMS_ANALYSIS_QUANTITATION_KDTREEDATA_H
