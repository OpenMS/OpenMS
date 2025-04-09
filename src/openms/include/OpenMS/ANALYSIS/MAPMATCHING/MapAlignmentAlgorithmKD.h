// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

namespace OpenMS
{

/**
    @brief An efficient reference-free feature map alignment algorithm for unlabeled data

    This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
    in a compatibility graph on feature data. This graph is comprised of nodes corresponding
    to features and edges connecting features f and f' iff both are within each other's tolerance
    windows (wrt. RT and m/z difference). CCCs are those CCs that do not contain multiple features
    from the same input map, and whose features all have the same charge state.

    All CCCs above a user-specified minimum size are considered true sets of corresponding features
    and based on these, LOWESS transformations are computed for each input map such that the average
    deviation from the mean retention time within all CCCs is minimized.

    @ingroup MapAlignment
*/

class OPENMS_DLLAPI MapAlignmentAlgorithmKD
{
public:

  /// Constructor
  MapAlignmentAlgorithmKD(Size num_maps, const Param& param);

  /// Default destructor
  virtual ~MapAlignmentAlgorithmKD();

  /// Compute data points needed for RT transformation in the current @p kd_data, add to fit_data_
  void addRTFitData(const KDTreeFeatureMaps& kd_data);

  /// Fit LOWESS to fit_data_, store final models in transformations_
  void fitLOWESS();

  /// Transform RTs for @p kd_data
  void transform(KDTreeFeatureMaps& kd_data) const;

protected:

  virtual void updateMembers_();

  /// Compute connected components, store CC indices in member cc_index. Return number of CCs.
  Size computeCCs_(const KDTreeFeatureMaps& kd_data, std::vector<Size>& cc_index) const;

  /// Return connected components
  void getCCs_(const KDTreeFeatureMaps& kd_data, std::map<Size, std::vector<Size> >& result) const;

  /// Filter connected components (return conflict-free CCs of sufficiently large size and small diameter)
  void filterCCs_(const KDTreeFeatureMaps& kd_data, const std::map<Size, std::vector<Size> >& ccs, std::map<Size, std::vector<Size> >& filtered_ccs) const;

private:

  /// Default constructor is not supposed to be used.
  MapAlignmentAlgorithmKD();

  /// RT data for fitting the LOWESS
  std::vector<TransformationModel::DataPoints> fit_data_;

  /// LOWESS transformations
  std::vector<TransformationModelLowess*> transformations_;

  /// Parameters
  Param param_;

  /// Maximum absolute log10 fold change threshold between compatible features
  double max_pairwise_log_fc_;

  /// RT tolerance
  double rt_tol_secs_;

  /// m/z tolerance
  double mz_tol_;

  /// m/z unit ppm?
  bool mz_ppm_;

};

} // namespace OpenMS

