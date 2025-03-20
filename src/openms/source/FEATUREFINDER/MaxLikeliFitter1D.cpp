// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MaxLikeliFitter1D.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{

  void MaxLikeliFitter1D::updateMembers_()
  {
    Fitter1D::updateMembers_();
  }

  MaxLikeliFitter1D::QualityType MaxLikeliFitter1D::fitOffset_(std::unique_ptr<InterpolationModel>& model,
                                                               const RawDataArrayType & set,
                                                               const CoordinateType stdev1,
                                                               const CoordinateType stdev2,
                                                               const CoordinateType offset_step) const
  {
    const CoordinateType offset_min = model->getInterpolation().supportMin() - stdev1;
    const CoordinateType offset_max = model->getInterpolation().supportMin() + stdev2;

    CoordinateType offset;
    QualityType correlation;

    //test model with default offset
    std::vector<float> real_data;
    real_data.reserve(set.size());
    std::vector<float> model_data;
    model_data.reserve(set.size());

    for (Size i = 0; i < set.size(); ++i)
    {
      real_data.push_back(set[i].getIntensity());
      model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
    }

    CoordinateType max_offset = model->getInterpolation().getOffset();
    QualityType max_correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

    //test different offsets
    for (offset = offset_min; offset <= offset_max; offset += offset_step)
    {
      // set offset
      model->setOffset(offset);

      // get samples
      model_data.clear();
      for (Size i = 0; i < set.size(); ++i)
      {
        model_data.push_back(model->getIntensity(DPosition<1>(set[i].getPosition())));
      }

      correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

      if (correlation > max_correlation)
      {
        max_correlation = correlation;
        max_offset = offset;
      }
    }

    model->setOffset(max_offset);

    return max_correlation;
  }

} // namespace OpenMS
