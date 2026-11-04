// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>

#include <utility>

namespace OpenMS
{

  SpectrumAccessQuadMZTransforming::SpectrumAccessQuadMZTransforming(
      OpenSwath::SpectrumAccessPtr sptr,
      double a, double b, double c, bool ppm) :
        SpectrumAccessTransforming(std::move(sptr)), 
        a_(a), 
        b_(b), 
        c_(c), 
        ppm_(ppm)
    {}
        
    SpectrumAccessQuadMZTransforming::~SpectrumAccessQuadMZTransforming() = default;

    std::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessQuadMZTransforming::lightClone() const
    {
      // Create a light clone of *this by initializing a new
      // SpectrumAccessQuadMZTransforming with a light clone of the underlying
      // SpectrumAccess object and the parameters.
      return std::shared_ptr<SpectrumAccessQuadMZTransforming>(
          new SpectrumAccessQuadMZTransforming(sptr_->lightClone(), a_, b_, c_, ppm_));
    }

    OpenSwath::SpectrumPtr SpectrumAccessQuadMZTransforming::getSpectrumById(int id)
    {
      OpenSwath::SpectrumPtr s = sptr_->getSpectrumById(id);
      // The detached spectrum owns this array throughout the transformation.
      // Fetch it once rather than copying its shared owner for every operand.
      auto& mz = s->getMZArray()->data;
      for (double& value : mz)
      {
        // mz = a + b * mz + c * mz^2
        double predict = 
          a_ + 
          b_ * value +
          c_ * value * value;

        // If ppm is true, we predicted the ppm deviation, not the actual new mass
        if (ppm_)
        {
          value = value - predict*value/1000000;
        }
        else
        {
          value = predict;
        }
      }
      return s;
    }

}
