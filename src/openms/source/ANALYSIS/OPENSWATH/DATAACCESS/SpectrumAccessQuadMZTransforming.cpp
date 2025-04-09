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

    boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessQuadMZTransforming::lightClone() const
    {
      // Create a light clone of *this by initializing a new
      // SpectrumAccessQuadMZTransforming with a light clone of the underlying
      // SpectrumAccess object and the parameters.
      return boost::shared_ptr<SpectrumAccessQuadMZTransforming>(
          new SpectrumAccessQuadMZTransforming(sptr_->lightClone(), a_, b_, c_, ppm_));
    }

    OpenSwath::SpectrumPtr SpectrumAccessQuadMZTransforming::getSpectrumById(int id)
    {
      OpenSwath::SpectrumPtr s = sptr_->getSpectrumById(id);
      for (size_t i = 0; i < s->getMZArray()->data.size(); i++)
      {
        // mz = a + b * mz + c * mz^2
        double predict = 
          a_ + 
          b_ * s->getMZArray()->data[i] +
          c_ * s->getMZArray()->data[i] * s->getMZArray()->data[i];

        // If ppm is true, we predicted the ppm deviation, not the actual new mass
        if (ppm_)
        {
          s->getMZArray()->data[i] = s->getMZArray()->data[i] - predict*s->getMZArray()->data[i]/1000000;
        }
        else
        {
          s->getMZArray()->data[i] = predict;
        }
      }
      return s;
    }

}
