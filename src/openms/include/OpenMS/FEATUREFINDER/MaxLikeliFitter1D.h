// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/Fitter1D.h>

namespace OpenMS
{
  class InterpolationModel;

  /**
  @brief Abstract base class for all 1D-model fitters using maximum likelihood optimization.
  */
  class OPENMS_DLLAPI MaxLikeliFitter1D : public Fitter1D
  {
  public:
    /// default constructor
    MaxLikeliFitter1D() : Fitter1D()
    {
    }

    /// copy constructor
    MaxLikeliFitter1D(const MaxLikeliFitter1D& source) : Fitter1D(source)
    {
    }

    /// destructor
    ~MaxLikeliFitter1D() override
    {
    }

    /// assignment operator
    MaxLikeliFitter1D& operator=(const MaxLikeliFitter1D& source)
    {
      if (&source == this)
        return *this;

      Fitter1D::operator=(source);

      return *this;
    }

  protected:
    /// fit an offset on the basis of the Pearson correlation coefficient
    QualityType fitOffset_(std::unique_ptr<InterpolationModel>& model, const RawDataArrayType& set, const CoordinateType stdev1, const CoordinateType stdev2, const CoordinateType offset_step) const;

    void updateMembers_() override;
  };
} // namespace OpenMS
