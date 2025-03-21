// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>

using namespace std;

namespace OpenMS
{

  SplinePackage::SplinePackage(std::vector<double> pos, const std::vector<double>& intensity) :
    spline_(pos, intensity)
  {
    if (!(pos.size() == intensity.size() && pos.size() > 1))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z (or RT) and intensity vectors either not of the same size or too short.");
    }

    pos_min_ = pos.front();
    pos_max_ = pos.back();
    pos_step_width_ = (pos_max_ - pos_min_) / (pos.size() - 1);
  }

  SplinePackage::~SplinePackage() = default;

  double SplinePackage::getPosMin() const
  {
    return pos_min_;
  }

  double SplinePackage::getPosMax() const
  {
    return pos_max_;
  }

  double SplinePackage::getPosStepWidth() const
  {
    return pos_step_width_;
  }

  bool SplinePackage::isInPackage(double pos) const
  {
    return pos >= pos_min_ && pos <= pos_max_;
  }

  double SplinePackage::eval(double pos) const
  {
    if (this->isInPackage(pos))
    {
      return max(0.0, spline_.eval(pos));
    }
    else
    {
      return 0;
    }
  }

}
