// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MobilityPeak1D.h>

namespace OpenMS
{
  // Moved function definitions

  bool MobilityPeak1D::operator==(const MobilityPeak1D& rhs) const
  {
    return position_ == rhs.position_ && intensity_ == rhs.intensity_;
  }

  MobilityPeak1D::IntensityType MobilityPeak1D::getIntensity() const
  {
    return intensity_;
  }

  void MobilityPeak1D::setIntensity(MobilityPeak1D::IntensityType intensity)
  {
    intensity_ = intensity;
  }

  MobilityPeak1D::CoordinateType MobilityPeak1D::getMobility() const
  {
    return position_[0];
  }

  void MobilityPeak1D::setMobility(MobilityPeak1D::CoordinateType mobility)
  {
    position_[0] = mobility;
  }

  MobilityPeak1D::CoordinateType MobilityPeak1D::getPos() const
  {
    return position_[0];
  }

  void MobilityPeak1D::setPos(MobilityPeak1D::CoordinateType pos)
  {
    position_[0] = pos;
  }

  MobilityPeak1D::PositionType const& MobilityPeak1D::getPosition() const
  {
    return position_;
  }

  MobilityPeak1D::PositionType& MobilityPeak1D::getPosition()
  {
    return position_;
  }

  void MobilityPeak1D::setPosition(MobilityPeak1D::PositionType const& position)
  {
    position_ = position;
  }

  bool MobilityPeak1D::operator!=(const MobilityPeak1D& rhs) const
  {
    return !(operator==(rhs));
  }

  std::ostream & operator<<(std::ostream & os, const MobilityPeak1D & point)
  {
    os << "POS: " << point.getMobility() << " INT: " << point.getIntensity();
    return os;
  }
}
