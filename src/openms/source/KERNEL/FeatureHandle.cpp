// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureHandle.h>

#include <OpenMS/KERNEL/BaseFeature.h>

namespace OpenMS
{
  FeatureHandle::FeatureHandle() :
    Peak2D(),
    UniqueIdInterface(),
    map_index_(0),
    charge_(0),
    width_(0)
  {
  }

  FeatureHandle::FeatureHandle(UInt64 map_index, const Peak2D& point, UInt64 element_index) :
    Peak2D(point),
    map_index_(map_index),
    charge_(0),
    width_(0)
  {
    setUniqueId(element_index);
  }

  FeatureHandle::FeatureHandle(UInt64 map_index, const BaseFeature& feature) :
    Peak2D(feature),
    UniqueIdInterface(feature),
    map_index_(map_index),
    charge_(feature.getCharge()),
    width_(feature.getWidth())
  {
  }

  FeatureHandle::FeatureHandle(const FeatureHandle& rhs) = default;

  FeatureHandle& FeatureHandle::operator=(const FeatureHandle& rhs) = default;

  FeatureHandle::~FeatureHandle() = default;

  UInt64 FeatureHandle::getMapIndex() const
  {
    return map_index_;
  }

  void FeatureHandle::setMapIndex(UInt64 i)
  {
    map_index_ = i;
  }

  void FeatureHandle::setCharge(FeatureHandle::ChargeType charge)
  {
    charge_ = charge;
  }

  FeatureHandle::ChargeType FeatureHandle::getCharge() const
  {
    return charge_;
  }

  void FeatureHandle::setWidth(FeatureHandle::WidthType width)
  {
    width_ = width;
  }

  FeatureHandle::WidthType FeatureHandle::getWidth() const
  {
    return width_;
  }

  bool FeatureHandle::operator==(const FeatureHandle& i) const
  {
    return (Peak2D::operator==(i))
           && (UniqueIdInterface::operator==(i))
           && (map_index_ == i.map_index_)
           && (charge_ == i.charge_)
           && (width_ == i.width_);
  }

  bool FeatureHandle::operator!=(const FeatureHandle& i) const
  {
    return !(operator==(i));
  }

  std::ostream& operator<<(std::ostream& os, const FeatureHandle& cons)
  {
    os << "---------- FeatureHandle -----------------\n"
       << "RT: " << cons.getRT() << std::endl
       << "m/z: " << cons.getMZ() << std::endl
       << "Intensity: " << cons.getIntensity() << std::endl
       << "Map Index: " << cons.getMapIndex() << std::endl
       << "Element Id: " << cons.getUniqueId() << std::endl;
    return os;
  }

}
