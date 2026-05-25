// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Feature.h>

using namespace std;

namespace OpenMS
{
  Feature::Feature() :
    BaseFeature(),
    bounding_box_(),
    convex_hulls_(),
    subordinates_()
  {
    std::fill(qualities_, qualities_ + 2, QualityType(0.0));
  }

  Feature::Feature(const BaseFeature& base) :
    BaseFeature(base)
  {
    std::fill(qualities_, qualities_ + 2, QualityType(0.0));
  }

  Feature::Feature(const Feature& feature) :
    BaseFeature(feature),
    bounding_box_(feature.bounding_box_),
    convex_hulls_(feature.convex_hulls_),
    subordinates_(feature.subordinates_)
  {
    std::copy(feature.qualities_, feature.qualities_ + 2, qualities_);
  }

  Feature::Feature(Feature&& feature) noexcept :
    BaseFeature(std::move(feature)),
    bounding_box_(feature.bounding_box_),
    convex_hulls_(std::move(feature.convex_hulls_)),
    subordinates_(std::move(feature.subordinates_))
  {
    std::copy(feature.qualities_, feature.qualities_ + 2, qualities_);
  }

  Feature::~Feature() = default;

  Feature::QualityType Feature::getOverallQuality() const
  {
    return quality_;
  }

  void Feature::setOverallQuality(Feature::QualityType q)
  {
    quality_ = q;
  }

  Feature::QualityType Feature::getQuality(Size index) const
  {
    OPENMS_PRECONDITION(index < 2, "Feature<2>::getQuality(Size): index overflow!");
    return qualities_[index];
  }

  void Feature::setQuality(Size index, Feature::QualityType q)
  {
    OPENMS_PRECONDITION(index < 2, "Feature<2>::setQuality(Size): index overflow!");
    qualities_[index] = q;
  }

  const DBoundingBox<2>& Feature::getBoundingBox() const
  {
    return bounding_box_;
  }

  void Feature::setBoundingBox(const DBoundingBox<2>& box)
  {
    bounding_box_ = box;
  }

  void Feature::setBoundingBox(double rt_min, double mz_min, double rt_max, double mz_max)
  {
    bounding_box_ = DBoundingBox<2>(
      DBoundingBox<2>::PositionType(rt_min, mz_min),
      DBoundingBox<2>::PositionType(rt_max, mz_max));
  }

  bool Feature::hasBoundingBox() const
  {
    return !bounding_box_.isEmpty();
  }

  bool Feature::updateBoundingBoxFromConvexHulls()
  {
    if (convex_hulls_.empty())
    {
      return false;
    }
    DBoundingBox<2> box;
    for (const auto& hull : convex_hulls_)
    {
      DBoundingBox<2> hb = hull.getBoundingBox();
      if (!hb.isEmpty())
      {
        box.enlarge(hb.minPosition());
        box.enlarge(hb.maxPosition());
      }
    }
    bounding_box_ = box;
    return true;
  }

  const std::vector<ConvexHull2D>& Feature::getConvexHulls() const
  {
    return convex_hulls_;
  }

  void Feature::setConvexHulls(const std::vector<ConvexHull2D>& hulls)
  {
    convex_hulls_ = hulls;
  }

  void Feature::addConvexHull(const ConvexHull2D& hull)
  {
    convex_hulls_.push_back(hull);
  }

  void Feature::clearConvexHulls()
  {
    convex_hulls_.clear();
  }

  bool Feature::encloses(double rt, double mz) const
  {
    if (!convex_hulls_.empty())
    {
      ConvexHull2D::PointType tmp(rt, mz);
      for (const ConvexHull2D& hull : convex_hulls_)
      {
        if (hull.encloses(tmp))
        {
          return true;
        }
      }
      return false;
    }
    if (hasBoundingBox())
    {
      return bounding_box_.encloses(rt, mz);
    }
    return false;
  }

  Feature& Feature::operator=(const Feature& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    BaseFeature::operator=(rhs);
    std::copy(rhs.qualities_, rhs.qualities_ + 2, qualities_);
    bounding_box_  = rhs.bounding_box_;
    convex_hulls_  = rhs.convex_hulls_;
    subordinates_  = rhs.subordinates_;

    return *this;
  }

  Feature& Feature::operator=(Feature&& rhs) & noexcept
  {
    if (this == &rhs)
    {
      return *this;
    }

    BaseFeature::operator=(std::move(rhs));
    std::copy(rhs.qualities_, rhs.qualities_ + 2, qualities_);
    bounding_box_  = rhs.bounding_box_;
    convex_hulls_  = std::move(rhs.convex_hulls_);
    subordinates_  = std::move(rhs.subordinates_);

    return *this;
  }

  bool Feature::operator==(const Feature& rhs) const
  {
    return BaseFeature::operator==(rhs)
           && equal(qualities_, qualities_ + 2, rhs.qualities_)
           && (bounding_box_ == rhs.bounding_box_)
           && (convex_hulls_ == rhs.convex_hulls_)
           && (subordinates_  == rhs.subordinates_);
  }

  const std::vector<Feature>& Feature::getSubordinates() const
  {
    return subordinates_;
  }

  std::vector<Feature>& Feature::getSubordinates()
  {
    return subordinates_;
  }

  void Feature::setSubordinates(const std::vector<Feature>& rhs)
  {
    subordinates_ = rhs;
  }

  void Feature::updateAllIDReferences(const IdentificationData::RefTranslator& trans)
  {
    updateIDReferences(trans); // update the feature itself (via BaseFeature method)
    for (Feature& sub : subordinates_) // recursively update subordinate features
    {
      sub.updateAllIDReferences(trans);
    }
  }

}
