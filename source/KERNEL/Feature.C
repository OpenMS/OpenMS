// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Feature.h>

using namespace std;

namespace OpenMS
{
  Feature::Feature() :
    BaseFeature(),
    convex_hulls_(),
    convex_hulls_modified_(true),
    convex_hull_(),
    subordinates_()
  {
    std::fill(qualities_, qualities_ + 2, QualityType(0.0));
  }

  Feature::Feature(const Feature & feature) :
    BaseFeature(feature),
    model_desc_(feature.model_desc_),
    convex_hulls_(feature.convex_hulls_),
    convex_hulls_modified_(feature.convex_hulls_modified_),
    convex_hull_(feature.convex_hull_),
    subordinates_(feature.subordinates_)
  {
    std::copy(feature.qualities_, feature.qualities_ + 2, qualities_);
  }

  Feature::~Feature()
  {}

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

  const ModelDescription<2> & Feature::getModelDescription() const
  {
    return model_desc_;
  }

  ModelDescription<2> & Feature::getModelDescription()
  {
    return model_desc_;
  }

  void Feature::setModelDescription(const ModelDescription<2> & q)
  {
    model_desc_ = q;
  }

  const std::vector<ConvexHull2D> & Feature::getConvexHulls() const
  {
    return convex_hulls_;
  }

  std::vector<ConvexHull2D> & Feature::getConvexHulls()
  {
    convex_hulls_modified_ = true;
    return convex_hulls_;
  }

  void Feature::setConvexHulls(const std::vector<ConvexHull2D> & hulls)
  {
    convex_hulls_modified_ = true;
    convex_hulls_ = hulls;
  }

  ConvexHull2D & Feature::getConvexHull() const
  {
    //recalculate convex hull if necessary
    if (convex_hulls_modified_)
    {
      //only one mass trace convex hull => use it as overall convex hull
      if (convex_hulls_.size() == 1)
      {
        convex_hull_ = convex_hulls_[0];
      }
      else
      {
        convex_hull_.clear();
        if (convex_hulls_.size() > 0)
        {
          /*
          -- this does not work with our current approach of "non-convex"hull computation as the mass traces of features cannot be combined
          -- meaningfully. We thus print only the bounding box of the traces (for now)

          for (Size hull=0; hull<convex_hulls_.size(); ++hull)
          {
              convex_hull_.addPoints(convex_hulls_[hull].getHullPoints());
          }
          */

          DBoundingBox<2> box;
          for (Size hull = 0; hull < convex_hulls_.size(); ++hull)
          {
            box.enlarge(convex_hulls_[hull].getBoundingBox().minPosition()[0], convex_hulls_[hull].getBoundingBox().minPosition()[1]);
            box.enlarge(convex_hulls_[hull].getBoundingBox().maxPosition()[0], convex_hulls_[hull].getBoundingBox().maxPosition()[1]);
          }
          convex_hull_.addPoint(ConvexHull2D::PointType(box.minX(), box.minY()));
          convex_hull_.addPoint(ConvexHull2D::PointType(box.maxX(), box.minY()));
          convex_hull_.addPoint(ConvexHull2D::PointType(box.minX(), box.maxY()));
          convex_hull_.addPoint(ConvexHull2D::PointType(box.maxX(), box.maxY()));
        }

      }

      convex_hulls_modified_ = false;
    }

    return convex_hull_;
  }

  bool Feature::encloses(DoubleReal rt, DoubleReal mz) const
  {
    ConvexHull2D::PointType tmp(rt, mz);
    for (vector<ConvexHull2D>::const_iterator   it = convex_hulls_.begin(); it != convex_hulls_.end(); ++it)
    {
      if (it->encloses(tmp))
        return true;
    }
    return false;
  }

  Feature & Feature::operator=(const Feature & rhs)
  {
    if (this == &rhs)
      return *this;

    BaseFeature::operator=(rhs);
    copy(rhs.qualities_, rhs.qualities_ + 2, qualities_);
    model_desc_                         = rhs.model_desc_;
    convex_hulls_                       = rhs.convex_hulls_;
    convex_hulls_modified_  = rhs.convex_hulls_modified_;
    convex_hull_                = rhs.convex_hull_;
    subordinates_                       = rhs.subordinates_;

    return *this;
  }

  bool Feature::operator==(const Feature & rhs) const
  {
    return BaseFeature::operator==(rhs)
           && equal(qualities_, qualities_ + 2, rhs.qualities_)
           && (model_desc_ == rhs.model_desc_)
           && (convex_hulls_ == rhs.convex_hulls_)
           && (subordinates_  == rhs.subordinates_);
  }

  const std::vector<Feature> & Feature::getSubordinates() const
  {
    return subordinates_;
  }

  std::vector<Feature> & Feature::getSubordinates()
  {
    return subordinates_;
  }

  void Feature::setSubordinates(const std::vector<Feature> & rhs)
  {
    subordinates_ = rhs;
  }
}
