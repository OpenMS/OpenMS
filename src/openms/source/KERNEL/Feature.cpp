// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/CONCEPT/Macros.h> // for OPENMS_PRECONDITION
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h> // for ConvexHull2D, etc
#include <OpenMS/KERNEL/BaseFeature.h> // for BaseFeature::QualityType, etc

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
    convex_hulls_(feature.convex_hulls_),
    convex_hulls_modified_(feature.convex_hulls_modified_),
    convex_hull_(feature.convex_hull_),
    subordinates_(feature.subordinates_)
  {
    std::copy(feature.qualities_, feature.qualities_ + 2, qualities_);
  }

  Feature::~Feature()
  {
  }

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

  bool Feature::encloses(double rt, double mz) const
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
