// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

namespace OpenMS
{

  ConvexHull2D::ConvexHull2D() :
    map_points_(),
    outer_points_()
  {
  }

  /// assignment operator
  ConvexHull2D& ConvexHull2D::operator=(const ConvexHull2D& rhs)
  {
    if (&rhs == this)
    {
      return *this;
    }
    map_points_ = rhs.map_points_;
    outer_points_ = rhs.outer_points_;

    return *this;
  }

  /// equality operator
  bool ConvexHull2D::operator==(const ConvexHull2D& rhs) const
  {
    // different size => return false
    if (map_points_.size() != rhs.map_points_.size())
    {
      return false;
    }
    if (outer_points_.size() != rhs.outer_points_.size())
    {
      return false;
    }
    //different points now => return false
    for (HullPointType::const_iterator it = rhs.map_points_.begin(); it !=   rhs.map_points_.end(); ++it)
    {
      if (map_points_.find(it->first) != map_points_.end())
      {
        if (map_points_.at(it->first) != it->second)
        {
          return false;
        }
      }
      else
      {
        return false;
      }
    }
    //different points now => return false
    for (Size i = 0; i < rhs.outer_points_.size(); ++i)
    {
      if (outer_points_[i] != rhs.outer_points_[i])
      {
        return false;
      }
    }

    return true;
  }

  /// removes all points
  void ConvexHull2D::clear()
  {
    map_points_.clear();
    outer_points_.clear();
  }

  /// accessor for the points
  const ConvexHull2D::PointArrayType& ConvexHull2D::getHullPoints() const
  {
    // construct outer hull if required
    if (outer_points_.empty() && !map_points_.empty())
    {
      // walk the outer hull
      outer_points_.reserve(map_points_.size() * 2);

      // traverse lower m/z's of RT scans
      for (HullPointType::const_iterator it = map_points_.begin(); it != map_points_.end(); ++it)
      {
        PointType p;
        p.setX(it->first);
        p.setY(it->second.minPosition()[0]);
        outer_points_.push_back(p);
      }

      // traverse higher m/z's of RT scans
      for (HullPointType::const_reverse_iterator it = map_points_.rbegin(); it != map_points_.rend(); ++it)
      {
        PointType p;
        p.setX(it->first);
        p.setY(it->second.maxPosition()[0]);
        // turning point (avoid listing it twice if last scan only has a single point)
        if ((it == map_points_.rbegin()) && (it->second.width() == 0))
        {
          continue;
        }
        // do not list first scan again if it's only a single point
        else if (it == --map_points_.rend() && (it->second.width() == 0))
        {
          continue;
        }
        outer_points_.push_back(p);
      }
    }
    return outer_points_;
  }

  void ConvexHull2D::setHullPoints(const ConvexHull2D::PointArrayType& points)
  {
    map_points_.clear();
    outer_points_ = points;
  }

  void ConvexHull2D::expandToBoundingBox()
  {
    DBoundingBox<2> bb(getBoundingBox());
    typedef DBoundingBox<2>::PositionType Point;
    clear();
    addPoint(Point(bb.minPosition()[0], bb.minPosition()[1]));
    addPoint(Point(bb.minPosition()[0], bb.maxPosition()[1]));
    addPoint(Point(bb.maxPosition()[0], bb.minPosition()[1]));
    addPoint(Point(bb.maxPosition()[0], bb.maxPosition()[1]));
  }

  /// returns the bounding box of the convex hull points
  DBoundingBox<2> ConvexHull2D::getBoundingBox() const
  {
    DBoundingBox<2> bb;

    // the internal structure might not be defined, but we try it first
    if (!map_points_.empty())
    {
      for (HullPointType::const_iterator it = map_points_.begin(); it != map_points_.end(); ++it)
      {
        bb.enlarge(it->first, it->second.minPosition()[0]);
        bb.enlarge(it->first, it->second.maxPosition()[0]);
      }
    }
    else if (!outer_points_.empty())
    {
      for (PointArrayType::const_iterator it = outer_points_.begin(); it != outer_points_.end(); ++it)
      {
        bb.enlarge((*it)[0], (*it)[1]);
      }
    }

    return bb;
  }

  bool ConvexHull2D::addPoint(const PointType& point)
  {
    outer_points_.clear();

    if (map_points_.find(point[0]) != map_points_.end())
    {
      if (map_points_.at(point[0]).encloses(point[1]))
      {
        return false;
      }
      map_points_[point[0]].enlarge(point[1]);
    }
    else
    {
      map_points_[point[0]] = DBoundingBox<1>(point[1], point[1]);
    }

    return true;
  }

  void ConvexHull2D::addPoints(const PointArrayType& points)
  {
    for (PointArrayTypeConstIterator it = points.begin(); it != points.end(); ++it)
    {
      addPoint(*it);
    }
  }

  Size ConvexHull2D::compress()
  {
    // iterate over rt scans and check if the m/z span is always the same in consecutive scans
    // keep the min&max scan only
    //
    if (map_points_.size() < 3)
    {
      return 0; // we need at least one "middle" scan
    }
    HullPointType compressed_map;

    compressed_map[map_points_.begin()->first] = map_points_.begin()->second; // copy first scan
    HullPointType::const_iterator pred_it = map_points_.begin();
    HullPointType::const_iterator middle_it = pred_it; middle_it++;
    HullPointType::const_iterator succ_it = pred_it; succ_it++; succ_it++;

    for (Size p = 1; p < map_points_.size() - 1; ++p)
    {
      if (pred_it->second == middle_it->second && middle_it->second == succ_it->second)
      {
        // middle is identical in m/z range .. do not add to the compressed_map
      }
      else
      {
        compressed_map[middle_it->first] = middle_it->second;
      }
      ++succ_it;
      ++middle_it;
      ++pred_it;
    }
    compressed_map[middle_it->first] = middle_it->second; // copy last scan
    if (succ_it != map_points_.end())
    {
      throw Exception::BufferOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    //std::cout << "compressed CH from " << map_points_.size() << " to " << compressed_map.size() << "\n";
    Size saved_points = map_points_.size() - compressed_map.size();
    //copy
    map_points_.clear();
    map_points_.insert(compressed_map.begin(), compressed_map.end());
    return saved_points;
  }

  bool ConvexHull2D::encloses(const PointType& point) const
  {
    if ((map_points_.empty()) && !outer_points_.empty()) // we cannot answer the query as we lack the internal data structure
    { // (if you need this you need to augment encloses() to work on outer_points_ only)
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    if (map_points_.find(point[0]) != map_points_.end())
    {
      if (map_points_.at(point[0]).encloses(point[1]))
      {
        return true;
      }
    }

    // find the two RT scans surrounding the point:
    HullPointType::const_iterator it_upper = map_points_.end(), it_lower = map_points_.end();
    // iterate over keys (which are sorted by ascending RT)
    for (HullPointType::const_iterator it = map_points_.begin(); it != map_points_.end(); ++it)
    {
      // lower bound
      if (((it->first) < (point[0])))
      {
        it_lower = it;
      }
      // upper bound
      if ((it_upper == map_points_.end()) && ((it->first) > (point[0])))
        {
          it_upper = it;
        }
    }

    // point is not between two scans
    if ((it_lower == map_points_.end()) || (it_upper == map_points_.end()))
    {
      return false;
    }
    // check if point is within bounds
    double mz_low = it_lower->second.minPosition()[0] // m/z offset
                    + ((point[0] - (it_lower->first)) / ((it_upper->first) - (it_lower->first)))      // factor (0-1)
                    * (it_upper->second.minPosition()[0] - it_lower->second.minPosition()[0]);                           // m/z range

    double mz_high = it_lower->second.maxPosition()[0] // m/z offset
                     + ((point[0] - (it_lower->first)) / ((it_upper->first) - (it_lower->first)))     // factor (0-1)
                     * (it_upper->second.maxPosition()[0] - it_lower->second.maxPosition()[0]);                          // m/z range


    DBoundingBox<1> range(mz_low, mz_high);
    return range.encloses(point[1]);
  }

} // namespace OpenMS
