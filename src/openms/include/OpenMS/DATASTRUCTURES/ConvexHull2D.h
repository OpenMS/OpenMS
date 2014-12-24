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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CONVEXHULL2D_H
#define OPENMS_DATASTRUCTURES_CONVEXHULL2D_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief A 2-dimensional hull representation in [counter]clockwise direction - depending on axis labelling.

      The current implementation does not guarantee to produce convex hulls.
      It can still store 'old' convex hulls from featureXML without problems, but does not support the enclose() query in this case,
      and you will get an exception. As an alternative, you can use my_hull.getBoundingBox().encloses(), which yields similar results,
      and will always work.

      If you are creating new hull from peaks (e.g. during FeatureFinding), the generated hulls of a feature are defined as
      a range in m/z dimension for each RT scan (thus might be non-convex). This has the advantage that one can clearly see
      where points range within each scan (although missing points within this range are still not shown).
      When hulls are created like this, the encloses() function is supported, and works as expected, i.e.
      for the shape defined by this hull (view it in TOPPView) it answers whether the point is inside the shape.
  However, once you store the hull in featureXML and load it again, the encloses() function is not supported any longer, because
  the old convex hulls did not save min&max for each scan.
  (to support encloses() at least for the new hulls, one would need to check if there exists a min&max value for each scan
   --> then the query would be valid and the inner representation can be filled. Old featureXML's are not supported in any case.)

      The outer hullpoints can be queried by getHullPoints().

      @improvement For chromatograms we could postprocess the input and remove points in intermediate RT scans,
      which are currently reported but make the number of points rather large.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI ConvexHull2D
  {
public:
    typedef DPosition<2> PointType;
    typedef std::vector<PointType> PointArrayType;
    typedef PointArrayType::size_type SizeType;
    typedef PointArrayType::const_iterator PointArrayTypeConstIterator;

    typedef Map<PointType::CoordinateType, DBoundingBox<1> > HullPointType;

    /// default constructor
    ConvexHull2D();

    /// assignment operator
    ConvexHull2D & operator=(const ConvexHull2D & rhs);

    /// equality operator
    bool operator==(const ConvexHull2D & rhs) const;

    /// removes all points
    void clear();

    /// accessor for the outer points
    const PointArrayType & getHullPoints() const;

    /// accessor for the outer(!) points (no checking is performed if this is actually a convex hull)
    void setHullPoints(const PointArrayType & points);

    /// returns the bounding box of the feature hull points
    DBoundingBox<2> getBoundingBox() const;

    /// adds a point to the hull if it is not already contained. Returns if the point was added.
    /// this will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
    bool addPoint(const PointType & point);

    /// adds points to the hull if it is not already contained.
    /// this will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
    void addPoints(const PointArrayType & points);

    /**
      @brief Allows to reduce the disk/memory footprint of a hull

      Removes points from the hull which lie on a straight line and do not contribute to
      the hulls shape. Should be called before saving to disk.

      Example: Consider a series of 3 scans with the same dimension in m/z. After calling
      compress, the points from the second scan will be removed, since they do not contribute
      to the convex hull.

      @returns Number of removed scans
    **/
    Size compress();

    /**
              @brief Expand a convex hull to its bounding box.

              This reduces the size of a convex hull to four points, its
      bounding box, thus reducing size when storing the information.
              Note that this leads to an enclosed area that can be significantly
              larger than the original convex hull.
          **/
    void expandToBoundingBox();


    /**
              @brief returns if the @p point lies in the feature hull

              This function is only supported if the hull is created using addPoint() or addPoints(),
              but not then suing setHullPoints().
              If you require the latter functionality, then augment this function.

              @throws Exception::NotImplemented if only hull points (outer_points_), but no internal structure (map_points_) is given
          **/
    bool encloses(const PointType & point) const;

protected:
    /// internal structure maintaining the hull and enabling queries to encloses()
    HullPointType map_points_;

    /// just the list of points of the outer hull (derived from map_points_ or given by user)
    mutable PointArrayType outer_points_;

  };
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DCONVEXHULL_H
