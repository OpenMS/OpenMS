// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /**
      @brief A D-dimensional half-open interval.

      This class describes a range in D-dimensional space delimited
      by two points (i.e. a D-dimensional hyper-rectangle). The
      two points define the lower left and the upper right corner
      in 2D and analogous points in higher dimensions.

      A range is a pair of positions in D-space represented by DPosition.
      The two limiting points are accessed as minPosition() and maxPosition().

      A range denotes a semi-open interval. A lower coordinate of each
      dimension is part of the range, the higher coordinate is not.

      @ingroup Datastructures
  */
  template <UInt D>
  class DRange :
    public Internal::DIntervalBase<D>
  {
public:

    /**
        @name Type definitions
    */
    //@{
    /// Dimensions
    enum {DIMENSION = D};
    /// Base class type
    typedef Internal::DIntervalBase<D> Base;
    /// Position type
    typedef typename Base::PositionType PositionType;
    /// Coordinate type of the positions
    typedef typename Base::CoordinateType CoordinateType;
    /// Types that describe the kind of intersection between two ranges
    enum DRangeIntersection
    {
      Disjoint, ///< No intersection
      Intersects, ///< Intersection
      Inside ///< One contains the other
    };

    //@}

    using Base::min_;
    using Base::max_;

    /**@name Constructors and Destructor */
    //@{
    /**
        @brief Default constructor.

        Creates a range with all coordinates zero.
    */
    DRange() :
      Base()
    {
    }

    /// Constructor that takes two Points and constructs a range.
    DRange(const PositionType& lower, const PositionType& upper) :
      Base(lower, upper)
    {
    }

    /// Copy constructor
    DRange(const DRange& range) :
      Base(range)
    {
    }

    /// Move constructor
    DRange(DRange&&) noexcept = default;

    /// Copy constructor for the base class
    DRange(const Base& range) :
      Base(range)
    {
    }

    /// Convenient constructor for DRange<2>
    DRange(CoordinateType minx, CoordinateType miny, CoordinateType maxx, CoordinateType maxy)
    {
      static_assert(D == 2);
      min_[0] = minx;
      min_[1] = miny;
      max_[0] = maxx;
      max_[1] = maxy;
      Base::normalize_();
    }

    /// Assignment operator
    DRange& operator=(const DRange& rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Assignment operator for the base class
    DRange& operator=(const Base& rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Destructor
    ~DRange()
    {
    }

    //@}

    /**@name Predicates */
    //@{
    ///Equality operator
    bool operator==(const DRange& rhs) const
    {
      return Base::operator==(rhs);
    }

    /// Equality operator
    bool operator==(const Base& rhs) const
    {
      return Base::operator==(rhs);
    }

    /**
         @brief Checks whether this range contains a certain point.

         @param position The point's position.
         @returns true if point lies inside this area.
    */
    bool encloses(const PositionType& position) const
    {
      for (UInt i = 0; i != D; i++)
      {
        if (position[i] < min_[i]) return false;

        if (position[i] >= max_[i]) return false;
      }
      return true;
    }

    ///@brief 2D-version of encloses for convenience only
    bool encloses(CoordinateType x, CoordinateType y) const
    {
      if (x < min_[0]) return false;

      if (x >= max_[0]) return false;

      if (y < min_[1]) return false;

      if (y >= max_[1]) return false;

      return true;
    }

    /// Returns the smallest range containing this range and @p other_range
    DRange united(const DRange<D>& other_range) const
    {
      PositionType united_min;
      PositionType united_max;
      DRange<D> united_range = DRange<D>::empty;

      PositionType other_min = other_range.minPosition();
      PositionType other_max = other_range.maxPosition();

      for (Size i = 0; i != D; ++i)
      {
        united_min[i] = min_[i] < other_min[i] ? min_[i] : other_min[i];
        united_max[i] = max_[i] > other_max[i] ? max_[i] : other_max[i];
      }
      united_range.setMinMax(united_min, united_max);

      return united_range;
    }

    /**
         @brief Checks how this range intersects with another @p range.

         @param range The max_ range.
    */
    DRangeIntersection intersects(const DRange& range) const
    {
      //check if r.min_ is in this area
      if (encloses(range.min_))
      {
        //check if r.max_ in this area => Inside / Intersects
        for (Size i = 0; i != D; i++)
        {
          if (range.max_[i] > max_[i])
          {
            return Intersects;
          }
        }
        return Inside;
      }
      // => r.min_ is not inside this area
      //check if any r.min_ >= max_ => Disjoint
      for (Size i = 0; i != D; i++)
      {
        if (range.min_[i] >= max_[i])
        {
          return Disjoint;
        }
      }
      // => some coordinate of r.min_ has to be smaller than the one of min_
      //check if all coords of r are smaller than the those of the range
      for (Size i = 0; i != D; i++)
      {
        if (range.max_[i] <= min_[i])
        {
          return Disjoint;
        }
      }
      return Intersects;
    }

    /**
         @brief Checks whether this range intersects with another @p range.

         @param range The max_ range.
         @returns True if the areas intersect (i.e. they intersect or one contains the other).
    */
    bool isIntersected(const DRange& range) const
    {
      //check if r.min_ is in this area
      if (encloses(range.min_))
      {
        return true;
      }

      // => r.min_ is not inside this area
      //check if any r.min_ >= max_ => Disjoint
      for (Size i = 0; i != D; i++)
      {
        if (range.min_[i] >= max_[i])
        {
          return false;
        }
      }
      // => some coordinate of r.min_ has to be smaller than the one of min_
      //check if all coords of r are smaller than the those of the range
      for (Size i = 0; i != D; i++)
      {
        if (range.max_[i] <= min_[i])
        {
          return false;
        }
      }
      return true;
    }

    /**
         @brief Extends the range in all dimensions by a certain multiplier.

         Extends the range, while maintaining the original center position.

         Examples (for D=1):
           factor = 1.01 extends the range by 1% in total, i.e. 0.5% left and right.
           factor = 2.00 doubles the total range, e.g. from [0,100] to [-50,150]

         @param factor Multiplier (allowed is [0, inf)).
         @return A reference to self
    */
    DRange<D>& extend(double factor)
    {
      if (factor < 0)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DRange::extend(): factor must not be negative!");
      }

      for (UInt i = 0; i != D; ++i)
      {
        Internal::DIntervalBase<1>::CoordinateType extra = (max_[i] - min_[i]) / 2.0 * (factor - 1);
        min_[i] -= extra;
        max_[i] += extra;
      }
      return *this;
    }

    /**
     @brief Extends the range in all dimensions by a certain amount.

     Extends the range, while maintaining the original center position.
     If a negative @p addition is given, the range shrinks and may result in min==max (but never min>max).

     Examples (for D=1):
       addition = 0.5 extends the range by 1 in total, i.e. 0.5 left and right.
   
     @param addition Additive for each dimension (can be negative). Resulting invalid min/max are not fixed automatically!
     @return A reference to self
    */
    DRange<D>& extend(typename Base::PositionType addition)
    {
      addition /= 2;
      min_ -= addition;
      max_ += addition;
      for (UInt i = 0; i != D; ++i)
      {
        // invalid range --> reduce to single center point
        if (min_[i] > max_[i]) min_[i] = max_[i] = (min_[i] + max_[i]) / 2;
      }
      return *this;
    }

    DRange<D>& ensureMinSpan(typename Base::PositionType min_span)
    {
      typename Base::PositionType extend_by {};
      for (UInt i = 0; i != D; ++i)
      {
        // invalid range --> reduce to single center point
        if (max_[i] - min_[i] < min_span[i])
        {
          extend_by[i] = min_span[i] - (max_[i] - min_[i]); // add whatever is missing to get to min_span
        }
      }
      extend(extend_by);
      return *this;
    }

    /// swaps dimensions for 2D data (i.e. x and y coordinates)
    DRange<D>& swapDimensions()
    {
      static_assert(D==2);
      std::swap(min_[0], min_[1]);
      std::swap(max_[0], max_[1]);
      return *this;
    }

    /**
     * @brief Make sure @p point is inside the current area
     * @param point A point potentially outside the current range, which will be pulled into the current range.
     */
    void pullIn(DPosition<D>& point) const
    {
      for (UInt i = 0; i != D; ++i)
      {
        point[i] = std::max(min_[i], std::min(point[i], max_[i]));
      }
    }

    //@}
  };

  ///Print the contents to a stream.
  template <UInt D>
  std::ostream& operator<<(std::ostream& os, const DRange<D>& area)
  {
    os << "--DRANGE BEGIN--\n";
    os << "MIN --> " << area.min_ << '\n';
    os << "MAX --> " << area.max_ << '\n';
    os << "--DRANGE END--\n";
    return os;
  }

} // namespace OpenMS

