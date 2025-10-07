// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

  /**
      @brief A D-dimensional bounding box.

      A DBoundingBox denotes a closed interval.  Upper and lower margins are both contained.

      @ingroup Datastructures
  */
  template <UInt D>
  class DBoundingBox :
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
    //@}


    // for convenience
    using Base::min_;
    using Base::max_;

    /**@name Constructors and Destructor */
    //@{
    ///Default constructor.
    DBoundingBox() :
      Base()
    {
    }

    /// Copy constructor
    DBoundingBox(const DBoundingBox& rhs) :
      Base(rhs)
    {
    }

    /// Assignment operator
    DBoundingBox& operator=(const DBoundingBox& rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Assignment operator for the base class
    DBoundingBox& operator=(const Base& rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Destructor
    ~DBoundingBox()
    {
    }

    ///Constructor from two positions
    DBoundingBox(const PositionType& minimum, const PositionType& maximum) :
      Base(minimum, maximum)
    {
    }

    //@}

    /**@name Accessors */
    //@{

    /// Enlarges the bounding box such that it contains a position.
    void enlarge(const PositionType& p)
    {
      for (UInt i = 0; i < DIMENSION; ++i)
      {
        if (p[i] < min_[i]) min_[i] = p[i];
        if (p[i] > max_[i]) max_[i] = p[i];
      }
    }

    ///Enlarges the bounding box such that it contains a position specified by two coordinates
    void enlarge(CoordinateType x, CoordinateType y)
    {
      enlarge(PositionType(x, y));
    }

    //}@

    /**@name Predicates */
    //@{

    /// Equality operator
    bool operator==(const DBoundingBox& rhs) const
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
      for (UInt i = 0; i < DIMENSION; ++i)
      {
        if (position[i] < min_[i] || position[i] > max_[i])
        {
          return false;
        }
      }
      return true;
    }

    ///2D-version encloses(x,y) is for convenience only
    bool encloses(CoordinateType x, CoordinateType y) const
    {
      return encloses(PositionType(x, y));
    }

    /**
         Checks whether this bounding box intersects with another bounding box
    */
    bool intersects(const DBoundingBox& bounding_box) const
    {
      for (UInt i = 0; i < DIMENSION; ++i)
      {
        if (bounding_box.min_[i] > max_[i]) return false;

        if (bounding_box.max_[i] <  min_[i]) return false;
      }
      return true;
    }

    /// Test if bounding box is empty
    bool isEmpty() const
    {
      for (UInt i = 0; i != D; i++)
      {
        if (max_[i] <= min_[i])
        {
          return true;
        }
      }
      return false;
    }

    //@}


  };

  /**@brief Print the contents to a stream.

  @relatesalso DBoundingBox
  */
  template <UInt D>
  std::ostream& operator<<(std::ostream& os, const DBoundingBox<D>& bounding_box)
  {
    os << "--DBOUNDINGBOX BEGIN--" << std::endl;
    os << "MIN --> " << bounding_box.minPosition() << std::endl;
    os << "MAX --> " << bounding_box.maxPosition() << std::endl;
    os << "--DBOUNDINGBOX END--" << std::endl;
    return os;
  }

} // namespace OpenMS

