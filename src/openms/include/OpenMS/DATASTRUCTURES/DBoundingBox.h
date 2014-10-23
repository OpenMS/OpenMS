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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DBOUNDINGBOX_H
#define OPENMS_DATASTRUCTURES_DBOUNDINGBOX_H

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
    DBoundingBox(const DBoundingBox & rhs) :
      Base(rhs)
    {
    }

    /// Assignment operator
    DBoundingBox & operator=(const DBoundingBox & rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Assignment operator for the base class
    DBoundingBox & operator=(const Base & rhs)
    {
      Base::operator=(rhs);
      return *this;
    }

    /// Destructor
    ~DBoundingBox()
    {
    }

    ///Constructor from two positions
    DBoundingBox(const PositionType & minimum, const PositionType & maximum) :
      Base(minimum, maximum)
    {
    }

    //@}

    /**@name Accessors */
    //@{

    /// Enlarges the bounding box such that it contains a position.
    void enlarge(const PositionType & p)
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
    bool operator==(const DBoundingBox & rhs) const
    {
      return Base::operator==(rhs);
    }

    /// Equality operator
    bool operator==(const Base & rhs) const
    {
      return Base::operator==(rhs);
    }

    /**
        @brief Checks whether this range contains a certain point.

        @param position The point's position.
        @returns true if point lies inside this area.
    */
    bool encloses(const PositionType & position) const
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
    bool intersects(const DBoundingBox & bounding_box) const
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
  std::ostream & operator<<(std::ostream & os, const DBoundingBox<D> & bounding_box)
  {
    os << "--DBOUNDINGBOX BEGIN--" << std::endl;
    os << "MIN --> " << bounding_box.minPosition() << std::endl;
    os << "MAX --> " << bounding_box.maxPosition() << std::endl;
    os << "--DBOUNDINGBOX END--" << std::endl;
    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_DBOUNDINGBOX_H
