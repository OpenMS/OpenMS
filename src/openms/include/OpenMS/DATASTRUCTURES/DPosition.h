// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <algorithm>
#include <limits>
#include <ostream>

namespace OpenMS
{
  /**
    @brief Representation of a coordinate in D-dimensional space.

    @ingroup Datastructures
  */
  template <UInt D, typename TCoordinateType = double>
  class DPosition
  {
public:

    /// Coordinate type
    typedef TCoordinateType CoordinateType;
    /// Mutable iterator
    typedef CoordinateType* Iterator;
    /// Non-mutable iterator
    typedef const CoordinateType* ConstIterator;
    /// Dimensions
    enum
    {
      DIMENSION = D
    };
    /**
      @name STL compatibility type definitions
    */
    //@{
    typedef CoordinateType value_type;
    typedef CoordinateType& reference;
    typedef CoordinateType* pointer;
    typedef CoordinateType* iterator;
    typedef const CoordinateType* const_iterator;
    //@}

    /**
      @name Constructors and Destructor
    */
    //@{
    /**
      @brief Default constructor.

      Creates a position with all coordinates zero.
    */
    DPosition()
    {
      clear();
    }

    /// Destructor (not-virtual as this will save a lot of space!)
    ~DPosition()
    {
    }

    /// Constructor that fills all dimensions with the value @p x
    DPosition(CoordinateType x)
    {
      std::fill(&(coordinate_[0]), &(coordinate_[D]), x);
    }

    /// Copy constructor
    DPosition(const DPosition& pos)
    {
      std::copy(&(pos.coordinate_[0]), &(pos.coordinate_[D]),
                &(coordinate_[0]));
    }

    /// Move constructor
    DPosition(DPosition&& rhs) noexcept
    {
      // NOTE: do not change this before testing with nightly Windows builds ( = default causes segfault)
      std::move(std::begin(rhs.coordinate_), std::end(rhs.coordinate_), &coordinate_[0]);
    }

    /// Constructor only for DPosition<2> that takes two Coordinates.
    DPosition(CoordinateType x, CoordinateType y)
    {
      OPENMS_PRECONDITION(D == 2, "DPosition<D, TCoordinateType>:DPosition(x,y): index overflow!");
      coordinate_[0] = x;
      coordinate_[1] = y;
    }

    /// Assignment operator
    DPosition& operator=(const DPosition& source)
    {
      if (&source == this) return *this;

      std::copy(&(source.coordinate_[0]), &(source.coordinate_[D]),
                &(coordinate_[0]));

      return *this;
    }

    //@}

    /**@name Accessors */
    //@{

    ///Const accessor for the dimensions
    CoordinateType operator[](Size index) const
    {
      OPENMS_PRECONDITION(index < D, "DPosition<D,TCoordinateType>:operator [] (Position): index overflow!");
      return coordinate_[index];
    }

    ///Accessor for the dimensions
    CoordinateType& operator[](Size index)
    {
      OPENMS_PRECONDITION(index < D, "DPosition<D,TCoordinateType>:operator [] (Position): index overflow!");
      return coordinate_[index];
    }

    ///Name accessor for the first dimension. Only for DPosition<2>, for visualization.
    CoordinateType getX() const
    {
      OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:getX(): index overflow!");
      return coordinate_[0];
    }

    ///Name accessor for the second dimension. Only for DPosition<2>, for visualization.
    CoordinateType getY() const
    {
      OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:getY(): index overflow!");
      return coordinate_[1];
    }

    ///Name mutator for the first dimension. Only for DPosition<2>, for visualization.
    void setX(CoordinateType c)
    {
      OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:setX(): index overflow!");
      coordinate_[0] = c;
    }

    ///Name mutator for the second dimension. Only for DPosition<2>, for visualization.
    void setY(CoordinateType c)
    {
      OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:setY(): index overflow!");
      coordinate_[1] = c;
    }

    /// Equality operator
    bool operator==(const DPosition& point) const
    {
      for (Size i = 0; i < D; i++)
      {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
        if (coordinate_[i] != point.coordinate_[i]) return false;

#pragma clang diagnostic pop
      }
      return true;
    }

    /// Equality operator
    bool operator!=(const DPosition& point) const
    {
      return !(operator==(point));
    }

    /**
      @brief Lexicographical less than operator.
      Lexicographical comparison from dimension 0 to dimension D-1 is done.
    */
    bool operator<(const DPosition& point) const
    {
      for (Size i = 0; i < D; i++)
      {
        if (coordinate_[i] < point.coordinate_[i]) return true;

        if (coordinate_[i] > point.coordinate_[i]) return false;
      }
      return false;
    }

    /// Lexicographical greater less or equal operator.
    bool operator<=(const DPosition& point) const
    {
      for (Size i = 0; i < D; i++)
      {
        if (coordinate_[i] < point.coordinate_[i]) return true;

        if (coordinate_[i] > point.coordinate_[i]) return false;
      }
      return true;
    }

    /// Spatially (geometrically) less or equal operator. All coordinates must be "<=".
    bool spatiallyLessEqual(const DPosition& point) const
    {
      for (Size i = 0; i < D; i++)
      {
        if (coordinate_[i] > point.coordinate_[i]) return false;
      }
      return true;
    }

    /// Spatially (geometrically) greater or equal operator. All coordinates must be ">=".
    bool spatiallyGreaterEqual(const DPosition& point) const
    {
      for (Size i = 0; i < D; i++)
      {
        if (coordinate_[i] < point.coordinate_[i]) return false;
      }
      return true;
    }

    /// Lexicographical greater than operator.
    bool operator>(const DPosition& point) const
    {
      return !(operator<=(point));
    }

    /// Lexicographical greater or equal operator.
    bool operator>=(const DPosition& point) const
    {
      return !operator<(point);
    }

    /// Addition (a bit inefficient)
    DPosition operator+(const DPosition& point) const
    {
      DPosition result(*this);
      for (Size i = 0; i < D; ++i)
      {
        result.coordinate_[i] += point.coordinate_[i];
      }
      return result;
    }

    /// Addition
    DPosition& operator+=(const DPosition& point)
    {
      for (Size i = 0; i < D; ++i)
      {
        coordinate_[i] += point.coordinate_[i];
      }
      return *this;
    }

    /// Subtraction (a bit inefficient)
    DPosition operator-(const DPosition& point) const
    {
      DPosition result(*this);
      for (Size i = 0; i < D; ++i)
      {
        result.coordinate_[i] -= point.coordinate_[i];
      }
      return result;
    }

    /// Subtraction
    DPosition& operator-=(const DPosition& point)
    {
      for (Size i = 0; i < D; ++i)
      {
        coordinate_[i] -= point.coordinate_[i];
      }
      return *this;
    }

    /// Negation (a bit inefficient)
    DPosition   operator-() const
    {
      DPosition<D, CoordinateType> result(*this);
      for (Size i = 0; i < D; ++i)
      {
        result.coordinate_[i] = -result.coordinate_[i];
      }
      return result;
    }

    /// Inner product
    CoordinateType operator*(const DPosition& point) const
    {
      CoordinateType prod(0);
      for (Size i = 0; i < D; ++i)
      {
        prod += (point.coordinate_[i] * coordinate_[i]);
      }
      return prod;
    }

    /// Scalar multiplication
    DPosition& operator*=(CoordinateType scalar)
    {
      for (Size i = 0; i < D; ++i)
      {
        coordinate_[i] *= scalar;
      }
      return *this;
    }

    /// Scalar division
    DPosition& operator/=(CoordinateType scalar)
    {
      for (Size i = 0; i < D; ++i)
      {
        coordinate_[i] /= scalar;
      }
      return *this;
    }

    /// Returns the number of dimensions
    static Size size()
    {
      return D;
    }

    /// Set all dimensions to zero
    void clear()
    {
      for (Size i = 0; i < D; ++i)
      {
        coordinate_[i] = static_cast<CoordinateType>(0);
      }
    }

    //@}

    /** @name Static values */
    //@{
    /// all zero
    inline static const DPosition zero()
    {
      return DPosition(0);
    }

    /// smallest positive
    inline static const DPosition minPositive()
    {
      return DPosition((std::numeric_limits<typename DPosition::CoordinateType>::min)());
    }

    /// smallest negative
    inline static const DPosition minNegative()
    {
      return DPosition(-(std::numeric_limits<typename DPosition::CoordinateType>::max)());
    }

    /// largest positive
    inline static const DPosition maxPositive()
    {
      return DPosition((std::numeric_limits<typename DPosition::CoordinateType>::max)());
    }

    //@}

    /** @name Iteration */
    //@{
    /// Non-mutable begin iterator
    ConstIterator begin() const
    {
      return &(coordinate_[0]);
    }

    /// Non-mutable end iterator
    ConstIterator end() const
    {
      return &(coordinate_[0]) + D;
    }

    /// Mutable begin iterator
    Iterator begin()
    {
      return &(coordinate_[0]);
    }

    /// Mutable end iterator
    Iterator end()
    {
      return &(coordinate_[0]) + D;
    }

    //@}

protected:
    CoordinateType coordinate_[D];

  }; // DPosition

  /// Scalar multiplication (a bit inefficient)
  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator*(DPosition<D, TCoordinateType> position, typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] *= scalar;
    }
    return position;
  }

  /// Scalar multiplication (a bit inefficient)
  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator*(typename DPosition<D, TCoordinateType>::CoordinateType scalar, DPosition<D, TCoordinateType> position)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] *= scalar;
    }
    return position;
  }

  /// Scalar multiplication (a bit inefficient)
  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator/(DPosition<D, TCoordinateType> position, typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] /= scalar;
    }
    return position;
  }

  /// Print the contents to a stream.
  template <UInt D, typename TCoordinateType>
  std::ostream& operator<<(std::ostream& os, const DPosition<D, TCoordinateType>& pos)
  {
    os << precisionWrapper(pos[0]);
    for (UInt i = 1; i < D; ++i)
    {
      os << ' ' << precisionWrapper(pos[i]);
    }
    return os;
  }

} // namespace OpenMS

