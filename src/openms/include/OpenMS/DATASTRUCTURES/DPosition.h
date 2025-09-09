// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <array>
#include <cmath>  // for std::abs on integrals and floats
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
  using CoordinateType =  TCoordinateType;

  using DataType = std::array<CoordinateType, D>;

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
  DPosition() = default;

  /// Constructor that fills all dimensions with the value @p x
  DPosition(CoordinateType x);

  /// Constructor only for DPosition<2> that takes two Coordinates.
  DPosition(CoordinateType x, CoordinateType y);

  /// Constructor only for DPosition<3> that takes three Coordinates.
  DPosition(CoordinateType x, CoordinateType y, CoordinateType z);

  /// Copy constructor
  DPosition(const DPosition& pos) = default;

  /// Move constructor
  DPosition(DPosition&& rhs) noexcept = default;

  /// Assignment operator
  DPosition& operator=(const DPosition& source) = default;

  /// Move Assignment operator
  DPosition& operator=(DPosition&& source) noexcept = default;

  /// Destructor (not-virtual as this will save a lot of space!)
  ~DPosition() noexcept = default;

  //@}

  /// Swap the two objects
  void swap(DPosition& rhs) noexcept
  {
    std::swap(coordinate_, rhs.coordinate_);
  }

  /// Make all dimension values positive
  DPosition& abs() noexcept
  {
    for (Size i = 0; i < D; ++i)
    {
      coordinate_[i] = std::abs(coordinate_[i]);
    }
    return *this;
  }

  /**@name Accessors */
  //@{

  /// Const accessor for the dimensions
  CoordinateType operator[](Size index) const;

  /// Accessor for the dimensions
  CoordinateType& operator[](Size index);

  /// Name accessor for the first dimension. Only for DPosition<2>, for visualization.
  CoordinateType getX() const;

  /// Name accessor for the second dimension. Only for DPosition<2>, for visualization.
  CoordinateType getY() const;

  /// Name mutator for the first dimension. Only for DPosition<2>, for visualization.
  void setX(CoordinateType c);

  /// Name mutator for the second dimension. Only for DPosition<2>, for visualization.
  void setY(CoordinateType c);

  /// Equality operator
  bool operator==(const DPosition& point) const = default;

  /// Equality operator
  bool operator!=(const DPosition& point) const = default;
  /**
    @brief Lexicographical less than operator.
    Lexicographical comparison from dimension 0 to dimension D-1 is done.
  */
  bool operator<(const DPosition& point) const;

  /// Lexicographical greater less or equal operator.
  bool operator<=(const DPosition& point) const;

  /// Spatially (geometrically) less or equal operator. All coordinates must be "<=".
  bool spatiallyLessEqual(const DPosition& point) const;

  /// Spatially (geometrically) greater or equal operator. All coordinates must be ">=".
  bool spatiallyGreaterEqual(const DPosition& point) const;

  /// Lexicographical greater than operator.
  bool operator>(const DPosition& point) const;

  /// Lexicographical greater or equal operator.
  bool operator>=(const DPosition& point) const;

    /// Addition (a bit inefficient)
    DPosition operator+(const DPosition& point) const;

    /// Addition
    DPosition& operator+=(const DPosition& point);

    /// Subtraction (a bit inefficient)
    DPosition operator-(const DPosition& point) const;

    /// Subtraction
    DPosition& operator-=(const DPosition& point);

    /// Negation (a bit inefficient)
    DPosition operator-() const;

    /// Inner product
    CoordinateType operator*(const DPosition& point) const;

    /// Scalar multiplication
    DPosition& operator*=(CoordinateType scalar);

    /// Scalar division
    DPosition& operator/=(CoordinateType scalar);

    /// Returns the number of dimensions
    constexpr static Size size()
    {
      return D;
    }

    /// Set all dimensions to zero
    void clear()
    {
      coordinate_.fill(0);
    }

    //@}

    /** @name Static values */
    //@{
    /// all zero
    inline static constexpr DPosition zero()
    {
      return DPosition(0);
    }

    /// smallest positive
    inline static constexpr DPosition minPositive()
    {
      return DPosition((std::numeric_limits<typename DPosition::CoordinateType>::min)());
    }

    /// smallest negative
    inline static constexpr DPosition minNegative()
    {
      return DPosition(std::numeric_limits<typename DPosition::CoordinateType>::lowest());
    }

    /// largest positive
    inline static constexpr DPosition maxPositive()
    {
      return DPosition((std::numeric_limits<typename DPosition::CoordinateType>::max)());
    }

    //@}

    /** @name Iteration */
    //@{
    /// Non-mutable begin iterator
    const_iterator begin() const
    {
      return &coordinate_[0];
    }

    /// Non-mutable end iterator
    const_iterator end() const
    {
      return &coordinate_[0] + D;
    }

    /// Mutable begin iterator
    iterator begin()
    {
      return &coordinate_[0];
    }

    /// Mutable end iterator
    iterator end()
    {
      return &coordinate_[0] + D;
    }

    //@}

protected:
    DataType coordinate_{};
  }; // DPosition

  // Member template definitions for DPosition

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>::DPosition(typename DPosition<D, TCoordinateType>::CoordinateType x)
  {
    std::fill(coordinate_.begin(), coordinate_.end(), x);
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>::DPosition(typename DPosition<D, TCoordinateType>::CoordinateType x, typename DPosition<D, TCoordinateType>::CoordinateType y)
  {
    static_assert(D == 2, "DPosition<D, TCoordinateType>:DPosition(x,y): index overflow!");
    coordinate_[0] = x;
    coordinate_[1] = y;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>::DPosition(typename DPosition<D, TCoordinateType>::CoordinateType x, typename DPosition<D, TCoordinateType>::CoordinateType y, typename DPosition<D, TCoordinateType>::CoordinateType z)
  {
    static_assert(D == 3, "DPosition<D, TCoordinateType>:DPosition(x,y,z): index overflow!");
    coordinate_[0] = x;
    coordinate_[1] = y;
    coordinate_[2] = z;
  }

  template <UInt D, typename TCoordinateType>
  typename DPosition<D, TCoordinateType>::CoordinateType DPosition<D, TCoordinateType>::operator[](Size index) const
  {
    OPENMS_PRECONDITION(index < D, "DPosition<D,TCoordinateType>:operator [] (Position): index overflow!");
    return coordinate_[index];
  }

  template <UInt D, typename TCoordinateType>
  typename DPosition<D, TCoordinateType>::CoordinateType& DPosition<D, TCoordinateType>::operator[](Size index)
  {
    OPENMS_PRECONDITION(index < D, "DPosition<D,TCoordinateType>:operator [] (Position): index overflow!");
    return coordinate_[index];
  }

  template <UInt D, typename TCoordinateType>
  typename DPosition<D, TCoordinateType>::CoordinateType DPosition<D, TCoordinateType>::getX() const
  {
    OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:getX(): index overflow!");
    return coordinate_[0];
  }

  template <UInt D, typename TCoordinateType>
  typename DPosition<D, TCoordinateType>::CoordinateType DPosition<D, TCoordinateType>::getY() const
  {
    OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:getY(): index overflow!");
    return coordinate_[1];
  }

  template <UInt D, typename TCoordinateType>
  void DPosition<D, TCoordinateType>::setX(typename DPosition<D, TCoordinateType>::CoordinateType c)
  {
    OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:setX(): index overflow!");
    coordinate_[0] = c;
  }

  template <UInt D, typename TCoordinateType>
  void DPosition<D, TCoordinateType>::setY(typename DPosition<D, TCoordinateType>::CoordinateType c)
  {
    OPENMS_PRECONDITION(D == 2, "DPosition<D,TCoordinateType>:setY(): index overflow!");
    coordinate_[1] = c;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::operator<(const DPosition& point) const
  {
    return coordinate_ < point.coordinate_;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::operator<=(const DPosition& point) const
  {
    return coordinate_ <= point.coordinate_;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::spatiallyLessEqual(const DPosition& point) const
  {
    for (Size i = 0; i < D; i++)
    {
      if (coordinate_[i] > point.coordinate_[i]) return false;
    }
    return true;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::spatiallyGreaterEqual(const DPosition& point) const
  {
    for (Size i = 0; i < D; i++)
    {
      if (coordinate_[i] < point.coordinate_[i]) return false;
    }
    return true;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::operator>(const DPosition& point) const
  {
    return coordinate_ > point.coordinate_;
  }

  template <UInt D, typename TCoordinateType>
  bool DPosition<D, TCoordinateType>::operator>=(const DPosition& point) const
  {
    return coordinate_ >= point.coordinate_;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> DPosition<D, TCoordinateType>::operator+(const DPosition& point) const
  {
    DPosition result(*this);
    for (Size i = 0; i < D; ++i)
    {
      result.coordinate_[i] += point.coordinate_[i];
    }
    return result;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>& DPosition<D, TCoordinateType>::operator+=(const DPosition& point)
  {
    for (Size i = 0; i < D; ++i)
    {
      coordinate_[i] += point.coordinate_[i];
    }
    return *this;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> DPosition<D, TCoordinateType>::operator-(const DPosition& point) const
  {
    DPosition result(*this);
    for (Size i = 0; i < D; ++i)
    {
      result.coordinate_[i] -= point.coordinate_[i];
    }
    return result;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>& DPosition<D, TCoordinateType>::operator-=(const DPosition& point)
  {
    for (Size i = 0; i < D; ++i)
    {
      coordinate_[i] -= point.coordinate_[i];
    }
    return *this;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> DPosition<D, TCoordinateType>::operator-() const
  {
    DPosition result(*this);
    for (Size i = 0; i < D; ++i)
    {
      result.coordinate_[i] = -result.coordinate_[i];
    }
    return result;
  }

  template <UInt D, typename TCoordinateType>
  typename DPosition<D, TCoordinateType>::CoordinateType DPosition<D, TCoordinateType>::operator*(const DPosition& point) const
  {
    CoordinateType prod(0);
    for (Size i = 0; i < D; ++i)
    {
      prod += (point.coordinate_[i] * coordinate_[i]);
    }
    return prod;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>& DPosition<D, TCoordinateType>::operator*=(typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      coordinate_[i] *= scalar;
    }
    return *this;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType>& DPosition<D, TCoordinateType>::operator/=(typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      coordinate_[i] /= scalar;
    }
    return *this;
  }

  // Free function template definitions

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator*(DPosition<D, TCoordinateType> position, typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] *= scalar;
    }
    return position;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator*(typename DPosition<D, TCoordinateType>::CoordinateType scalar, DPosition<D, TCoordinateType> position)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] *= scalar;
    }
    return position;
  }

  template <UInt D, typename TCoordinateType>
  DPosition<D, TCoordinateType> operator/(DPosition<D, TCoordinateType> position, typename DPosition<D, TCoordinateType>::CoordinateType scalar)
  {
    for (Size i = 0; i < D; ++i)
    {
      position[i] /= scalar;
    }
    return position;
  }

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

