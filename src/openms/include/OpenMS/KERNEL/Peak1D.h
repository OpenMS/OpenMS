// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <iosfwd>

namespace OpenMS
{

  /**
    @brief A 1-dimensional raw data point or peak.

    This data structure is intended for continuous data or peak data.
    If you want to annotated single peaks with meta data, use RichPeak1D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI Peak1D
  {
public:

    ///@name Type definitions
    ///@{
    /// Dimension
    enum {DIMENSION = 1};
    /// Intensity type
    using IntensityType = float;
    /// Position type
    using PositionType = DPosition<1>;
    /// Coordinate type
    using CoordinateType = double;
    ///@}

    ///@name Constructors and Destructor
    ///@{
    /// Default constructor
    inline Peak1D() = default;

    /// construct with position and intensity
    inline Peak1D(PositionType a, IntensityType b) :
      position_(a),
      intensity_(b)
    {}

    /// Copy constructor
    Peak1D(const Peak1D & p) = default;

    Peak1D(Peak1D&&) noexcept = default;

    /// Assignment operator
    Peak1D& operator=(const Peak1D& rhs) = default;

    /// Move assignment operator
    Peak1D& operator=(Peak1D&&) noexcept = default;

    /**
      @brief Destructor

      @note The destructor is non-virtual although many classes are derived from
      Peak1D.  This is intentional, since otherwise we would "waste"
      space for a vtable pointer in each instance. Normally you should not derive other classes from
      Peak1D (unless you know what you are doing, of course).
    */
    ~Peak1D() = default;

    ///@}

    /**
      @name Accessors
    */
    ///@{
    /// Non-mutable access to the data point intensity (height)
    inline IntensityType getIntensity() const { return intensity_; }
    /// Mutable access to the data point intensity (height)
    inline void setIntensity(IntensityType intensity) { intensity_ = intensity; }

    /// Non-mutable access to m/z
    inline CoordinateType getMZ() const
    {
      return position_[0];
    }

    /// Mutable access to m/z
    inline void setMZ(CoordinateType mz)
    {
      position_[0] = mz;
    }

    /// Alias for getMZ()
    inline CoordinateType getPos() const
    {
      return position_[0];
    }

    /// Alias for setMZ()
    inline void setPos(CoordinateType pos)
    {
      position_[0] = pos;
    }

    /// Non-mutable access to the position
    inline PositionType const & getPosition() const
    {
      return position_;
    }

    /// Mutable access to the position
    inline PositionType & getPosition()
    {
      return position_;
    }

    /// Mutable access to the position
    inline void setPosition(PositionType const & position)
    {
      position_ = position;
    }

    ///@}

    /// Equality operator
    bool operator==(const Peak1D & rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return intensity_ == rhs.intensity_ && position_ == rhs.position_;
#pragma clang diagnostic pop
    }

    /// Equality operator
    bool operator!=(const Peak1D & rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @name Comparator classes. These classes implement binary predicates that can be used to
      compare two peaks with respect to their intensities, positions.
    */
    ///@{
    /// Comparator by intensity
    struct IntensityLess
    {
      inline bool operator()(Peak1D const & left, Peak1D const & right) const
      {
        return left.getIntensity() < right.getIntensity();
      }

      inline bool operator()(Peak1D const & left, IntensityType right) const
      {
        return left.getIntensity() < right;
      }

      inline bool operator()(IntensityType left, Peak1D const & right) const
      {
        return left < right.getIntensity();
      }

      inline bool operator()(IntensityType left, IntensityType right) const
      {
        return left < right;
      }

    };

    /// Comparator by m/z position.
    struct MZLess
    {
      inline bool operator()(const Peak1D & left, const Peak1D & right) const
      {
        return left.getMZ() < right.getMZ();
      }

      inline bool operator()(Peak1D const & left, CoordinateType right) const
      {
        return left.getMZ() < right;
      }

      inline bool operator()(CoordinateType left, Peak1D const & right) const
      {
        return left < right.getMZ();
      }

      inline bool operator()(CoordinateType left, CoordinateType right) const
      {
        return left < right;
      }

    };

    /// Comparator by position. As this class has dimension 1, this is basically an alias for MZLess.
    struct PositionLess
    {
      inline bool operator()(const Peak1D & left, const Peak1D & right) const
      {
        return left.getPosition() < right.getPosition();
      }

      inline bool operator()(const Peak1D & left, const PositionType & right) const
      {
        return left.getPosition() < right;
      }

      inline bool operator()(const PositionType & left, const Peak1D & right) const
      {
        return left < right.getPosition();
      }

      inline bool operator()(const PositionType & left, const PositionType & right) const
      {
        return left < right;
      }

    };
    ///@}

protected:
    /// The data point position
    PositionType  position_;
    /// The data point intensity
    IntensityType intensity_ = 0.0;
  };

  /// Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Peak1D & point);

} // namespace OpenMS

