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
#include <functional>

namespace OpenMS
{

  /**
    @brief A 2-dimensional raw data point or peak.

    This data structure is intended for continuous data or peak data.
    If you want to annotated single peaks with meta data, use RichPeak2D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI Peak2D
  {
public:

    ///@name Type definitions
    ///@{

    /// Intensity type
    typedef float IntensityType;
    /// Coordinate type (of the position)
    typedef double CoordinateType;
    /// Position type
    typedef DPosition<2> PositionType;
    ///@}

    /// @name Dimension descriptions
    ///@{

    /// This enum maps the symbolic names of the dimensions to numbers
    enum DimensionDescription
    {
      RT = 0,       ///< Retention time dimension id (0 if used as a const int)
      MZ = 1,       ///< Mass-to-charge dimension id (1 if used as a const int)
      DIMENSION = 2 ///< Number of dimensions
    };

    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionName(UInt const dim);
    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionNameRT();
    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionNameMZ();

    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionName(UInt const dim);
    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionNameRT();
    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionNameMZ();

    /// Unit of measurement (abbreviated form)
    static char const * shortDimensionUnit(UInt const dim);
    /// Unit of measurement (abbreviated form)
    static char const * shortDimensionUnitRT();
    /// Unit of measurement (abbreviated form)
    static char const * shortDimensionUnitMZ();

    /// Unit of measurement (self-explanatory form)
    static char const * fullDimensionUnit(UInt const dim);
    /// Unit of measurement (self-explanatory form)
    static char const * fullDimensionUnitRT();
    /// Unit of measurement (self-explanatory form)
    static char const * fullDimensionUnitMZ();

    ///@}

protected:

    /// @name Dimension descriptions
    ///@{

    /// Short name of the dimension (abbreviated form)
    static char const * const dimension_name_short_[DIMENSION];

    /// Full name of the dimension (self-explanatory form)
    static char const * const dimension_name_full_[DIMENSION];

    /// Unit of measurement (abbreviated form)
    static char const * const dimension_unit_short_[DIMENSION];

    /// Unit of measurement (self-explanatory form)
    static char const * const dimension_unit_full_[DIMENSION];

    ///@}

public:

    ///@name Constructors and Destructor
    ///@{
    /// Default constructor
    Peak2D() = default;

    /// Member constructor
    explicit Peak2D(const PositionType& pos, const IntensityType in) :
      position_(pos),
      intensity_(in)
    {}

    /// Copy constructor
    Peak2D(const Peak2D & p) = default;

    /// Move constructor
    Peak2D(Peak2D&&) noexcept = default;

    /// Assignment operator
    Peak2D& operator=(const Peak2D& rhs) = default;

    /// Move assignment operator
    Peak2D& operator=(Peak2D&&) noexcept = default;
    /**
      @brief Destructor

      @note The destructor is non-virtual although many classes are derived from
      Peak2D.  This is intentional, since otherwise we would "waste"
      space for a vtable pointer in each instance. Normally you should not derive other classes from
      Peak2D (unless you know what you are doing, of course).
    */
    ~Peak2D() noexcept = default;
    
    ///@}

    ///@name Accessors
    ///@{
    /// Non-mutable access to the data point intensity (height)
    IntensityType getIntensity() const
    {
      return intensity_;
    }

    /// Sets data point intensity (height)
    void setIntensity(IntensityType intensity)
    {
      intensity_ = intensity;
    }

    /// Non-mutable access to the position
    PositionType const & getPosition() const
    {
      return position_;
    }

    /// Mutable access to the position
    PositionType & getPosition()
    {
      return position_;
    }

    /// Mutable access to the position
    void setPosition(const PositionType & position)
    {
      position_ = position;
    }

    /// Returns the m/z coordinate (index 1)
    CoordinateType getMZ() const
    {
      return position_[MZ];
    }

    /// Mutable access to the m/z coordinate (index 1)
    void setMZ(CoordinateType coordinate)
    {
      position_[MZ] = coordinate;
    }

    /// Returns the RT coordinate (index 0)
    CoordinateType getRT() const
    {
      return position_[RT];
    }

    /// Mutable access to the RT coordinate (index 0)
    void setRT(CoordinateType coordinate)
    {
      position_[RT] = coordinate;
    }

    ///@}

    /// Equality operator
    bool operator==(const Peak2D & rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return intensity_ == rhs.intensity_ && position_ == rhs.position_;
#pragma clang diagnostic pop
    }

    /// Equality operator
    bool operator!=(const Peak2D & rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @name	Comparator classes.
      These classes implement binary predicates that can be used
      to compare two peaks with respect to their intensities, positions, etc.
    */
    ///@{
    /// Comparator by intensity
    struct IntensityLess
    {
      bool operator()(const Peak2D & left, const Peak2D & right) const
      {
        return left.getIntensity() < right.getIntensity();
      }

      bool operator()(const Peak2D & left, IntensityType right) const
      {
        return left.getIntensity() < right;
      }

      bool operator()(IntensityType left, const Peak2D & right) const
      {
        return left < right.getIntensity();
      }

      bool operator()(IntensityType left, IntensityType right) const
      {
        return left < right;
      }

    };

    /// Comparator by RT position
    struct RTLess
    {
      bool operator()(const Peak2D & left, const Peak2D & right) const
      {
        return left.getRT() < right.getRT();
      }

      bool operator()(const Peak2D & left, CoordinateType right) const
      {
        return left.getRT() < right;
      }

      bool operator()(CoordinateType left, const Peak2D & right) const
      {
        return left < right.getRT();
      }

      bool operator()(CoordinateType left, CoordinateType right) const
      {
        return left < right;
      }

    };

    /// Comparator by m/z position
    struct MZLess
    {
      bool operator()(const Peak2D & left, const Peak2D & right) const
      {
        return left.getMZ() < right.getMZ();
      }

      bool operator()(const Peak2D & left, CoordinateType right) const
      {
        return left.getMZ() < right;
      }

      bool operator()(CoordinateType left, const Peak2D & right) const
      {
        return left < right.getMZ();
      }

      bool operator()(CoordinateType left, CoordinateType right) const
      {
        return left < right;
      }

    };

    /// Comparator by position. Lexicographical comparison (first RT then m/z) is done.
    struct PositionLess
    {
      bool operator()(const Peak2D & left, const Peak2D & right) const
      {
        return left.getPosition() < right.getPosition();
      }

      bool operator()(const Peak2D & left, const PositionType & right) const
      {
        return left.getPosition() < right;
      }

      bool operator()(const PositionType & left, const Peak2D & right) const
      {
        return left < right.getPosition();
      }

      bool operator()(const PositionType & left, const PositionType & right) const
      {
        return left < right;
      }

    };
    ///@}

    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Peak2D & point);

protected:

    /// The data point position
    PositionType position_{};
    /// The data point intensity
    IntensityType intensity_{};
  };

  /// Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Peak2D & point);

} // namespace OpenMS

