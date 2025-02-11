// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <iosfwd>

namespace OpenMS
{
  /**
    @brief A 1-dimensional raw data mobility point or peak. The unit (ms, 1/K_0, etc) is implicit.

    This data structure is intended for continuous mobility data or centroided mobility data.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MobilityPeak1D
  {
  public:
    ///@name Type definitions
    ///@{
    /// Dimension
    enum{ DIMENSION = 1 };
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
    MobilityPeak1D() = default;

    /// construct with position and intensity
    MobilityPeak1D(PositionType a, IntensityType b) : position_(a), intensity_(b)
    {
    }

    /// Copy constructor
    MobilityPeak1D(const MobilityPeak1D& p) = default;

    // Move constructor
    MobilityPeak1D(MobilityPeak1D&&) noexcept = default;
    
    /// Assignment operator
    MobilityPeak1D& operator=(const MobilityPeak1D& rhs) = default;

    /// Move assignment operator
    MobilityPeak1D& operator=(MobilityPeak1D&&) noexcept = default;


    /**
      @brief Destructor

      @note The destructor is non-virtual although many classes are derived from
      MobilityPeak1D.  This is intentional, since otherwise we would "waste"
      space for a vtable pointer in each instance. Normally you should not derive other classes from
      MobilityPeak1D (unless you know what you are doing, of course).
    */
    ~MobilityPeak1D() noexcept = default;

    ///@}

    /**
      @name Accessors
    */
    ///@{
    /// Non-mutable access to the data point intensity (height)
    IntensityType getIntensity() const
    {
      return intensity_;
    }
    /// Mutable access to the data point intensity (height)
    void setIntensity(IntensityType intensity)
    {
      intensity_ = intensity;
    }

    /// Non-mutable access to m/z
    inline CoordinateType getMobility() const
    {
      return position_[0];
    }

    /// Mutable access to mobility
    inline void setMobility(CoordinateType mobility)
    {
      position_[0] = mobility;
    }

    /// Alias for getMobility()
    inline CoordinateType getPos() const
    {
      return position_[0];
    }

    /// Alias for setMobility()
    inline void setPos(CoordinateType pos)
    {
      position_[0] = pos;
    }

    /// Non-mutable access to the position
    inline PositionType const& getPosition() const
    {
      return position_;
    }

    /// Mutable access to the position
    inline PositionType& getPosition()
    {
      return position_;
    }

    /// Mutable access to the position
    inline void setPosition(PositionType const& position)
    {
      position_ = position;
    }

    ///@}

    /// Equality operator
    bool operator==(const MobilityPeak1D& rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return intensity_ == rhs.intensity_ && position_ == rhs.position_;
#pragma clang diagnostic pop
    }

    /// Equality operator
    bool operator!=(const MobilityPeak1D& rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @name Comparator classes. These classes implement binary predicates that can be used to
      compare two peaks with respect to their intensities, positions.
    */
    ///@{
    /// Comparator by intensity
    struct IntensityLess {
      inline bool operator()(MobilityPeak1D const& left, MobilityPeak1D const& right) const
      {
        return left.getIntensity() < right.getIntensity();
      }

      inline bool operator()(MobilityPeak1D const& left, IntensityType right) const
      {
        return left.getIntensity() < right;
      }

      inline bool operator()(IntensityType left, MobilityPeak1D const& right) const
      {
        return left < right.getIntensity();
      }

      inline bool operator()(IntensityType left, IntensityType right) const
      {
        return left < right;
      }
    };

    /// Comparator by mobility position.
    struct MobilityLess {
      inline bool operator()(const MobilityPeak1D& left, const MobilityPeak1D& right) const
      {
        return left.getMobility() < right.getMobility();
      }

      inline bool operator()(MobilityPeak1D const& left, CoordinateType right) const
      {
        return left.getMobility() < right;
      }

      inline bool operator()(CoordinateType left, MobilityPeak1D const& right) const
      {
        return left < right.getMobility();
      }

      inline bool operator()(CoordinateType left, CoordinateType right) const
      {
        return left < right;
      }
    };

    /// Comparator by position. As this class has dimension 1, this is basically an alias for MobilityLess.
    struct PositionLess {
      inline bool operator()(const MobilityPeak1D& left, const MobilityPeak1D& right) const
      {
        return left.getPosition() < right.getPosition();
      }

      inline bool operator()(const MobilityPeak1D& left, const PositionType& right) const
      {
        return left.getPosition() < right;
      }

      inline bool operator()(const PositionType& left, const MobilityPeak1D& right) const
      {
        return left < right.getPosition();
      }

      inline bool operator()(const PositionType& left, const PositionType& right) const
      {
        return left < right;
      }
    };
    ///@}

  protected:
    /// The data point position
    PositionType position_{};
    /// The data point intensity
    IntensityType intensity_ = 0.0;
  };

  /// Print the contents to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const MobilityPeak1D& point);

} // namespace OpenMS
