// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CHROMATOGRAMPEAK_H
#define OPENMS_KERNEL_CHROMATOGRAMPEAK_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

  /**
    @brief A 1-dimensional raw data point or peak for chromatograms.

    This datastructure is intended for chromatograms.
    If you want to annotated single peaks with meta data, use RichChromatogramPeak instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI ChromatogramPeak
  {
public:
    /**
      @name Type definitions
    */
    //@{
    /// Dimension
    enum {DIMENSION = 1};
    /// Intensity type
    typedef double IntensityType;
    /// Position type
    typedef DPosition<1> PositionType;
    /// Coordinate type
    typedef double CoordinateType;
    //@}

    /**
      @name Constructors and Destructor
     */
    //@{
    /// Default constructor
    inline ChromatogramPeak() :
      position_(),
      intensity_(0)
    {}

    /// Copy constructor
    inline ChromatogramPeak(const ChromatogramPeak & p) :
      position_(p.position_),
      intensity_(p.intensity_)
    {}

    /// Constructor with position and intensity
    inline ChromatogramPeak(const PositionType retention_time, const IntensityType intensity) :
      position_(retention_time),
      intensity_(intensity)
    {}

    /**
      @brief Destructor

      @note The destructor is non-virtual although many classes are derived from
      ChromatogramPeak.  This is intentional, since otherwise we would "waste"
      space for a vtable pointer in each instance. Normally you should not derive other
      classes from ChromatogramPeak (unless you know what you are doing, of course).
    */
    ~ChromatogramPeak()
    {}
    //@}

    /**
      @name Accessors
    */
    //@{

    /// Non-mutable access to the data point intensity (height)
    inline IntensityType getIntensity() const { return intensity_; }
    /// Mutable access to the data point intensity (height)
    inline void setIntensity(IntensityType intensity) { intensity_ = intensity; }

    /// Non-mutable access to RT
    inline CoordinateType getRT() const
    {
      return position_[0];
    }

    /// Mutable access to RT
    inline void setRT(CoordinateType rt)
    {
      position_[0] = rt;
    }

    /// Alias for getRT()
    inline CoordinateType getPos() const
    {
      return position_[0];
    }

    /// Alias for setRT()
    inline void setPos(CoordinateType pos)
    {
      position_[0] = pos;
    }

    /// Alias for getRT()
    inline CoordinateType getMZ() const
    {
      return position_[0];
    }

    /// Alias for setRT()
    inline void setMZ(CoordinateType rt)
    {
      position_[0] = rt;
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

    //@}

    /// Assignment operator
    inline ChromatogramPeak & operator=(const ChromatogramPeak & rhs)
    {
      if (this == &rhs) return *this;

      intensity_ = rhs.intensity_;
      position_ = rhs.position_;

      return *this;
    }

    /// Equality operator
    inline bool operator==(const ChromatogramPeak & rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return intensity_ == rhs.intensity_ && position_ == rhs.position_;
#pragma clang diagnostic pop
    }

    /// Equality operator
    inline bool operator!=(const ChromatogramPeak & rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @name	Comparator classes.

      These classes implement binary predicates that can be used
      to compare two peaks with respect to their intensities, positions.
    */
    //@{
    /// Comparator by intensity
    struct IntensityLess :
      std::binary_function<ChromatogramPeak, ChromatogramPeak, bool>
    {
      inline bool operator()(ChromatogramPeak const & left, ChromatogramPeak const & right) const
      {
        return left.getIntensity() < right.getIntensity();
      }

      inline bool operator()(ChromatogramPeak const & left, IntensityType right) const
      {
        return left.getIntensity() < right;
      }

      inline bool operator()(IntensityType left, ChromatogramPeak const & right) const
      {
        return left < right.getIntensity();
      }

      inline bool operator()(IntensityType left, IntensityType right) const
      {
        return left < right;
      }

    };

    /// Comparator by RT position.
    struct RTLess :
      public std::binary_function<ChromatogramPeak, ChromatogramPeak, bool>
    {
      inline bool operator()(const ChromatogramPeak & left, const ChromatogramPeak & right) const
      {
        return left.getRT() < right.getPos();
      }

      inline bool operator()(ChromatogramPeak const & left, CoordinateType right) const
      {
        return left.getRT() < right;
      }

      inline bool operator()(CoordinateType left, ChromatogramPeak const & right) const
      {
        return left < right.getRT();
      }

      inline bool operator()(CoordinateType left, CoordinateType right) const
      {
        return left < right;
      }

    };

    /// Comparator by position. As this class has dimension 1, this is basically an alias for RTLess.
    struct PositionLess :
      public std::binary_function<ChromatogramPeak, ChromatogramPeak, bool>
    {
      inline bool operator()(const ChromatogramPeak & left, const ChromatogramPeak & right) const
      {
        return left.getPosition() < right.getPosition();
      }

      inline bool operator()(const ChromatogramPeak & left, const PositionType & right) const
      {
        return left.getPosition() < right;
      }

      inline bool operator()(const PositionType & left, const ChromatogramPeak & right) const
      {
        return left < right.getPosition();
      }

      inline bool operator()(const PositionType & left, const PositionType & right) const
      {
        return left < right;
      }
    };
    //@}
protected:
    /// The data point position
    PositionType    position_;
    /// The data point intensity
    IntensityType intensity_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ChromatogramPeak & point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CHROMATOGRAMPEAK_H
