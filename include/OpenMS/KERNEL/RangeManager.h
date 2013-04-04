// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_RANGEMANAGER_H
#define OPENMS_KERNEL_RANGEMANAGER_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

namespace OpenMS
{
  /**
    @brief Handles the managment of a position and intensity range.

    This is needed for all peak and feature container like Spectrum, MSExperiment and FeatureMap.
  */
  template <UInt D>
  class RangeManager
  {
public:
    /// Dimension of the position range
    enum {DIMENSION = D};
    /// Position range type
    typedef DRange<D> PositionRangeType;
    /// Position Type
    typedef DPosition<D> PositionType;
    /// Intensity range type
    typedef DRange<1> IntensityRangeType;

    /// Default constructor
    RangeManager() :
      int_range_(),
      pos_range_()
    {}

    /// Copy constructor
    RangeManager(const RangeManager & rhs) :
      int_range_(rhs.int_range_),
      pos_range_(rhs.pos_range_)
    {}

    /// Destructor
    virtual ~RangeManager()
    {}

    /// Assignment operator
    RangeManager & operator=(const RangeManager & rhs)
    {
      if (this == &rhs) return *this;

      int_range_ = rhs.int_range_;
      pos_range_ = rhs.pos_range_;

      return *this;
    }

    /// Equality operator
    bool operator==(const RangeManager & rhs) const
    {
      return int_range_ == rhs.int_range_ &&
             pos_range_ == rhs.pos_range_;
    }

    /// Equality operator
    bool operator!=(const RangeManager & rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @name Range methods

      @note The range values are not updated automatically. Call updateRanges() to update the
      values!
    */
    ///@{

    /// Returns the minimum position
    const PositionType & getMin() const
    {
      return pos_range_.minPosition();
    }

    /// Returns the maximum position
    const PositionType & getMax() const
    {
      return pos_range_.maxPosition();
    }

    /// Returns the minimum intensity
    DoubleReal getMinInt() const
    {
      return int_range_.minPosition()[0];
    }

    /// Returns the maximum intensity
    DoubleReal getMaxInt() const
    {
      return int_range_.maxPosition()[0];
    }

    /**
      @brief Updates minimum and maximum position/intensity.

      This method is usually implemented by calling clearRanges() and updateRanges_().
    */
    virtual void updateRanges() = 0;

    /// Resets the ranges
    void clearRanges()
    {
      int_range_ = IntensityRangeType::empty;
      pos_range_ = PositionRangeType::empty;
    }

    ///@}
protected:
    /// Intensity range (1-dimensional)
    IntensityRangeType int_range_;
    /// Position range (D-dimensional)
    PositionRangeType pos_range_;

    /// Updates the range using data points in the iterator range.
    template <class PeakIteratorType>
    void updateRanges_(const PeakIteratorType & begin, const PeakIteratorType & end)
    {
      //prevent invalid range by empty container
      if (begin == end)
      {
        return;
      }

      PositionType min = pos_range_.minPosition();
      PositionType max = pos_range_.maxPosition();

      DoubleReal it_min = int_range_.minPosition()[0];
      DoubleReal it_max = int_range_.maxPosition()[0];

      for (PeakIteratorType it = begin; it != end; ++it)
      {
        //update position
        for (UInt i = 0; i < D; ++i)
        {
          DoubleReal tmp = it->getPosition()[i];
          if (tmp < min[i])
          {
            min[i] = tmp;
          }
          if (tmp > max[i])
          {
            max[i] = tmp;
          }
        }

        //update intensity
        DoubleReal tmp = it->getIntensity();
        if (tmp < it_min)
        {
          it_min = tmp;
        }
        if (tmp > it_max)
        {
          it_max = tmp;
        }
      }

      pos_range_.setMin(min);
      pos_range_.setMax(max);

      int_range_.setMinX(it_min);
      int_range_.setMaxX(it_max);
    }

  };
}  // namespace OpenMS

#endif  // OPENMS_KERNEL_DRANGE_H
