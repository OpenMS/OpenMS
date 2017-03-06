// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_CUMULATIVEHISTOGRAM_H
#define OPENMS_MATH_STATISTICS_CUMULATIVEHISTOGRAM_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>

//STL
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>

namespace OpenMS
{
  namespace Math
  {

    /**
        @brief Representation of a histogram, which counts the cumulative number of observations.

        The first template argument gives the Type of the
        values that are stored in the bins. The second argument
        gives the type for the bin size and range.

        @ingroup Math
    */
    template <typename ValueType = UInt, typename BinSizeType = double>
    class CumulativeHistogram :
      public Histogram<ValueType, BinSizeType>
    {
    public:

      ///default constructor
      CumulativeHistogram() :
        Histogram< ValueType, BinSizeType >(),
        complementary_(false),
        inclusive_(true)
      {
      }


      ///Copy constructor
      CumulativeHistogram(const CumulativeHistogram & cum_histogram) :
        Histogram< ValueType, BinSizeType>(cum_histogram),
        complementary_(cum_histogram.complementary_),
        inclusive_(cum_histogram.inclusive_)
      {
      }


      /**
        @brief constructor with min, max and bin size

        @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
      */
      CumulativeHistogram(BinSizeType min, BinSizeType max, BinSizeType bin_size):
        Histogram< ValueType, BinSizeType >(min, max, bin_size),
        complementary_(false),
        inclusive_(true)
      {
      }


      /**
        @brief constructor with min, max, bin size, complementary, and inclusive

        @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
      */
      CumulativeHistogram(BinSizeType min, BinSizeType max, BinSizeType bin_size, bool complementary, bool inclusive):
        Histogram< ValueType, BinSizeType >(min, max, bin_size),
        complementary_(complementary),
        inclusive_(inclusive)
      {
      }


      /**
        @brief constructor with data iterator and min, max, bin_size parameters

        @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
      */
      template <typename DataIterator>
      CumulativeHistogram(DataIterator begin, DataIterator end, BinSizeType min, BinSizeType max, BinSizeType bin_size) :
        Histogram< ValueType, BinSizeType >(begin, end, min, max, bin_size),
        complementary_(false),
        inclusive_(true)
      {
      }


      /**
        @brief constructor with data iterator and min, max, bin_size, and complementary parameters

        @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
      */
      template <typename DataIterator>
      CumulativeHistogram(DataIterator begin, DataIterator end, BinSizeType min, BinSizeType max, BinSizeType bin_size,
                          bool complementary, bool inclusive) :
        Histogram< ValueType, BinSizeType >(min, max, bin_size),
        complementary_(complementary),
        inclusive_(inclusive)
      {
        for (DataIterator it = begin; it != end; ++it)
        {
          this->inc((BinSizeType) *it);
        }
      }


      ///returns the highest value of all bins
      ValueType maxValue() const
      {
        return this->complementary_ ? this->bins_.front() : this->bins_.back();
      }

      ///returns the lowest value of all bins
      ValueType minValue() const
      {
        return this->complementary_ ? this->bins_.back() : this->bins_.front();
      }


      /**
        @brief increases all the bins up to value @p val (or from value @p val if complementary) by @p increment

          @return The index of the bin the given value @p val belongs to.

        @exception Exception::OutOfRange is thrown if the value is out of valid range
      */
      Size inc(BinSizeType val, ValueType increment = 1)
      {
        Size bin_index = this->valToBin_(val);
        if (this->complementary_)
        {
          for (size_t i = 0; i < bin_index; ++i)
          {
           this->bins_[i] += increment;
          }
        }
        else
        {
          for (size_t i = bin_index + 1; i < this->bins_.size(); ++i)
          {
            this->bins_[i] += increment;
          }
        }
        if (this->inclusive_)
        {
          this->bins_[bin_index] += increment;
        }
        return bin_index;
      }

    protected:

      bool complementary_;  // Whether the complementary cumulative histogram should be computed
      bool inclusive_;  // Whether the bin corresponding to a particular value should be increased when the value is encountered
    };
  }   // namespace Math
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_HISTOGRAM_H
