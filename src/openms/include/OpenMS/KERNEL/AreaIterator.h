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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS includes
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/PeakIndex.h>

// STL includes
#include <iterator>

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief Forward iterator for an area of peaks in an experiment

        This iterator allows us to move through the data structure in a linear
        manner i.e. we don't need to jump to the next spectrum manually.

        @note This iterator iterates over spectra with MS level 1 only!
    */
    template <class ValueT, class ReferenceT, class PointerT, class SpectrumIteratorT, class PeakIteratorT>
    class AreaIterator :
      public std::iterator<std::forward_iterator_tag, ValueT>
    {
public:
      typedef double CoordinateType;
      typedef ValueT PeakType;
      typedef SpectrumIteratorT SpectrumIteratorType;
      typedef PeakIteratorT PeakIteratorType;

      /** @name Typedefs for STL compliance
      */
      //@{
      /// The iterator's value type
      typedef ValueT value_type;
      /// The reference type as returned by operator*()
      typedef ReferenceT reference;
      /// The pointer type as returned by operator->()
      typedef PointerT pointer;
      /// The difference type
      typedef unsigned int difference_type;
      //@}

      /// Constructor for the begin iterator
      AreaIterator(SpectrumIteratorType first, SpectrumIteratorType begin, SpectrumIteratorType end, CoordinateType low_mz, CoordinateType high_mz) :
        first_(first),
        current_scan_(begin),
        end_scan_(end),
        low_mz_(low_mz),
        high_mz_(high_mz),
        is_end_(false)
      {
        nextScan_();
      }

      /// Default constructor (for the end iterator)
      AreaIterator() :
        first_(),
        current_scan_(),
        end_scan_(),
        current_peak_(),
        end_peak_(),
        low_mz_(0.0),
        high_mz_(0.0),
        is_end_(true)
      {}

      /// Destructor
      ~AreaIterator()
      {}

      /// Copy constructor
      AreaIterator(const AreaIterator & rhs) :
        first_(rhs.first_),
        current_scan_(rhs.current_scan_),
        end_scan_(rhs.end_scan_),
        current_peak_(rhs.current_peak_),
        end_peak_(rhs.end_peak_),
        low_mz_(rhs.low_mz_),
        high_mz_(rhs.high_mz_),
        is_end_(rhs.is_end_)
      {}

      /// Assignment operator
      AreaIterator & operator=(const AreaIterator & rhs)
      {
        if (&rhs == this) return *this;

        is_end_ = rhs.is_end_;
        //only copy iterators, if the assigned iterator is not the end iterator
        if (!is_end_)
        {
          first_ = rhs.first_;
          current_scan_ = rhs.current_scan_;
          end_scan_ = rhs.end_scan_;
          current_peak_ = rhs.current_peak_;
          end_peak_ = rhs.end_peak_;
          low_mz_ = rhs.low_mz_;
          high_mz_ = rhs.high_mz_;
        }

        return *this;
      }

      /// Test for equality
      bool operator==(const AreaIterator & rhs) const
      {
        //Both end iterators => equal
        if (is_end_ && rhs.is_end_) return true;

        //Normal and end iterator => not equal
        if (!is_end_ && rhs.is_end_) return false;

        if (is_end_ && !rhs.is_end_) return false;

        //Equality of pointed to peak addresses
        return &(*current_peak_) == &(*(rhs.current_peak_));
      }

      /// Test for inequality
      bool operator!=(const AreaIterator & rhs) const
      {
        return !(*this == rhs);
      }

      /// Step forward by one (prefix operator)
      AreaIterator & operator++()
      {
        //no increment if this is the end iterator
        if (is_end_) return *this;

        ++current_peak_;
        // test whether we arrived at the end of the current scan
        if (current_peak_ == end_peak_)
        {
          ++current_scan_;
          nextScan_();
        }
        return *this;
      }

      /// Step forward by one (postfix operator)
      AreaIterator operator++(int)
      {
        AreaIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      /// Dereferencing of this pointer yields the underlying peak
      reference operator*() const
      {
        return current_peak_.operator*();
      }

      /// Dereferencing of this pointer yields the underlying peak
      pointer operator->() const
      {
        return current_peak_.operator->();
      }

      /// returns the retention time of the current scan
      CoordinateType getRT() const
      {
        return current_scan_->getRT();
      }

      /// returns the PeakIndex corresponding to the current iterator position
      inline PeakIndex getPeakIndex() const
      {
        if (is_end_)
        {
          return PeakIndex();
        }
        else
        {
          return PeakIndex(current_scan_ - first_, current_peak_ - current_scan_->begin());
        }
      }

private:
      //Advances to the iterator to the next valid peak in the next valid spectrum
      void nextScan_()
      {
        while (true)
        {
          //if (current_scan_ != end_scan_) std::cout << "RT: " << current_scan_->getRT() << std::endl;
          while (current_scan_ != end_scan_ && current_scan_->getMSLevel() != 1)
          {
            ++current_scan_;
          }
          if (current_scan_ == end_scan_)
          {
            is_end_ = true;
            return;
          }
          current_peak_ = current_scan_->MZBegin(low_mz_);
          end_peak_ = current_scan_->MZEnd(high_mz_);
          if (current_peak_ != end_peak_)
          {
            return;
          }
          ++current_scan_;
        }
      }

      /// Iterator to the first scan of the map (needed to calculate the index)
      SpectrumIteratorType first_;
      /// Iterator to the current spectrum
      SpectrumIteratorType current_scan_;
      /// Past-the-end iterator of spectra
      SpectrumIteratorType end_scan_;
      /// Iterator to the current peak
      PeakIteratorType current_peak_;
      /// Past-the-end iterator of peaks in the current spectrum
      PeakIteratorType end_peak_;
      /// low m/z boundary
      CoordinateType low_mz_;
      /// high m/z boundary
      CoordinateType high_mz_;
      /// Flag that indicates that this iterator is the end iterator
      bool is_end_;

    };

  }
}

