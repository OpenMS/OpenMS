// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS includes
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/PeakIndex.h>
#include <OpenMS/KERNEL/RangeManager.h>

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

        Ion mobility can also be filtered for: the low/high range for IM are used to skip over spectra in the given RT range whose drift time
        is not within the given range. I.e. the RT range could contain multiple IM frames.

        @note This iterator iterates over spectra with same MS level as the MS level of the begin() spectrum in Param! You can explicitly set another MS level as well.
    */
    template<class ValueT, class ReferenceT, class PointerT, class SpectrumIteratorT, class PeakIteratorT>
    class AreaIterator
    {
    public:
      typedef double CoordinateType;
      typedef ValueT PeakType;
      typedef SpectrumIteratorT SpectrumIteratorType;
      typedef PeakIteratorT PeakIteratorType;
      using SpectrumT = typename std::iterator_traits<SpectrumIteratorType>::value_type;

      /// Parameters for the AreaIterator
      /// Required values must be set in the C'tor. Optional values can be set via member functions (which allow chaining).
      class Param
      {
        friend AreaIterator; // allow access to private members (avoids writing get-accessors)
      public:
        /**
         * \brief C'tor with mandatory parameters
         * \param first The very first spectrum of the experiment
         * \param begin The first spectrum with a valid RT/IM time
         * \param end The last spectrum with a valid RT/IM time
         * \param ms_level Only peaks from spectra with this ms_level are used
         */
        Param(SpectrumIteratorType first, SpectrumIteratorType begin, SpectrumIteratorType end, uint8_t ms_level) : first_(first), current_scan_(begin), end_scan_(end), ms_level_(ms_level)
        {
        }

        /// return the end-iterator
        static Param end()
        {
          static Param p;
          p.is_end_ = true;
          return p;
        }

        /// Assignment operator
        Param& operator=(const Param& rhs) = default;

        /** @name Named parameter idiom for chaining
         */
        //@{
        /// low m/z boundary
        Param& lowMZ(CoordinateType low_mz)
        {
          low_mz_ = low_mz;
          return *this;
        }
        /// high m/z boundary
        Param& highMZ(CoordinateType high_mz)
        {
          high_mz_ = high_mz;
          return *this;
        }
        /// low ion mobility boundary
        Param& lowIM(CoordinateType low_im)
        {
          low_im_ = low_im;
          return *this;
        }
        /// high ion mobility boundary
        Param& highIM(CoordinateType high_im)
        {
          high_im_ = high_im;
          return *this;
        }
        /// Only scans of this MS level are iterated over
        Param& msLevel(int8_t ms_level)
        {
          ms_level_ = ms_level;
          return *this;
        }
        //@}

      protected:
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
        /* optional parameters */
        /// low m/z boundary
        CoordinateType low_mz_ = std::numeric_limits<CoordinateType>::lowest();
        /// high m/z boundary
        CoordinateType high_mz_ = std::numeric_limits<CoordinateType>::max();
        /// low mobility boundary
        CoordinateType low_im_ = std::numeric_limits<CoordinateType>::lowest();
        /// high mobility boundary
        CoordinateType high_im_ = std::numeric_limits<CoordinateType>::max();
        /// Only scans of this MS level are iterated over
        int8_t ms_level_ {};
        /// Flag that indicates that this iterator is the end iterator
        bool is_end_ = false;

      private:
        /// only used internally for end()
        Param() = default;
      };


      /** @name Typedefs for STL compliance, these replace std::iterator
       */
      //@{
      /// The iterator's category type
      typedef std::forward_iterator_tag iterator_category;
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
      explicit AreaIterator(const Param& p) : p_(p)
      {
        nextScan_();
      }

      /// Default constructor (for the end iterator)
      AreaIterator() : p_(Param::end())
      {
      }

      /// Destructor
      ~AreaIterator() = default;

      /// Copy constructor
      AreaIterator(const AreaIterator& rhs) = default;

      /// Assignment operator
      AreaIterator& operator=(const AreaIterator& rhs)
      {
        p_.is_end_ = rhs.p_.is_end_;
        // only copy iterators, if the assigned iterator is not the end iterator
        if (!p_.is_end_)
        {
          p_ = rhs.p_;
        }

        return *this;
      }

      /// Test for equality
      bool operator==(const AreaIterator& rhs) const
      {
        // Both end iterators => equal
        if (p_.is_end_ && rhs.p_.is_end_)
          return true;

        // Normal and end iterator => not equal
        if (p_.is_end_ ^ rhs.p_.is_end_)
          return false;

        // Equality of pointed to peak addresses
        return &(*(p_.current_peak_)) == &(*(rhs.p_.current_peak_));
      }

      /// Test for inequality
      bool operator!=(const AreaIterator& rhs) const
      {
        return !(*this == rhs);
      }

      /// Step forward by one (prefix operator)
      AreaIterator& operator++()
      {
        // no increment if this is the end iterator
        if (p_.is_end_)
          return *this;

        ++p_.current_peak_;
        // test whether we arrived at the end of the current scan
        if (p_.current_peak_ == p_.end_peak_)
        {
          ++p_.current_scan_;
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
        return p_.current_peak_.operator*();
      }

      /// Dereferencing of this pointer yields the underlying peak
      pointer operator->() const
      {
        return p_.current_peak_.operator->();
      }

      /// returns the retention time of the current scan
      CoordinateType getRT() const
      {
        return p_.current_scan_->getRT();
      }

      /// returns the ion mobility time of the current scan
      CoordinateType getDriftTime() const
      {
        return p_.current_scan_->getDriftTime();
      }

      /// returns the current scan into which the iterator points
      const SpectrumT& getSpectrum() const
      {
        return *p_.current_scan_;
      }

      /// returns the PeakIndex corresponding to the current iterator position
      inline PeakIndex getPeakIndex() const
      {
        if (p_.is_end_)
        {
          return {};
        }
        else
        {
          return PeakIndex(p_.current_scan_ - p_.first_, p_.current_peak_ - p_.current_scan_->begin());
        }
      }

    private:
      /// advances the iterator to the next valid peak in the next valid spectrum
      void nextScan_()
      {
        using MSLevelType = decltype(p_.current_scan_->getMSLevel());
        RangeMobility mb {p_.low_im_, p_.high_im_};
        while (true)
        {
          // skip over invalid MS levels and Mobility
          while (p_.current_scan_ != p_.end_scan_ && (p_.current_scan_->getMSLevel() != (MSLevelType)p_.ms_level_ || !mb.containsMobility(p_.current_scan_->getDriftTime())))
          {
            ++p_.current_scan_;
          }
          if (p_.current_scan_ == p_.end_scan_)
          {
            p_.is_end_ = true;
            return;
          }
          p_.current_peak_ = p_.current_scan_->MZBegin(p_.low_mz_);
          p_.end_peak_ = p_.current_scan_->MZEnd(p_.high_mz_);
          if (p_.current_peak_ != p_.end_peak_)
          {
            return;
          }
          ++p_.current_scan_;
        }
      }

      /// holds spectra iterators and area limits
      Param p_;
    };

  } // namespace Internal
} // namespace OpenMS
