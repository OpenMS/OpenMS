// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MobilityPeak1D.h>

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/RangeManager.h>

namespace OpenMS
{
  enum class DriftTimeUnit;
  /**
    @brief The representation of a 1D ion mobilogram.

    It contains peak data of type MobilityPeak1D.

    @note For range operations, see \ref RangeUtils "RangeUtils module"!

    @ingroup Kernel
  */
  class OPENMS_DLLAPI Mobilogram final : public RangeManagerContainer<RangeMobility, RangeIntensity>
  {
  public:
    /// Comparator for the RT of the mobilogram.
    struct OPENMS_DLLAPI RTLess {
      bool operator()(const Mobilogram& a, const Mobilogram& b) const;
    };


    ///@name Base type definitions
    //@{
    /// Peak type
    using PeakType = MobilityPeak1D;
    /// Coordinate (mobility) type
    using CoordinateType = PeakType::CoordinateType;
    /// Mobilogram base type
    using ContainerType = std::vector<PeakType>;
    /// RangeManager
    using RangeManagerContainerType = RangeManagerContainer<RangeMobility, RangeIntensity>;
    using RangeManagerType = RangeManager<RangeMobility, RangeIntensity>;
    //@}

    ///@name Peak container iterator type definitions
    //@{
    /// Mutable iterator
    using Iterator = ContainerType::iterator;
    using iterator = Iterator;
    /// Non-mutable iterator
    using ConstIterator = ContainerType::const_iterator;
    using const_iterator = ConstIterator;
    /// Mutable reverse iterator
    using ReverseIterator = ContainerType::reverse_iterator;
    using reverse_iterator = ReverseIterator;
    /// Non-mutable reverse iterator
    using ConstReverseIterator = ContainerType::const_reverse_iterator;
    using const_reverse_iterator = ConstReverseIterator;
    //@}
    /*using typename ContainerType::const_reference;
    using typename ContainerType::difference_type;
    using typename ContainerType::pointer;
    using typename ContainerType::reference;
    using typename ContainerType::size_type;
    using typename ContainerType::value_type;*/

    // rule of 6

    /// Constructor
    Mobilogram() = default;

    /// Copy constructor
    Mobilogram(const Mobilogram& source) = default;

    /// Move constructor
    Mobilogram(Mobilogram&&) noexcept = default;

    /// Assignment operator
    Mobilogram& operator=(const Mobilogram& source) = default;

    /// Move assignment operator
    Mobilogram& operator=(Mobilogram&&) noexcept = default;

    /// Destructor
    ~Mobilogram() = default;


    /// Equality operator
    bool operator==(const Mobilogram& rhs) const;

    /// Equality operator
    bool operator!=(const Mobilogram& rhs) const
    {
      return !(operator==(rhs));
    }

    ///@name Export methods for std::vector<MobilityPeak1D>
    //@{
    MobilityPeak1D& operator[](Size i) noexcept
    {
      return data_[i];  
    }
    const MobilityPeak1D& operator[](Size i) const noexcept
    {
      return data_[i];
    }


    MobilityPeak1D& front() noexcept
    {
      return data_.front();
    }
    const MobilityPeak1D& front() const noexcept
    {
      return data_.front();
    }

    MobilityPeak1D& back() noexcept
    {
      return data_.back();
    }
    const MobilityPeak1D& back() const noexcept
    {
      return data_.back();
    }

    Iterator begin() noexcept
    {
      return data_.begin();
    }
    ConstIterator begin() const noexcept
    {
      return data_.begin();
    }
    ConstIterator cbegin() const noexcept
    {
      return data_.cbegin();
    }

    Iterator end() noexcept
    {
      return data_.end();
    }
    ConstIterator end() const noexcept
    {
      return data_.end();
    }
    ConstIterator cend() const noexcept
    {
      return data_.cend();
    }

    ReverseIterator rbegin() noexcept
    {
      return data_.rbegin();
    }
    ConstReverseIterator crbegin() const
    {
      return data_.crbegin();
    }
    ReverseIterator rend() noexcept
    {
      return data_.rend();
    }
    ConstReverseIterator crend() const
    {
      return data_.crend();
    }

    bool empty() const noexcept
    {
      return data_.empty();
    }
    ConstIterator erase(ConstIterator where) noexcept
    {
      return data_.erase(where);
    }

    void push_back(MobilityPeak1D mb)
    {
      data_.push_back(mb);
    }
    MobilityPeak1D& emplace_back(MobilityPeak1D mb)
    {
      return data_.emplace_back(mb);
    }
    template<class... Args>
    void emplace_back(Args&&... args)
    {
      data_.emplace_back(args...);
    }

    void pop_back()
    {
      data_.pop_back();
    }

    Iterator insert(ConstIterator where, ConstIterator first, ConstIterator last)
    {
      return data_.insert(where, first, last);
    }

    void resize(size_t new_size)
    {
      return data_.resize(new_size);
    }
    void reserve(size_t new_size)
    {
      return data_.reserve(new_size);
    }

    size_t size() const noexcept
    {
      return data_.size();
    }

    void swap(Mobilogram& mb) noexcept
    {
      data_.swap(mb.data_);
      std::swap(retention_time_, mb.retention_time_);
      std::swap(drift_time_unit_, mb.drift_time_unit_);
    }
    //@}

    // Docu in base class (RangeManager)
    void updateRanges() override;

    ///@name Accessors for meta information
    ///@{
    /// Returns the retention time (in seconds)
    double getRT() const noexcept
    {
      return retention_time_;
    }

    /// Sets the retention time (in seconds)
    void setRT(double rt) noexcept
    {
      retention_time_ = rt; 
    }

    /**
      @brief Returns the ion mobility drift time unit
    */
    DriftTimeUnit getDriftTimeUnit() const noexcept
    {
      return drift_time_unit_;
    }

    /// returns the ion mobility drift time unit as string
    String getDriftTimeUnitAsString() const;

    /**
      @brief Sets the ion mobility drift time unit
    */
    void setDriftTimeUnit(DriftTimeUnit dt) noexcept;

    //@}


    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity.
    */
    void sortByIntensity(bool reverse = false);

    /**
      @brief Lexicographically sorts the peaks by their position (mobility).

      The mobilogram is sorted with respect to position (mobility).
    */
    void sortByPosition();

    /// Checks if all peaks are sorted with respect to ascending mobility
    bool isSorted() const;

    /// Checks if container is sorted by a certain user-defined property.
    /// You can pass any lambda function with <tt>[](Size index_1, Size index_2) --> bool</tt>
    /// which given two indices into Mobilogram (either for peaks or data arrays) returns a weak-ordering.
    /// (you need to capture the Mobilogram in the lambda and operate on it, based on the indices)
    template<class Predicate>
    bool isSorted(const Predicate& lambda) const
    {
      auto value_2_index_wrapper = [this, &lambda](const PeakType& value1, const PeakType& value2) {
        // translate values into indices (this relies on no copies being made!)
        const Size index1 = (&value1) - (&this->front());
        const Size index2 = (&value2) - (&this->front());
        // just make sure the pointers above are actually pointing to a Peak inside our container
        assert(index1 < this->size());
        assert(index2 < this->size());
        return lambda(index1, index2);
      };
      return std::is_sorted(this->begin(), this->end(), value_2_index_wrapper);
    }

    //@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific mobility

      @param mb The target mobility value
      @return Returns the index of the peak.

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the mobilogram is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType mb) const;

    /**
      @brief Binary search for the peak nearest to a specific mobility given a +/- tolerance windows

      @param mb The target mobility value
      @param tolerance The non-negative tolerance applied to both sides of @p mb

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if mobilogram is empty

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findNearest(CoordinateType mb, CoordinateType tolerance) const;

    /**
      @brief Search for the peak nearest to a specific mobility given two +/- tolerance windows

      @param mb The target mobility value
      @param tolerance_left The non-negative tolerance applied left of @p mb
      @param tolerance_right The non-negative tolerance applied right of @p mb

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if mobilogram is empty

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
      @note Search for the left border is done using a binary search followed by a linear scan
    */
    Int findNearest(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const;

    /**
      @brief Search for the peak with highest intensity among the peaks near to a specific mobility given two +/- tolerance windows in Th

      @param mb The target mobility value
      @param tolerance_left The non-negative tolerance applied left of @p mb
      @param tolerance_right The non-negative tolerance applied right of @p mb

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if mobilogram is empty

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findHighestInWindow(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    Iterator MBBegin(CoordinateType mb);

    /**
      @brief Binary search for peak range begin

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    Iterator MBBegin(Iterator begin, CoordinateType mb, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    Iterator MBEnd(CoordinateType mb);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    Iterator MBEnd(Iterator begin, CoordinateType mb, Iterator end);

    /**
      @brief Binary search for peak range begin

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    ConstIterator MBBegin(CoordinateType mb) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    ConstIterator MBBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    ConstIterator MBEnd(CoordinateType mb) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    ConstIterator MBEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const;

    /**
      @brief Binary search for peak range begin

      Alias for MBBegin()

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    Iterator PosBegin(CoordinateType mb);

    /**
      @brief Binary search for peak range begin

      Alias for MBBegin()

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    Iterator PosBegin(Iterator begin, CoordinateType mb, Iterator end);

    /**
      @brief Binary search for peak range begin

      Alias for MBBegin()

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    ConstIterator PosBegin(CoordinateType mb) const;

    /**
      @brief Binary search for peak range begin

      Alias for MBBegin()

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
    */
    ConstIterator PosBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MBEnd()

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    Iterator PosEnd(CoordinateType mb);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MBEnd()

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    Iterator PosEnd(Iterator begin, CoordinateType mb, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MBEnd()

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    ConstIterator PosEnd(CoordinateType mb) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MBEnd()

      @note Make sure the mobilogram is sorted with respect to mobility. Otherwise the result is undefined.
    */
    ConstIterator PosEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const;

    //@}


    /**
      @brief Clears all data and ranges

      Will delete (clear) all peaks contained in the mobilogram 
    */
    void clear() noexcept;

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the mobilogram is unsorted.
    ConstIterator getBasePeak() const;

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the mobilogram is unsorted.
    Iterator getBasePeak();

    /// compute the total ion count (sum of all peak intensities)
    PeakType::IntensityType calculateTIC() const;

  protected:
    /// the actual peaks
    std::vector<MobilityPeak1D> data_;

    /// Retention time
    double retention_time_ = -1;

    /// Drift time unit
    DriftTimeUnit drift_time_unit_ = DriftTimeUnit::NONE;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Mobilogram& mb);
} // namespace OpenMS
