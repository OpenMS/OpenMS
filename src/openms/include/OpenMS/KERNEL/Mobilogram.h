// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MobilityPeak1D.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/METADATA/DataArrays.h>

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <algorithm>
#include <numeric>

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
    /// Float data array vector type
    typedef OpenMS::DataArrays::FloatDataArray FloatDataArray ;
    typedef std::vector<FloatDataArray> FloatDataArrays;
    /// std::string data array vector type
    typedef OpenMS::DataArrays::StringDataArray StringDataArray ;
    typedef std::vector<StringDataArray> StringDataArrays;
    /// Integer data array vector type
    typedef OpenMS::DataArrays::IntegerDataArray IntegerDataArray ;
    typedef std::vector<IntegerDataArray> IntegerDataArrays;
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
    std::string getDriftTimeUnitAsString() const;

    /**
      @brief Sets the ion mobility drift time unit
    */
    void setDriftTimeUnit(DriftTimeUnit dt) noexcept;

    //@}

    /**
      @name Peak data array methods

      These methods are used to annotate each peak in a chromatogram with meta information.
      It is an intermediate way between storing the information in the peak's MetaInfoInterface
      and deriving a new peak type with members for this information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    ///@{
    /// Returns a const reference to the float meta data arrays
    const FloatDataArrays& getFloatDataArrays() const;

    /// Returns a mutable reference to the float meta data arrays
    FloatDataArrays& getFloatDataArrays();

    /// Sets the float meta data arrays
    void setFloatDataArrays(const FloatDataArrays& fda)
    {
      float_data_arrays_ = fda;
    }

    /// Returns a const reference to the string meta data arrays
    const StringDataArrays& getStringDataArrays() const;

    /// Returns a mutable reference to the string meta data arrays
    StringDataArrays& getStringDataArrays();

    /// Sets the string meta data arrays
    void setStringDataArrays(const StringDataArrays& sda)
    {
      string_data_arrays_ = sda;
    }

    /// Returns a const reference to the integer meta data arrays
    const IntegerDataArrays& getIntegerDataArrays() const;

    /// Returns a mutable reference to the integer meta data arrays
    IntegerDataArrays& getIntegerDataArrays();

    /// Sets the integer meta data arrays
    void setIntegerDataArrays(const IntegerDataArrays& ida)
    {
      integer_data_arrays_ = ida;
    }

    ///@}

    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
    */
    void sortByIntensity(bool reverse = false);

    /**
      @brief Lexicographically sorts the peaks by their position (mobility).

      The mobilogram is sorted with respect to position (mobility). Meta data arrays will be
      sorted accordingly.
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

    /**
      @brief Stably sorts the peaks (and all parallel data arrays) by a user-defined predicate.

      You can pass any @p lambda with signature <tt>bool(Size index_1, Size index_2)</tt> which,
      given two indices into the Mobilogram (either peaks or data arrays), returns a strict weak
      ordering. Capture the Mobilogram in the lambda and operate on it based on the indices. The
      sort is stable, so peaks with equivalent keys keep their relative order.

      @tparam Predicate Callable with signature <tt>bool(Size, Size)</tt> expressing a strict weak ordering.
      @param[in] lambda The comparison predicate.

      @note All data arrays are reordered alongside the peaks.
      @note Cached ranges are preserved (a permutation does not change them).

      @exception Exception::Precondition if a non-empty data array's size differs from the number
                 of peaks. This is checked up front, before @p lambda is ever invoked, so the
                 mobilogram is left unchanged.
    */
    template<class Predicate>
    void sort(const Predicate& lambda)
    {
      // Validate up front, before running the (possibly data-array-indexing) predicate, so a
      // mis-sized array throws instead of being read out of bounds during the sort.
      checkDataArraySizes_();
      std::vector<Size> indices(this->size());
      std::iota(indices.begin(), indices.end(), 0);
      std::stable_sort(indices.begin(), indices.end(), lambda);
      // 'indices' is a permutation of [0, size): in range and duplicate-free by construction,
      // so the unchecked path is safe and skips the redundant bounds re-scan.
      selectUnchecked(indices);
    }

    //@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific mobility

      @param[in] mb The target mobility value
      @return Returns the index of the peak.

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the mobilogram is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType mb) const;

    /**
      @brief Binary search for the peak nearest to a specific mobility given a +/- tolerance windows

      @param[in] mb The target mobility value
      @param[in] tolerance The non-negative tolerance applied to both sides of @p mb

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if mobilogram is empty

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findNearest(CoordinateType mb, CoordinateType tolerance) const;

    /**
      @brief Search for the peak nearest to a specific mobility given two +/- tolerance windows

      @param[in] mb The target mobility value
      @param[in] tolerance_left The non-negative tolerance applied left of @p mb
      @param[in] tolerance_right The non-negative tolerance applied right of @p mb

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if mobilogram is empty

      @note Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
      @note Search for the left border is done using a binary search followed by a linear scan
    */
    Int findNearest(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const;

    /**
      @brief Search for the peak with highest intensity among the peaks near to a specific mobility given two +/- tolerance windows in Th

      @param[in] mb The target mobility value
      @param[in] tolerance_left The non-negative tolerance applied left of @p mb
      @param[in] tolerance_right The non-negative tolerance applied right of @p mb

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

      Will delete (clear) all peaks contained in the mobilogram as well as any
      associated data arrays (FloatDataArrays, IntegerDataArrays,
      StringDataArrays) -- those arrays are parallel to the peaks, so retaining
      them would leave the mobilogram inconsistent.
    */
    void clear() noexcept;

    /**
      @brief Keep only the peaks at @p indices (and the matching data-array entries), in that order.

      All indices are range-checked before anything is modified, so a rejected call leaves the
      mobilogram unchanged. This is the safe default. Callers whose indices are in range
      by construction can use selectUnchecked().

      @param[in] indices Indices to keep, in the retained order. Must be < size() and unique
                 (see selectUnchecked() for the duplicate-index restriction).
      @return Reference to this Mobilogram
      @note Cached ranges are NOT recomputed. A permutation preserves them, but selecting a
            subset can leave them too wide -- call updateRanges() if you need them exact.
      @exception Exception::Precondition if an index is out of range, or if a non-empty data array's
                 size differs from the number of peaks.
    */
    Mobilogram& select(const std::vector<Size>& indices);

    /**
      @brief Like select(), but without the per-index range check.

      For callers whose indices are in range by construction (e.g. a sort() permutation).
      Data-array sizes are still validated, as that check is cheap.

      @note Duplicate indices are undefined behaviour: entries are moved, so a repeated index would
            read a moved-from value. Range and uniqueness are asserted via OPENMS_PRECONDITION in
            debug builds only.
      @note No strong exception guarantee: peaks are reordered before the data arrays, so a
            std::bad_alloc while reordering an array can leave peaks permuted and arrays not.

      @param[in] indices Indices to keep, in the retained order. Must be < size() and unique.
      @return Reference to this Mobilogram
      @exception Exception::Precondition if a non-empty data array's size differs from the number of peaks.
    */
    Mobilogram& selectUnchecked(const std::vector<Size>& indices);

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the mobilogram is unsorted.
    ConstIterator getBasePeak() const;

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the mobilogram is unsorted.
    Iterator getBasePeak();

    /// compute the total ion count (sum of all peak intensities)
    PeakType::IntensityType calculateTIC() const;

  protected:
    /**
      @brief Throws Exception::Precondition if any non-empty data array's size differs from the peak count.

      Shared by sort() and select() so that permuting the peaks can never mis-associate a parallel
      data array or index one out of bounds.

      @exception Exception::Precondition if a non-empty float, string or integer data array's size
                 differs from the number of peaks.
    */
    void checkDataArraySizes_() const;

    /// the actual peaks
    std::vector<MobilityPeak1D> data_;

    /// Retention time
    double retention_time_ = -1;

    /// Drift time unit
    DriftTimeUnit drift_time_unit_ = DriftTimeUnit::NONE;

    /// Float data arrays
    FloatDataArrays float_data_arrays_;

    /// std::string data arrays
    StringDataArrays string_data_arrays_;

    /// Integer data arrays
    IntegerDataArrays integer_data_arrays_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Mobilogram& mb);
} // namespace OpenMS
