// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm> // for min/max
#include <cassert>
#include <iosfwd>  // for std::ostream


namespace OpenMS
{
  /// Dimensions of data acquisition for MS data
  enum class MSDim
  {
    RT,
    MZ,
    INT,
    IM
  };

  /// Base class for a simple range with minimum and maximum
  struct OPENMS_DLLAPI RangeBase
  {
  public:
    /// Ctor: initialize with empty range
    RangeBase() = default;
    
    /// set min and max
    /// @throws Exception::InvalidRange if min > max
    RangeBase(const double min, const double max) :
        min_(min), max_(max)
    {
      if (min_ > max_)
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid initialization of range");
    }

    /// make the range empty, i.e. isEmpty() will be true
    void clear()
    {
      *this = RangeBase(); // overwrite with fresh instance
    }

    /// is the range empty (i.e. default constructed or cleared using clear())?
    bool isEmpty() const
    { // invariant: only possible when default constructed or clear()'ed
      return min_ > max_; 
    }

    /// is @p value within [min, max]?
    bool contains(const double value) const
    {
      return min_ <= value && value <= max_;
    }

    /// is the range @p inner_range within [min, max]?
    bool contains(const RangeBase& inner_range) const
    {
      return contains(inner_range.min_) && contains(inner_range.max_);
    }

    /** @name Accessors for min and max
      
        We use accessors, to keep range consistent (i.e. ensure that min <= max)
    */
    ///@{  
    
    /// sets the minimum (and the maximum, if uninitialized)
    void setMin(const double min)
    {
      min_ = min;
      if (max_ < min)
        max_ = min;
    }

    /// sets the maximum (and the minimum, if uninitialized)
    void setMax(const double max)
    {
      max_ = max;
      if (min_ > max)
        min_ = max;
    }

    /// only useful if isEmpty() returns false
    double getMin() const
    {
      return min_;
    }

    /// only useful if isEmpty() returns false
    double getMax() const
    {
      return max_;
    }
    ///@}

    /// ensure the range includes the range of @p other
    void extend(const RangeBase& other)
    {
      min_ = std::min(min_, other.min_);
      max_ = std::max(max_, other.max_);
    }

    /// extend the range such that it includes the given @p value
    void extend(const double value)
    {
      min_ = std::min(min_, value);
      max_ = std::max(max_, value);
    }

    /**
       @brief Scale the range of the dimension by a @p factor. A factor > 1 increases the range; factor < 1 decreases it.

       Let d = max - min; then min = min - d*(factor-1)/2,
       i.e. scale(1.5) extends the range by 25% on each side.
 
       Scaling an empty range will not have any effect.

       @param factor The multiplier to increase the range by
    */
    void scaleBy(const double factor)
    {
      if (isEmpty()) return;
      const double dist = max_ - min_;
      const double extension = dist * (factor - 1) / 2;
      min_ -= extension;
      max_ += extension;
    }

    void assign(const RangeBase& rhs)
    {
      *this = rhs;
    }

    bool operator==(const RangeBase& rhs) const
    {
      return min_ == rhs.min_ && max_ == rhs.max_;
    }

  protected:
    // make members non-accessible to maintain invariant: min <= max  (unless uninitialized)
    double min_ = std::numeric_limits<double>::max();
    double max_ = -std::numeric_limits<double>::max();
  };
  
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeBase& b);

  struct OPENMS_DLLAPI RangeRT : public RangeBase {

    const static MSDim DIM = MSDim::RT;

    RangeRT() = default;
    RangeRT(const double min, const double max) :
        RangeBase(min, max)
    {
    }

    /** @name Accessors for min and max
      
        We use accessors, to keep range consistent (i.e. ensure that min <= max)
    */
    ///@{

    /// sets the minimum (and the maximum, if uninitialized)
    void setMinRT(const double min)
    {
      setMin(min);
    }

    /// sets the maximum (and the minimum, if uninitialized)
    void setMaxRT(const double max)
    {
      setMax(max);
    }

    /// only useful if isEmpty() returns false
    double getMinRT() const
    {
      return min_;
    }

    /// only useful if isEmpty() returns false
    double getMaxRT() const
    {
      return max_;
    }
    ///@}

    /// extend the range such that it includes the given @p value
    void extendRT(const double value)
    {
      extend(value);
    }

    /// is @p value within [min, max]?
    bool containsRT(const double value) const
    {
      return RangeBase::contains(value);
    }

    /// is the range @p inner_range within [min, max] of this range?
    bool containsRT(const RangeBase& inner_range) const
    {
      return RangeBase::contains(inner_range);
    }
  };
  
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeRT& range);

  struct OPENMS_DLLAPI RangeMZ : public RangeBase
  {

    const static MSDim DIM = MSDim::MZ;

    RangeMZ() = default;
    RangeMZ(const double min, const double max) :
        RangeBase(min, max)
    {
    }

    /** @name Accessors for min and max
      
        We use accessors, to keep range consistent (i.e. ensure that min <= max)
    */
    ///@{

    /// sets the minimum (and the maximum, if uninitialized)
    void setMinMZ(const double min)
    {
      setMin(min);
    }

    /// sets the maximum (and the minimum, if uninitialized)
    void setMaxMZ(const double max)
    {
      setMax(max);
    }

    /// only useful if isEmpty() returns false
    double getMinMZ() const
    {
      return min_;
    }

    /// only useful if isEmpty() returns false
    double getMaxMZ() const
    {
      return max_;
    }
    ///@}

    /// extend the range such that it includes the given @p value
    void extendMZ(const double value)
    {
      extend(value);
    }

    /// is @p value within [min, max]?
    bool containsMZ(const double value) const
    {
      return RangeBase::contains(value);
    }

    /// is the range @p inner_range within [min, max] of this range?
    bool containsMZ(const RangeBase& inner_range) const
    {
      return RangeBase::contains(inner_range);
    }
  };
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeMZ& range);

  struct OPENMS_DLLAPI RangeIntensity : public RangeBase {

    const static MSDim DIM = MSDim::INT;

    RangeIntensity() = default;
    RangeIntensity(const double min, const double max) :
        RangeBase(min, max)
    {
    }

    /** @name Accessors for min and max
      
        We use accessors, to keep range consistent (i.e. ensure that min <= max)
    */
    ///@{

    /// sets the minimum (and the maximum, if uninitialized)
    void setMinIntensity(const double min)
    {
      setMin(min);
    }

    /// sets the maximum (and the minimum, if uninitialized)
    void setMaxIntensity(const double max)
    {
      setMax(max);
    }

    /// only useful if isEmpty() returns false
    double getMinIntensity() const
    {
      return min_;
    }

    /// only useful if isEmpty() returns false
    double getMaxIntensity() const
    {
      return max_;
    }
    ///@}

    /// extend the range such that it includes the given @p value
    void extendIntensity(const double value)
    {
      extend(value);
    }

    /// is @p value within [min, max]?
    bool containsIntensity(const double value) const
    {
      return RangeBase::contains(value);
    }

    /// is the range @p inner_range within [min, max] of this range?
    bool containsIntensity(const RangeBase& inner_range) const
    {
      return RangeBase::contains(inner_range);
    }
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeIntensity& range);

  struct OPENMS_DLLAPI RangeMobility : public RangeBase
  {
    const static MSDim DIM = MSDim::IM;

    RangeMobility() = default;
    RangeMobility(const double min, const double max) :
        RangeBase(min, max)
    {
    }

    /** @name Accessors for min and max
      
        We use accessors, to keep range consistent (i.e. ensure that min <= max)
    */
    ///@{

    /// sets the minimum (and the maximum, if uninitialized)
    void setMinMobility(const double min)
    {
      setMin(min);
    }

    /// sets the maximum (and the minimum, if uninitialized)
    void setMaxMobility(const double max)
    {
      setMax(max);
    }

    /// only useful if isEmpty() returns false
    double getMinMobility() const
    {
      return min_;
    }

    /// only useful if isEmpty() returns false
    double getMaxMobility() const
    {
      return max_;
    }
    ///@}

    /// extend the range such that it includes the given @p value
    void extendMobility(const double value)
    {
      extend(value);
    }

    /// is @p value within [min, max]?
    bool containsMobility(const double value) const
    {
      return RangeBase::contains(value);
    }

    /// is the range @p inner_range within [min, max] of this range?
    bool containsMobility(const RangeBase& inner_range) const
    {
      return RangeBase::contains(inner_range);
    }
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeMobility& range);
  
  /// Enum listing state of dimensions (RangeBases)
  enum class HasRangeType
  {
    ALL, ///< all dimensions are filled
    SOME,///< some dimensions are empty, some are filled
    NONE ///< all dimensions are empty (=cleared)
  };

  /**
    @brief Handles the management of a multidimensional range, e.g. RangeMZ and RangeIntensity for spectra.

    Instanciate it with the dimensions which are supported/required, e.g.
    <tt>RangeManager<RangeRT, RangeMZ> range_spec</tt> for a spectrum and use the strongly typed features, such as
    range_spec.getMaxRT()/setMaxRT(500.0) or range_spec.extend(RangeMZ{100, 1500});

    Use RangeManagerContainer as a base class for all peak and feature containers like MSSpectrum, MSExperiment and FeatureMap.

    The implementation uses non-virtual multiple inheritance using variadic templates. Each dimension, e.g. RangeRT, is inherited from, thus
    all members of the base class become accessible in the RangeManager, e.g. ::getMaxRT().
    Operations (e.g. assignment, or extension of ranges) across RangeManagers with a different, yet overlapping set of base classes
    is enabled using fold expressions and constexpr evaluations, which are resolved at compile time (see for_each_base_ member function).

  */
  template<typename... RangeBases>
  class RangeManager : public RangeBases...
  {
  public:
    // rule of 0 -- no need for a virtual d'tor or anything fancy
    // ...

    bool operator==(const RangeManager& rhs) const
    {
      bool equal = true;
      for_each_base_([&](auto* base) {
        using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
        equal &= ((T_BASE&) rhs == (T_BASE&) *this);
      });
      return equal;
    }


    /// copy all overlapping dimensions from @p rhs to this instance.
    /// Dimensions which are not contained in @p rhs are left untouched.
    /// @param rhs Range to copy from
    /// @return true if one or more dimensions overlapped
    template <typename... RangeBasesOther>
    bool assignUnsafe(const RangeManager<RangeBasesOther...>& rhs)
    {
      bool found = false;
      for_each_base_([&](auto* base) {
        using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
        if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
        {
          base->assign((T_BASE&) rhs);
          found = true;
        }
      });

      return found;
    }

    /// copy all overlapping dimensions from @p rhs to this instance.
    /// Dimensions which are not contained in @p rhs are left untouched.
    /// @param rhs Range to copy from
    /// @throw Exception::InvalidRange if no dimensions overlapped
    template<typename... RangeBasesOther>
    void assign(const RangeManager<RangeBasesOther...>& rhs)
    {
      if (!assignUnsafe(rhs))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION , "No assignment took place (no dimensions in common!);");
      }
    }                 

    /// extend all dimensions which overlap with @p rhs to contain the range of @p rhs
    /// Dimensions which are not contained in @p rhs are left untouched.
    /// @param rhs Range to extend from
    /// @return false if no dimensions overlapped
    template<typename... RangeBasesOther>
    bool extendUnsafe(const RangeManager<RangeBasesOther...>& rhs)
    {
      bool found = false;
      for_each_base_([&](auto* base) {
        using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
        if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
        {
          base->extend((T_BASE&) rhs);
          found = true;
        }
      });
      return found;
    }

    /// extend all dimensions which overlap with @p rhs to contain the range of @p rhs
    /// Dimensions which are not contained in @p rhs are left untouched.
    /// @param rhs Range to extend from
    /// @throw Exception::InvalidRange if no dimensions overlapped
    template<typename... RangeBasesOther>
    void extend(const RangeManager<RangeBasesOther...>& rhs)
    {
      if (!extendUnsafe(rhs))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);");
      }
    }

    /// calls RangeBase::scale() for each dimension
    void scaleBy(const double factor)
    {
      for_each_base_([&](auto* base) {
        base->scaleBy(factor);
      });
    }

    /// obtain a range dimension at runtime using @p dim
    RangeBase& getRangeForDim(MSDim dim)
    {
      RangeBase* r_base = nullptr;

      static_for_each_base_([&](auto* base) {
        using Base = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
        if (base->DIM == dim)
          r_base = (Base*) this;
      });

      assert((r_base != nullptr) && "No base class has this MSDim!");
      return *r_base;
    }
    
    /// is any/some/all dimension in this range populated?
    HasRangeType hasRange() const
    {
      constexpr size_t total{sizeof...(RangeBases)};// total number of bases
      size_t count{0};
      for_each_base_([&](auto* base) {
        count += !base->isEmpty();
      });
      switch (count)
      {
        case 0:
          return HasRangeType::NONE;
        case total:
          return HasRangeType::ALL;
        default:
          return HasRangeType::SOME;
      }
    }

    /// Are all dimensions of @p rhs (which overlap with this Range) contained in this range?
    /// An empty dimension is considered contained in the other dimension (even if that one is empty as well).
    /// If only all overlapping dimensions are empty, true is returned.
    /// @throws Exception::InvalidRange if no dimensions overlap
    template<typename... RangeBasesOther>
    bool containsAll(const RangeManager<RangeBasesOther...>& rhs) const
    {
      bool contained = true; // assume rhs is contained, until proven otherwise
      bool has_overlap = false;
      for_each_base_([&](auto* base) {
        using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
        if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
        {
          has_overlap = true; // at least one dimension overlaps
          if (((T_BASE&)rhs).isEmpty()) return;
          if (base->contains((T_BASE&) rhs)) return;
          contained = false;
        }
      });
      if (!has_overlap) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

      return contained;
    }

    /// Resets all ranges
    void clearRanges()
    {
      for_each_base_([&](auto* base) {
        base->clear();
      });
    }

    /// print each dimension (base classes) to a stream
    void printRange(std::ostream& out) const
    {
      for_each_base_([&](auto* base) {
        out << *base;
      });
    }

  protected:
    /// use fold expression to iterate over all RangeBases of RangeManager and apply a lambda (Visitor) for each one
    template<typename Visitor>
    void for_each_base_(Visitor&& visitor)
    {
      (void(visitor(static_cast<RangeBases*>(this))), ...);
    }
    /// .. and a const version
    template<typename Visitor>
    void for_each_base_(Visitor&& visitor) const
    {
      (void(visitor(static_cast<const RangeBases*>(this))), ...);
    }

    /// use fold expression to iterate over all RangeBases of RangeManager and apply a lambda (Visitor) for each one (for static members)
    template<typename Visitor>
    static void static_for_each_base_(Visitor&& visitor)
    {
      (void(visitor(static_cast<const RangeBases*>(nullptr))), ...);
    }
  };

  template<typename... Range>
  std::ostream& operator<<(std::ostream& out, const RangeManager<Range...>& me)
  {
    me.printRange(out);
    return out;
  }

  /// use this class as a base class for containers, e.g. MSSpectrum. It forces them to implement 'updateRanges()' as a common interface 
  /// and provides a 'getRange()' member which saves casting to a range type manually
  template <typename ...RangeBases>
  class RangeManagerContainer
    : public RangeManager<RangeBases...>
  {
    public:
    /// implement this function to reflect the underlying data of the derived class (e.g. an MSSpectrum)
    /// Usually, call clearRanges() internally and then populate the dimensions.
    virtual void updateRanges() = 0;

    /// get range of current data (call updateRanges() before to ensure the range is accurate)
    const RangeManager<RangeBases...>& getRange() const
    {
      return (RangeManager<RangeBases...>&) *this;
    }

    /// get mutable range, provided for efficiency reasons (avoid updateRanges(), if only minor changes were made)
    RangeManager<RangeBases...>& getRange()
    {
      return (RangeManager<RangeBases...>&)*this;
    }
    
  };

}  // namespace OpenMS
