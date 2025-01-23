// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>
#include <algorithm> // for min/max
#include <cassert>
#include <cmath>  // for nan()
#include <iosfwd> // for std::ostream

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

struct RangeRT;
struct RangeMZ;
struct RangeIntensity;
struct RangeMobility;

/// Base class for a simple range with minimum and maximum
struct OPENMS_DLLAPI RangeBase
{
public:
  /// C'tor: initialize with empty range
  RangeBase() = default;

  /// Cutom C'tor which sets the range to a singular point
  RangeBase(const double single): min_(single), max_(single)
  {
  }

  /// Custom C'tor to set min and max
  /// @throws Exception::InvalidRange if min > max
  RangeBase(const double min, const double max): min_(min), max_(max)
  {
    if (min_ > max_) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid initialization of range");
  }
  /// Copy C'tor
  RangeBase(const RangeBase& rhs) = default;

  /// Move C'tor (seems useless, but is required for completeness in derived classes' move c'tor)
  RangeBase(RangeBase&& rhs) noexcept = default;

  /// Assignment operator
  RangeBase& operator=(const RangeBase& rhs) = default;

  /// Move assignment (seems useless, but is required for completeness in derived classes' move c'tor)
  RangeBase& operator=(RangeBase&& rhs) noexcept = default;

  /// D'tor
  ~RangeBase() noexcept = default;

  /// conversion operator to allow accepting a RangeBase (instead of RangeRT) for the implicitly defined special members, e.g. assignment operator
  /// (RangeRT& operator=(const RangeRT&))
  operator RangeRT() const;
  /// conversion operator to allow accepting a RangeBase (instead of RangeMZ) for the implicitly defined special members, e.g. assignment operator
  /// (RangeMZ& operator=(const RangeMZ&))
  operator RangeMZ() const;
  /// conversion operator to allow accepting a RangeBase (instead of RangeIntensity) for the implicitly defined special members, e.g. assignment
  /// operator (RangeIntensity& operator=(const RangeIntensity&))
  operator RangeIntensity() const;
  /// conversion operator to allow accepting a RangeBase (instead of RangeMobility) for the implicitly defined special members, e.g. assignment
  /// operator (RangeMobility& operator=(const RangeMobility&))
  operator RangeMobility() const;

  /// make the range empty, i.e. isEmpty() will be true
  void clear()
  {
    *this = RangeBase(); // overwrite with fresh instance
  }

  /// is the range empty (i.e. min > max)?
  bool isEmpty() const
  {
    return min_ > max_;
  }

  /// is @p value within [min, max]?
  bool contains(const double value) const
  {
    return uint8_t(min_ <= value) & uint8_t(value <= max_); // using && leads to branches on all compilers in Debug and in Release on MVSC
  }

  /// is the range @p inner_range within [min, max]?
  bool contains(const RangeBase& inner_range) const
  {
    return uint8_t(contains(inner_range.min_))
           & uint8_t(contains(inner_range.max_)); // using && leads to branches on all compilers in Debug and in Release on MVSC
  }

  /** @name Accessors for min and max

      We use accessors, to keep range consistent (i.e. ensure that min <= max)
  */
  ///@{

  /// sets the minimum (and the maximum, if uninitialized)
  void setMin(const double min)
  {
    min_ = min;
    if (max_ < min) max_ = min;
  }

  /// sets the maximum (and the minimum, if uninitialized)
  void setMax(const double max)
  {
    max_ = max;
    if (min_ > max) min_ = max;
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

  /// Extend the range by @p by units to left and right
  /// Using negative values will shrink the range. It may become empty.
  /// Calling this on an empty range will not have any effect.
  void extendLeftRight(const double by)
  {
    if (isEmpty()) return;
    min_ -= by;
    max_ += by;
  }

  /**
   * \brief If the current range is a single point, e.g. min==max, then extend the range by @p min_span / 2 on either side.
   *
   * Calling span() afterwards, returns @p min_span.
   */
  void minSpanIfSingular(const double min_span)
  {
    if (min_ == max_) extendLeftRight(min_span / 2);
  }

  /// Ensure the range of this does not exceed the range of @p other.
  /// If @p other already contains() this range, nothing changes.
  /// If this range is entirely outside the range of @p other,
  /// the resulting range is empty.
  /// Empty ranges are not modified.
  /// @throw Exception::InvalidRange if @p other is empty
  void clampTo(const RangeBase& other)
  {
    if (isEmpty()) return;
    if (other.isEmpty()) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

    min_ = std::max(min_, other.min_);
    max_ = std::min(max_, other.max_);
  }

  /// Move range of *this to min/max of @p sandbox, without changing the span, if possible.
  /// This does tighten the range unless @p sandbox's ranges are smaller than *this.
  /// Empty ranges are not modified.
  /// @param sandbox Range to translate/move the current range into
  /// @throw Exception::InvalidRange if @p sandbox is empty
  void pushInto(const RangeBase& sandbox)
  {
    if (isEmpty()) return;
    if (sandbox.isEmpty()) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

    if (! sandbox.contains(*this))
    {
      if (getSpan() > sandbox.getSpan())
      { // make interval fit into sandbox (= ensure full containment)
        max_ = min_ + sandbox.getSpan();
      }
      if (min_ < sandbox.min_)
      { // need to shift right (positive shift)
        shift(sandbox.min_ - min_);
      }
      else if (max_ > sandbox.max_)
      { // need to shift left (negative shift)
        shift(sandbox.max_ - max_);
      }
    }
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

  /// Move the range by @p distance (negative values shift left)
  /// Shifting an empty range will not have any effect.
  void shift(const double distance)
  {
    if (isEmpty()) return;
    min_ += distance;
    max_ += distance;
  }

  /// Compute the center point of the range
  /// If range is empty(), 'nan' will be returned
  double center() const
  {
    if (isEmpty()) return nan("");
    return min_ + (max_ - min_) / 2.0;
  }

  /// Get the 'width' of the range
  /// If range is empty(), 'nan' will be returned
  double getSpan() const
  {
    if (isEmpty()) return nan("");
    return max_ - min_;
  }

  bool operator==(const RangeBase& rhs) const
  {
    return min_ == rhs.min_ && max_ == rhs.max_;
  }

  /**
   * \brief Return the current range, or (if empty) a full range (-1e308, 1e308).
   * \return A range where always: min <= max
   */
  std::pair<double, double> getNonEmptyRange() const
  {
    // pair with full range
    if (isEmpty()) return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()};
    else
      return {min_, max_};
  }

protected:
  // make members non-accessible to maintain invariant: min <= max  (unless uninitialized)
  double min_ = std::numeric_limits<double>::max();
  double max_ = std::numeric_limits<double>::lowest();
};

OPENMS_DLLAPI std::ostream& operator<<(std::ostream& out, const RangeBase& b);

struct OPENMS_DLLAPI RangeRT : public RangeBase
{

  const static MSDim DIM = MSDim::RT;

  // Rule of 0!
  using RangeBase::RangeBase; // inherit C'tors from base
  using RangeBase::operator=; // inherit assignment operator from base

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

  // Rule of 0!
  using RangeBase::RangeBase; // inherit C'tors from base
  using RangeBase::operator=; // inherit assignment operator

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

struct OPENMS_DLLAPI RangeIntensity : public RangeBase
{

  const static MSDim DIM = MSDim::INT;

  // Rule of 0!
  using RangeBase::RangeBase; // inherit C'tors from base
  using RangeBase::operator=; // inherit assignment operator

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

  // Rule of 0!
  using RangeBase::RangeBase; // inherit C'tors from base
  using RangeBase::operator=; // inherit assignment operator

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
  ALL,  ///< all dimensions are filled
  SOME, ///< some dimensions are empty, some are filled
  NONE  ///< all dimensions are empty (=cleared)
};

/**
  @brief Handles the management of a multidimensional range, e.g. RangeMZ and RangeIntensity for spectra.

  Instantiate it with the dimensions which are supported/required, e.g.
  `RangeManager<RangeRT, RangeMZ> range_spec` for a spectrum and use the strongly typed features, such as
  `range_spec.getMaxRT()/setMaxRT(500.0)` or `range_spec.extend(RangeMZ{100, 1500})`.

  Use RangeManagerContainer as a base class for all peak and feature containers like MSSpectrum, MSExperiment and FeatureMap.

  The implementation uses non-virtual multiple inheritance using variadic templates. Each dimension, e.g. RangeRT, is inherited from, thus
  all members of the base class become accessible in the RangeManager, e.g. getMaxRT().
  Operations (e.g. assignment, or extension of ranges) across RangeManagers with a different, yet overlapping set of base classes
  is enabled using fold expressions and constexpr evaluations, which are resolved at compile time (see for_each_base_ member function).

*/
template<typename... RangeBases>
class RangeManager : public RangeBases...
{
public:
  using ThisRangeType = RangeManager<RangeBases...>;

  // rule of 0 -- no need for a virtual d'tor or anything fancy
  // ...

  bool operator==(const RangeManager& rhs) const
  {
    bool equal = true;
    for_each_base_([&](auto* base) {
      using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      equal &= ((T_BASE&)rhs == (T_BASE&)*this);
    });
    return equal;
  }

  bool operator!=(const RangeManager& rhs) const
  {
    return ! operator==(rhs);
  }

  /// copy all overlapping dimensions from @p rhs to this instance.
  /// Dimensions which are not contained in @p rhs are left untouched.
  /// @param rhs Range to copy from
  /// @return true if one or more dimensions overlapped
  template<typename... RangeBasesOther>
  bool assignUnsafe(const RangeManager<RangeBasesOther...>& rhs)
  {
    bool found = false;
    for_each_base_([&](auto* base) {
      using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
      {
        base->operator=((T_BASE&)rhs);
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
  auto& assign(const RangeManager<RangeBasesOther...>& rhs)
  {
    if (! assignUnsafe(rhs))
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);");
    }
    return *this;
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
        base->extend((T_BASE&)rhs);
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
    if (! extendUnsafe(rhs))
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);");
    }
  }

  /// calls RangeBase::scale() for each dimension
  void scaleBy(const double factor)
  {
    for_each_base_([&](auto* base) { base->scaleBy(factor); });
  }

  /**
   * \brief If any dimension is a single point, e.g. min==max, then extend this dimension by @p min_span / 2 on either side.
   *
   * Empty dimensions remain unchanged.
   *
   * @see DimBase::minSpanIfSingular
   *
   */
  void minSpanIfSingular(const double min_span)
  {
    for_each_base_([&](auto* base) { base->minSpanIfSingular(min_span); });
  }


  /// Move range of *this to min/max of @p rhs, without changing the span, if possible.
  /// This does tighten the range unless @p rhs's ranges are smaller than *this.
  /// Dimensions which are not contained in @p rhs or are empty are left untouched.
  /// @param rhs Range to translate/move the current range into
  /// @return true if dimensions overlapped, false otherwise
  template<typename... RangeBasesOther>
  bool pushIntoUnsafe(const RangeManager<RangeBasesOther...>& rhs)
  {
    bool found = false;
    for_each_base_([&](auto* base) {
      using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
      {
        const auto& rhs_base = (T_BASE&)rhs;
        if (! rhs_base.isEmpty()) base->pushInto(rhs_base);
        found = true;
      }
    });
    return found;
  }

  /// Move range of *this to min/max of @p sandbox, without changing the span, if possible.
  /// This does tighten the range unless @p sandbox's ranges are smaller than *this.
  /// Dimensions which are not contained in @p sandbox or are empty are left untouched.
  /// @param sandbox Range to translate/move the current range into
  /// @throw Exception::InvalidRange if no dimensions overlapped
  template<typename... RangeBasesOther>
  void pushInto(const RangeManager<RangeBasesOther...>& sandbox)
  {
    if (! pushIntoUnsafe(sandbox))
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);");
    }
  }


  /// Clamp min/max of all overlapping dimensions to min/max of @p rhs
  /// Dimensions which are not contained in @p rhs or where rhs is empty are left untouched.
  /// @param rhs Range to clamp to
  /// @return true if dimensions overlapped, false otherwise
  template<typename... RangeBasesOther>
  bool clampToUnsafe(const RangeManager<RangeBasesOther...>& rhs)
  {
    bool found = false;
    for_each_base_([&](auto* base) {
      using T_BASE = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      if constexpr (std::is_base_of_v<T_BASE, RangeManager<RangeBasesOther...>>)
      {
        const auto& rhs_base = (T_BASE&)rhs;
        if (! rhs_base.isEmpty()) base->clampTo(rhs_base);
        found = true;
      }
    });
    return found;
  }

  /// Clamp min/max of all overlapping dimensions to min/max of @p rhs.
  /// This may tighten the range (even to a single point).
  /// Dimensions which are not contained in @p rhs or where rhs is empty are left untouched.
  /// @param rhs Range to clamp to
  /// @throw Exception::InvalidRange if no dimensions overlapped
  template<typename... RangeBasesOther>
  void clampTo(const RangeManager<RangeBasesOther...>& rhs)
  {
    if (! clampToUnsafe(rhs))
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);");
    }
  }

  /// obtain a range dimension at runtime using @p dim
  const RangeBase& getRangeForDim(MSDim dim) const
  {
    RangeBase* r_base = nullptr;

    static_for_each_base_([&](auto* base) {
      using Base = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      if (base->DIM == dim) r_base = (Base*)this;
    });

    assert((r_base != nullptr) && "No base class has this MSDim!");
    return *r_base;
  }

  /// obtain a range dimension at runtime using @p dim
  RangeBase& getRangeForDim(MSDim dim)
  {
    RangeBase* r_base = nullptr;

    static_for_each_base_([&](auto* base) {
      using Base = std::decay_t<decltype(*base)>; // remove const/ref qualifiers
      if (base->DIM == dim) r_base = (Base*)this;
    });

    assert((r_base != nullptr) && "No base class has this MSDim!");
    return *r_base;
  }

  /// is any/some/all dimension in this range populated?
  HasRangeType hasRange() const
  {
    constexpr size_t total {sizeof...(RangeBases)}; // total number of bases
    size_t count {0};
    for_each_base_([&](auto* base) { count += ! base->isEmpty(); });
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
        if (base->contains((T_BASE&)rhs)) return;
        contained = false;
      }
    });
    if (! has_overlap) throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

    return contained;
  }

  /// Resets all ranges
  void clearRanges()
  {
    for_each_base_([&](auto* base) { base->clear(); });
  }

  /// Resets the dimension of the given @p range. Any type of ion mobility in @p range will clear the Mobility dimension.
  /// If the @p range is not contained in this class, then nothing happens.
  ThisRangeType& clear(const DIM_UNIT range)
  {
    switch (range)
    {
      case DIM_UNIT::RT:
        if constexpr (std::is_base_of_v<RangeRT, ThisRangeType>) this->RangeRT::clear();
        break;
      case DIM_UNIT::MZ:
        if constexpr (std::is_base_of_v<RangeMZ, ThisRangeType>) this->RangeMZ::clear();
        break;
      case DIM_UNIT::INT:
        if constexpr (std::is_base_of_v<RangeIntensity, ThisRangeType>) this->RangeIntensity::clear();
        break;
      // assume all ion mobility ranges are the same and never occur together. If this is violated at some point, then split RangeMobility into
      // subclasses...
      case DIM_UNIT::IM_MS:
      case DIM_UNIT::IM_VSSC:
      case DIM_UNIT::FAIMS_CV:
        if constexpr (std::is_base_of_v<RangeMobility, ThisRangeType>) this->RangeMobility::clear();
        break;
      default:
        // all cases should be covered above
        assert(false && "This should never be reached. Did you forget to implement a new DIM_UNIT?");
    }
    return *this;
  }

  /// print each dimension (base classes) to a stream
  void printRange(std::ostream& out) const
  {
    for_each_base_([&](auto* base) { out << *base; });
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
template<typename... RangeBases>
class RangeManagerContainer : public RangeManager<RangeBases...>
{
public:
  using ThisRangeType = typename RangeManager<RangeBases...>::ThisRangeType;

  /// D'tor
  virtual ~RangeManagerContainer() = default; // required since we have virtual methods

  /// implement this function to reflect the underlying data of the derived class (e.g. an MSSpectrum)
  /// Usually, call clearRanges() internally and then populate the dimensions.
  virtual void updateRanges() = 0;

  /// get range of current data (call updateRanges() before to ensure the range is accurate)
  const ThisRangeType& getRange() const
  {
    return (ThisRangeType&)*this;
  }

  /// get mutable range, provided for efficiency reasons (avoid updateRanges(), if only minor changes were made)
  ThisRangeType& getRange()
  {
    return (ThisRangeType&)*this;
  }
};

/// Range which contains all known dimensions
using RangeAllType = RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>;

} // namespace OpenMS
