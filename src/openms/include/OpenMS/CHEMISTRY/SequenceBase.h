// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <type_traits>

namespace OpenMS
{

  /**
   * @brief CRTP base class providing common interface for sequence types
   *
   * This template base class uses the Curiously Recurring Template Pattern (CRTP)
   * to provide a common interface for both AASequence and NASequence classes.
   *
   * @tparam Derived The derived sequence class (AASequence or NASequence)
   *
   * @ingroup Chemistry
   */
  template<typename Derived>
  class SequenceBase
  {
  public:
    /// Get the derived instance (CRTP pattern)
    const Derived& derived() const noexcept
    {
      return static_cast<const Derived&>(*this);
    }

    /// Get the derived instance (CRTP pattern)
    Derived& derived() noexcept
    {
      return static_cast<Derived&>(*this);
    }

    /// Check if sequence is empty
    bool empty() const
    {
      return derived().empty();
    }

    /// Get the size of the sequence
    Size size() const
    {
      return derived().size();
    }

    /// It's not feasible to support full functionality for all sequence types, without templating the fragment types
    /// instead we just expose the "Full" fragment type for both AASequence and NASequence

    /// Get monoisotopic weight (implemented in cpp file for each specialization)
    double getMonoWeight(Int charge = 0) const;
    
    /// Get average weight (implemented in cpp file for each specialization)
    double getAverageWeight(Int charge = 0) const;
    
    /// Get empirical formula (implemented in cpp file for each specialization)
    EmpiricalFormula getFormula(Int charge = 0) const;

    /// Get string representation
    String toString() const
    {
      return derived().toString();
    }

    /// Get prefix subsequence
    Derived getPrefix(Size length) const
    {
      return derived().getPrefix(length);
    }

    /// Get suffix subsequence
    Derived getSuffix(Size length) const
    {
      return derived().getSuffix(length);
    }

    /// Get subsequence
    Derived getSubsequence(Size start, Size length) const
    {
      return derived().getSubsequence(start, length);
    }

    /// Generic sequence processing function
    template<typename Func>
    auto processSequence(Func&& func) const -> decltype(func(derived()))
    {
      return func(derived());
    }

    /// Generic sequence processing function (mutable)
    template<typename Func>
    auto processSequence(Func&& func) -> decltype(func(derived()))
    {
      return func(derived());
    }

  protected:
    /// Protected constructor to prevent direct instantiation
    SequenceBase() = default;
    
    /// Protected destructor (CRTP pattern)
    ~SequenceBase() = default;
  };

  /**
   * @brief Type trait to check if a type is a sequence type
   */
  template<typename T>
  struct is_sequence_type : std::false_type {};

  template<typename Derived>
  struct is_sequence_type<SequenceBase<Derived>> : std::true_type {};

  template<typename T>
  inline constexpr bool is_sequence_type_v = is_sequence_type<T>::value;

} // namespace OpenMS