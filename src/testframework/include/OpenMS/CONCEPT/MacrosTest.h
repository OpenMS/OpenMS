// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Philippe M. Groarke, Hannes Roest $
// --------------------------------------------------------------------------

// See https://philippegroarke.com/posts/2018/easy_defensive_programming/
// https://github.com/p-groarke/defensive_cpp

#pragma once
#include <type_traits>

namespace OpenMS
{
  // no constexpr lambdas in C++11, therefore we have to use functions
  namespace Test
  {
    template <class T>
    constexpr bool fulfills_rule_of_5()
    {
      static_assert(std::is_destructible<T>::value, "T : must be destructible");
      static_assert(std::is_copy_constructible<T>::value, "T : must be copy constructible");
      static_assert(std::is_move_constructible<T>::value, "T : must be move constructible");
      static_assert(std::is_copy_assignable<T>::value, "T : must be copy assignable");
      static_assert(std::is_move_assignable<T>::value, "T : must be move assignable");

      return std::is_destructible<T>::value &&
             std::is_copy_constructible<T>::value &&
             std::is_move_constructible<T>::value &&
             std::is_copy_assignable<T>::value &&
             std::is_move_assignable<T>::value;
    }

    template <class T>
    constexpr bool fulfills_rule_of_6()
    {
      static_assert(fulfills_rule_of_5<T>(), "T : must fulfill rule of 5");
      static_assert(std::is_default_constructible<T>::value, "T : must be default constructible");

      return fulfills_rule_of_5<T>() &&
             std::is_default_constructible<T>::value;
    }

    template <class T>
    constexpr bool fulfills_fast_vector()
    {
      static_assert( (std::is_trivially_copy_constructible<T>::value && std::is_trivially_destructible<T>::value) ||
                      std::is_nothrow_move_constructible<T>::value,
                      "T : doesn't fulfill fast vector (trivially copy constructible " \
                        "and trivially destructible, or nothrow move constructible)");

      return (std::is_trivially_copy_constructible<T>::value && std::is_trivially_destructible<T>::value) ||
              std::is_nothrow_move_constructible<T>::value;
    };
  }
}

