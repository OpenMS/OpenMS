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
// $Maintainer: Hannes Roest $
// $Authors: Philippe M. Groarke, Hannes Roest $
// --------------------------------------------------------------------------

// See https://philippegroarke.com/posts/2018/easy_defensive_programming/
// https://github.com/p-groarke/defensive_cpp

#pragma once
#include <type_traits>

namespace OpenMS
{
  // no constexpr lamdas in C++11, therefore we have to use functions
  namespace Test
  {
#if __GNUC__ < 5
    // older gcc, will not compile
    bool fulfills_rule_of_5() {return true;}
    bool fulfills_rule_of_6() {return true;}
    bool fulfills_fast_vector() {return true;}
#else
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
#endif // check for gcc < 5.0

  }
}

