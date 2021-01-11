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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS
{
  /**
    @brief Returns the @c Type as as std::string.

    Have you ever spent a long time trying to find out what a @c typedef
    actually "points" to?  Then this can help.

    typeAsString is implemented as a function template.  There are two ways to us this:
    @code
      SomeType instance;
      string what_type_1 = typeAsString(instance);
      string what_type_2 = typeAsString< SomeType >();
    @endcode

    The %typeAsString< SomeType >() version seems to go a bit deeper.
    Sometimes the results
    depend on how the %typeAsString() is instantiated in the first place.
    The argument given to the function is never used, it only serves to infer the type.
    You can even supply function pointers, etc.

    Example (Tutorial_typeAsString.cpp):
    @dontinclude Tutorial_typeAsString.cpp
    @until end of Tutorial_typeAsString.cpp
    On a 64 bit platform running GCC 4.3.1, this produced the following output:
    @code
    int
    unsigned int
    double
    float

    int
    long unsigned int

    OpenMS::Peak1D
    OpenMS::Peak1D
    OpenMS::DPosition<1u>
    double
    float

    double ()(int, int*)
    WOW<const char* const*** const&, 5>
    Oink<double, 55, 666u, WOW>
    float ()(float&)
    double (WOW<char, 8>::*)(const double&)
    @endcode
  */
  template <typename Type>
  std::string typeAsString(const Type & /* unused */ = Type())
  {
#if ! defined(OPENMS_PRETTY_FUNCTION)
    return "[ Sorry, OpenMS::typeAsString() relies upon extension OPENMS_PRETTY_FUNCTION ]";
#else
    std::string pretty(OPENMS_PRETTY_FUNCTION);
    static char const context_left[] = "with Type =";
    static char const context_right[] = "]";
    size_t left = pretty.find(context_left);
    left += sizeof(context_left);
    size_t right = pretty.rfind(context_right);
    if (right <= left || right == std::string::npos || left == std::string::npos)
      return pretty;                      // different format as expected

    return pretty.substr(left, right - left);

#endif
  }

} // namespace OpenMS


