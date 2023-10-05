// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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


