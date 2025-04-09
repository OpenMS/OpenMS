// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <ctime>
#include <cstddef> // for size_t & ptrdiff_t
#include <limits>

#include <cstdint> // since C++11

namespace OpenMS
{
  /**
    @brief Signed integer type (32bit)

    @ingroup Concept
  */
  typedef int32_t Int32;

  /**
    @brief Unsigned integer type (32bit)

    @ingroup Concept
  */
  typedef uint32_t UInt32;

  /**
    @brief Signed integer type (64bit)

    @ingroup Concept
  */
  typedef int64_t Int64;

  /**
    @brief Unsigned integer type (64bit)

    @ingroup Concept
  */
  typedef uint64_t UInt64;

  /**
    @brief Time type

    Use this type to represent a point in time (as a synonym for time_t).

    @ingroup Concept
  */
  typedef time_t  Time;

  /**
    @brief Unsigned integer type

    @ingroup Concept
  */
  //typedef size_t UInt;
  typedef unsigned int UInt;

  /**
    @brief Signed integer type

    @ingroup Concept
  */
  //typedef OPENMS_SIZE_T_SIGNED Int;
  typedef int Int;

  /**
    @brief Byte type

    Use this type to represent byte data (8 bit length). A Byte is always unsigned.

    @ingroup Concept
  */
  typedef uint8_t Byte;

  /**
    @brief A unique object ID (as unsigned 64bit type).

    @see PersistentObject

    @ingroup Concept
  */
  typedef uint64_t UID;

  /**
    @brief Size type e.g. used as variable which can hold result of size()

    @ingroup Concept
  */
  typedef size_t Size;

  /**
    @brief Signed Size type e.g. used as pointer difference

    @ingroup Concept
  */
  typedef ptrdiff_t SignedSize;

  enum ASCII
  {
    ASCII__BACKSPACE        = '\b',
    ASCII__BELL             = '\a',
    ASCII__CARRIAGE_RETURN  = '\r',
    ASCII__HORIZONTAL_TAB   = '\t',
    ASCII__NEWLINE          = '\n',
    ASCII__RETURN           = ASCII__NEWLINE,
    ASCII__SPACE            = ' ',
    ASCII__TAB              = ASCII__HORIZONTAL_TAB,
    ASCII__VERTICAL_TAB     = '\v',

    ASCII__COLON            = ':',
    ASCII__COMMA            = ',',
    ASCII__EXCLAMATION_MARK = '!',
    ASCII__POINT            = '.',
    ASCII__QUESTION_MARK    = '?',
    ASCII__SEMICOLON        = ';'
  };

  //@}

  /**
    @name Numbers of digits used for writing floating point numbers (a.k.a. precision).

    These functions are provided to unify the handling of this issue throughout
    %OpenMS.  (So please don't use ad-hoc numbers ;-) )

    If you want to avoid side effects you can use precisionWrapper() to write a
    floating point number with appropriate precision; in this case the original
    state of the stream is automatically restored afterwards.  See
    precisionWrapper() for details.

    In practice, the number of decimal digits that the type can represent
    without loss of precision are 6 digits for single precision
    and 15 digits for double precision.
    We have \f$2^{24}/10^{6}=16.777216\f$ and \f$2^{53}/10^{15}=9.007199254740992\f$,
    so rounding will remove the remaining difference.

    Example:
    @code
    #define NUMBER 12345.67890123456789012345678901
    std::cout << NUMBER << '\n'; // default precision, writes: 12345.7

    double d = NUMBER;
    std::cout.precision(writtenDigits<double>(0.0)); // explicit template instantiation
    std::cout << writtenDigits<double>(0.0) << ": " << d << '\n'; // writes: 15: 12345.6789012346

    float r = NUMBER;
    std::cout.precision(writtenDigits(r)); // type deduced from argument
    std::cout << writtenDigits(r) << ": " << r << '\n'; // writes: 6: 12345.7

    long double l = NUMBER;
    std::cout.precision(writtenDigits(1L)); // argument is not used, but L suffix indicates a long double
    std::cout << writtenDigits(1L) << ": " << l << '\n'; // writes: 18: 12345.6789012345671

    double x = 88.99;
    std::cout.precision(15);
    std::cout << "15: " << x << '\n'; // writes: 15: 88.99
    std::cout.precision(16);
    std::cout << "16: " << x << '\n'; // writes: 16: 88.98999999999999
    @endcode
  */
  //@{


  /**
    @brief Number of digits commonly used for writing a floating point type
    (a.k.a. precision).  Specializations are defined for float, double, long
    double.
  */
  template <typename FloatingPointType>
  inline constexpr Int writtenDigits(const FloatingPointType& /* unused */ = FloatingPointType());

  /// Number of digits commonly used for writing a @c float (a.k.a. precision).
  template <>
  inline constexpr Int writtenDigits<float>(const float&)
  {
    return std::numeric_limits<float>::digits10;
  }

  /// Number of digits commonly used for writing a @c double (a.k.a. precision).
  template <>
  inline constexpr Int writtenDigits<double>(const double&)
  {
    return std::numeric_limits<double>::digits10;
  }

  /// We do not want to bother people who unintentionally provide an int argument to this.
  template <>
  inline constexpr Int writtenDigits<int>(const int&)
  {
    return std::numeric_limits<int>::digits10;
  }

  /// We do not want to bother people who unintentionally provide an unsigned int argument to this.
  template <>
  inline constexpr Int writtenDigits<unsigned int>(const unsigned int&)
  {
    return std::numeric_limits<unsigned int>::digits10;
  }

  /// We do not want to bother people who unintentionally provide a long int argument to this.
  template <>
  inline constexpr Int writtenDigits<long int>(const long int&)
  {
    return std::numeric_limits<int>::digits10;
  }

  /// We do not want to bother people who unintentionally provide an unsigned long int argument to this.
  template <>
  inline constexpr Int writtenDigits<unsigned long int>(const unsigned long int&)
  {
    return std::numeric_limits<unsigned int>::digits10;
  }

  class DataValue;
  /// DataValue will be printed like double.
  template <>
  inline constexpr Int writtenDigits<DataValue>(const DataValue&)
  {
    return std::numeric_limits<double>::digits10;
  }

  /*
    META-COMMENT:  DO NOT INTRODUCE ANY LINEBREAKS BELOW IN
    "<code>std::numeric_limits<long double>::digits10 == 18</code>".
    The doxygen parser (version 1.5.5) will get confused!  (Clemens)
  */

  /**
    @brief Number of digits commonly used for writing a @c long @c double (a.k.a. precision). ...

    Note: On Microsoft platforms, the I/O system seems to treat @c long @c double
    just like @c double.  We observed that
    <code>std::numeric_limits<long double>::digits10 == 18</code>
    with GCC 3.4 on MinGW, but this promise is
    <i>not</i> kept by the Microsoft I/O system libraries.  Therefore we use the
    value of @c digits10 for @c double also for @c long @c double.  See
    http://msdn.microsoft.com/ + search: "long double".
  */
  template <>
  inline constexpr Int writtenDigits<long double>(const long double&)
  {
#ifndef OPENMS_WINDOWSPLATFORM
    return std::numeric_limits<long double>::digits10;

#else
    return std::numeric_limits<double>::digits10;

#endif
  }

  /**
   @brief The general template definition will return the default precision of
   6 according to 27.4.4.1 basic_iosconstructors (C++ Standard).
   */
  template <typename FloatingPointType>
  inline constexpr Int writtenDigits(const FloatingPointType& /* unused */)
  {
    return 6;
  }

  namespace Internal
  {
    /**
      Used to set the locale to "C", to avoid
      problems on machines with incompatible
      locale settings (this overwrites the
      locale setting of the environment!)
    */
    extern OPENMS_DLLAPI const char* OpenMS_locale;
  }

} // namespace OpenMS

