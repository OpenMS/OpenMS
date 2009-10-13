// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_TYPES_H
#define OPENMS_CONCEPT_TYPES_H

#include <OpenMS/config.h>

#include <limits>
#include <cstddef> // for size_t
#include <ctime>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

// If possible use the ISO C99-compliant header stdint.h
// to define the portable integer types.
#ifdef OPENMS_HAS_STDINT_H
#include <stdint.h>
#endif

namespace OpenMS
{
	/**
		@brief Signed integer type (32bit)

		@ingroup Concept
  */
  typedef OPENMS_INT32_TYPE Int32;

	/**
		@brief Signed integer type (64bit)

		@ingroup Concept
  */
  typedef OPENMS_INT64_TYPE Int64;

	/**
		@brief Unsigned integer type (64bit)

		@ingroup Concept
  */
	typedef OPENMS_UINT64_TYPE UInt64;

	/**
		@brief Time type

		Use this type to represent a point in time (as a synonym for time_t).

		@ingroup Concept
	*/
	typedef time_t 	Time;

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
		@brief Real type

		Use this type to represent standard floating point numbers.

		@ingroup Concept
	*/
	typedef float Real;

	/**
		@brief Double-precision real type

		Use this type to represent double precision floating point numbers.

		@ingroup Concept
	*/
	typedef double DoubleReal;


	/**
		@brief Byte type

		Use this type to represent byte data (8 bit length). A Byte is always unsigned.

		@ingroup Concept
	*/
	typedef	OPENMS_BYTE_TYPE Byte;

	/**
		@brief A unique object ID (as unsigned 64bit type).

		@see PersistentObject

		@ingroup Concept
	*/
	typedef OPENMS_UINT64_TYPE UID;

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

  DoubleReal d = NUMBER;
  std::cout.precision(writtenDigits<DoubleReal>()); // explicit template instantiation
  std::cout << writtenDigits<DoubleReal>() << ": " << d << '\n'; // writes: 15: 12345.6789012346

  Real r = NUMBER;
  std::cout.precision(writtenDigits(r)); // type deduced from argument
  std::cout << writtenDigits(r) << ": " << r << '\n'; // writes: 6: 12345.7

  long double l = NUMBER;
  std::cout.precision(writtenDigits(1L)); // argument is not used, but L suffix indicates a long double
  std::cout << writtenDigits(1L) << ": " << l << '\n'; // writes: 18: 12345.6789012345671

  DoubleReal x = 88.99;
  std::cout.precision(15);
  std::cout << "15: " << x << '\n'; // writes: 15: 88.99
  std::cout.precision(16);
  std::cout << "16: " << x << '\n'; // writes: 16: 88.98999999999999
	@endcode
	*/
	//@{


	/**@brief Number of digits commonly used for writing a floating point type
  (a.k.a. precision).  Specializations are defined for float, double, long
  double.
	*/
 	template <typename FloatingPointType> inline Int writtenDigits(const FloatingPointType& /* unused */ = FloatingPointType());

	/// Number of digits commonly used for writing a @c float (a.k.a. precision).
	template <> inline Int writtenDigits <float> (const float & )
	{
		return std::numeric_limits<float>::digits10;
	}

	/// Number of digits commonly used for writing a @c double (a.k.a. precision).
	template <> inline Int writtenDigits <double> (const double & )
	{
		return std::numeric_limits<double>::digits10;
	}

	/// We do not want to bother people who unintentionally provide an int argument to this.
	template <> inline Int writtenDigits <int> (const int & )
	{
		return std::numeric_limits<int>::digits10;
	}

	/// We do not want to bother people who unintentionally provide an unsigned int argument to this.
	template <> inline Int writtenDigits <unsigned int> (const unsigned int & )
	{
		return std::numeric_limits<unsigned int>::digits10;
	}

	/// We do not want to bother people who unintentionally provide a long int argument to this.
	template <> inline Int writtenDigits <long int> (const long int & )
	{
		return std::numeric_limits<int>::digits10;
	}

	/// We do not want to bother people who unintentionally provide an unsigned long int argument to this.
	template <> inline Int writtenDigits <unsigned long int> (const unsigned long int & )
	{
		return std::numeric_limits<unsigned int>::digits10;
	}

	class DataValue;
	/// DataValue will be printed like double.
	template <> inline Int writtenDigits <DataValue> (const DataValue & )
	{
		return std::numeric_limits<double>::digits10;
	}

	/*
	META-COMMENT:  DO NOT INTRODUCE ANY LINEBREAKS BELOW IN
	"<code>std::numeric_limits<long double>::digits10 == 18</code>".
	The doxygen parser (version 1.5.5) will get confused!  (Clemens)
	*/

	/**@brief Number of digits commonly used for writing a @c long @c double (a.k.a. precision). ...

  Note: On Microsoft platforms, the I/O sytem seems to treat @c long @c double
  just like @c double.  We observed that
	<code>std::numeric_limits<long double>::digits10 == 18</code>
	with GCC 3.4 on MinGW, but this promise is
  <i>not</i> kept by the Microsoft I/O sytem libraries.  Therefore we use the
  value of @c digits10 for @c double also for @c long @c double.  See
  http://msdn.microsoft.com/ + search: "long double".
	*/
	template <> inline Int writtenDigits <long double> (const long double & )
	{
#ifndef OPENMS_WINDOWSPLATFORM
		return std::numeric_limits<long double>::digits10;
#else
		return std::numeric_limits<double>::digits10;
#endif
	}

	/// The general template definition will force a compile-time error if FloatingPointType is in fact not a floating point type.  Only the template specializations for float, double, long double shall be used.
 	template <typename FloatingPointType> inline Int writtenDigits(const FloatingPointType& /* unused */)
	{
		// Self-explanatory compile time error!
		return FloatingPointType::Sorry_but_writtenDigits_is_designed_to_work_for_floating_point_types_only;
	}

	// Note: I once tried to move PrecisionWrapper to namespace Internal, but oops! operator <<  won't be found (through ADL?) anymore.
	/// Wrapper class to implement output with appropriate precision.  See precisionWrapper().
	template <typename FloatingPointType >
	struct PrecisionWrapper
	{
		/// Constructor.  Note: Normally you will prefer to use the "make"-function precisionWrapper(), which see.
		PrecisionWrapper(const FloatingPointType rhs) : ref_(rhs) {}
		PrecisionWrapper(const PrecisionWrapper& rhs) : ref_(rhs.ref_) {}
		FloatingPointType const ref_;
	 private:
		PrecisionWrapper(); // intentionally not implemented
	};

	/**@brief Wrapper function that sets the appropriate precision for output
	temporarily.  The original precision is restored afterwards so that no side
	effects remain.  This is a "make"-function that deduces the typename
	FloatingPointType from its argument and returns a
	PrecisionWrapper<FloatingPointType>.

	Example:
	@code
	std::cout
	<< 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
	<< precisionWrapper(0.1234567890123456789f) << '\n' // float
	<< 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
	<< precisionWrapper(0.1234567890123456789) << '\n' // double
	<< 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
	<< precisionWrapper(0.1234567890123456789l) << '\n' // long double
	<< 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n';
	@endcode
	Result:
	@code
	0.123457 0.123457 0.123457
	0.123457
	0.123457 0.123457 0.123457
	0.123456789012346
	0.123457 0.123457 0.123457
	0.123456789012345679
	0.123457 0.123457 0.123457
	@endcode

	Note: Unfortunately we cannot return a const& - this will change when rvalue
	references become part of the new C++ standard.  In the meantime, we need a
	copy constructor for PrecisionWrapper.
  */
	template <typename FloatingPointType>
	inline const PrecisionWrapper<FloatingPointType> precisionWrapper(const FloatingPointType rhs)
	{
		return PrecisionWrapper<FloatingPointType>(rhs);
	}

	/// Output operator for a PrecisionWrapper.  Specializations are defined for float, double, long double.
	template <typename FloatingPointType >
	inline std::ostream & operator << ( std::ostream& os, const PrecisionWrapper<FloatingPointType> &rhs)
	{
	  // Same test as used by isnan(), spelled out here to avoid issues during overload resolution.
	  if ( rhs.ref_ != rhs.ref_ )
	  {
	    // That's what Linux GCC uses, and gnuplot understands.
	    // Windows would print stuff like 1.#QNAN which makes testing hard.
	    return os << "nan";
	  }
	  else
	  {
	    const std::streamsize prec_save = os.precision();
	    return os << std::setprecision(writtenDigits(FloatingPointType()))
	    << rhs.ref_ << std::setprecision(prec_save);
	  }
	}

	//@}

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

	Example (Tutorial_typeAsString.C):
	@dontinclude Tutorial_typeAsString.C
	@until end of Tutorial_typeAsString.C
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
	template < typename Type >
	std::string typeAsString(const Type & /* unused */ = Type() )
	{
#ifndef OPENMS_COMPILER_GXX
		return "[ Sorry, OpenMS::typeAsString() relies upon GNU extension __PRETTY_FUNCTION__ ]";
#else
		std::string pretty(__PRETTY_FUNCTION__);
		static char const context_left [] = "with Type =";
		static char const context_right [] = "]";
		size_t left = pretty.find(context_left);
		left+=sizeof(context_left);
		size_t right = pretty.rfind(context_right);
		if ( right <= left ) return pretty; // oops!
		return pretty.substr(left, right-left);
#endif
	}
	
	namespace Internal
  {
    /** Used to set the locale to "C", to avoid
        problems on machines with incompatible
        locale settings (this overwrites the
        locale setting of the environment!)
    */
    static const char* OpenMS_locale = setlocale(LC_ALL, "C");
  }

} // namespace OpenMS

#endif // OPENMS_CONCEPT_TYPES_H
