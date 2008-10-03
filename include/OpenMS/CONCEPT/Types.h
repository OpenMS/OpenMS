// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_TYPES_H
#define OPENMS_CONCEPT_TYPES_H

#include <OpenMS/config.h>

#include <limits>
#include <time.h>
#include <string>

#ifdef OPENMS_HAS_BASETSD_H
#include <basetsd.h>
#endif

// If possible use the ISO C99-compliant header stdint.h
// to define the portable integer types.
#ifdef OPENMS_HAS_STDINT_H
#include <stdint.h>
#endif

// Added to avoid warnings with MS Visual Studio .NET
#ifdef OPENMS_COMPILER_MSVC
#	pragma warning( disable : 4290 )
#endif

namespace OpenMS
{
	#ifndef OPENMS_HAS_STDINT_H

	/**
		@brief Signed integer type
		
		@ingroup Concept
  */
	typedef OPENMS_INDEX_TYPE Int;

	/**
		@brief Unsigned integer type
		
		@ingroup Concept
  */
	typedef OPENMS_SIZE_TYPE UInt;

	/**	
		@brief Time type
		
		Use this type to represent a point in time (as a replaecement for time_t).
		
		@ingroup Concept
	*/
	typedef time_t 	Time;

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
	typedef	unsigned char Byte;

	/**	
		@brief Pointer-sized unsigned int type
		
		This type holds unsigned 64 bit integer numbers and is used to store pointers
		in a portable fashion (for both 32-bit and 64-bit systems).
	
		@ingroup Concept
	*/
	typedef OPENMS_ULONG64_TYPE PointerSizeUInt;

	/**	
		@brief Pointer-sized signed int type
		
		This type holds signed 64 bit numbers and is used to store pointers
		in a portable fashion (for both 32-bit and 64-bit systems).
	
		@ingroup Concept
	*/
	typedef OPENMS_LONG64_TYPE PointerSizeInt;

	/**	
		@brief A unique object ID.
		
		@see PersistentObject
		
		@ingroup Concept
	*/
	typedef OPENMS_ULONG64_TYPE UID;

	/**	
		@brief Offset type (64-bit).
		
		This pointer holds signed 64 bit integers for large file access
		
		@ingroup Concept
	*/
	typedef OPENMS_LONG64_TYPE Offset64Int;
	
	
	#else
   	// the ISO C99 definitions
		typedef int32_t Int; 
		typedef uint32_t UInt;
		typedef time_t	Time;
   	typedef unsigned int UInt;
		typedef float Real;
		typedef double DoubleReal;
		typedef	uint8_t Byte;
		typedef int64_t	 PointerSizeInt;
		typedef uint64_t PointerSizeUInt;
		typedef uint64_t UID;
    typedef int64_t	 Offset64Int;

	#endif
  
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

	These static const integeres are provided to unify the handling of this
	issue throughout %OpenMS.  (So please don't use ad-hoc numbers ;-) )

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
  std::cout.precision(written_digits_doublereal);
  std::cout << written_digits_doublereal << ": " << d << '\n'; // writes: 15: 12345.6789012346
  
  Real r = NUMBER;
  std::cout.precision(written_digits_real);
  std::cout << written_digits_real << ": " << r << '\n'; // writes: 6: 12345.7
  
  long double l = NUMBER;
  std::cout.precision(written_digits_long_double);
  std::cout << written_digits_long_double << ": " << l << '\n'; // writes: 18: 12345.6789012345671

  DoubleReal x = 88.99;
  std::cout.precision(15);
  std::cout << "15: " << x << '\n'; // writes: 15: 88.99
  std::cout.precision(16);
  std::cout << "16: " << x << '\n'; // writes: 16: 88.98999999999999
	@endcode
	*/
	//@{
	/// Number of digits commonly used for writing a @c double (a.k.a. precision).
	static const Int written_digits_double = std::numeric_limits<double>::digits10;
	/// Number of digits commonly used for writing a @c DoubleReal (a.k.a. precision).
	static const Int written_digits_doublereal = std::numeric_limits<DoubleReal>::digits10;
	/// Number of digits commonly used for writing a @c float (a.k.a. precision).
	static const Int written_digits_float = std::numeric_limits<float>::digits10;
	/// Number of digits commonly used for writing a @c Real (a.k.a. precision).
	static const Int written_digits_real = std::numeric_limits<Real>::digits10;
	/// Number of digits commonly used for writing a @c long @c double (a.k.a. precision).
	static const Int written_digits_long_double = std::numeric_limits<long double>::digits10;
	//@}

	/**
	@brief Returns the @c Type as as std::string.

	Have you ever spent a long time trying to find out what a @c typedef
	actually "points" to?  Then this can help.

	typeAsString is implemented as a function template.  There are two ways to us this:
	@code
	SomeType instance;
	string what_type_1 = typeAsString(instance);
	string what_type_2 = typeAsString<SomeType>();
	@endcode
	The %typeAsString<SomeType>() version seems to go a bit deeper.
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
	
} // namespace OpenMS

#endif // OPENMS_CONCEPT_TYPES_H
