// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <limits.h>
#include <time.h>

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
	/**
		@defgroup Concept Concept
		
		@brief OpenMS concepts
	*/

	/**
		@defgroup BasicTypes Basic data types
		
		@brief Basic data types.
		
		@ingroup Concept
		
		@{
	*/
	
	#ifndef OPENMS_HAS_STDINT_H

	/// Signed integer type
	typedef OPENMS_INDEX_TYPE SignedInt;

	/// Unsigned integer type
	typedef OPENMS_SIZE_TYPE UnsignedInt;

	/**	
		@brief Time type
		
		Use this type to represent a point in time (as a replaecement for time_t).
	*/
	typedef time_t 	Time;

	/**	
		@brief Real type
		
		Use this type to represent standard floating point numbers.
	*/
	typedef float Real;

	/**	
		@brief Double-precision real type
		
		Use this type to represent double precision floating point numbers.
	*/
	typedef double DoubleReal;

	/**	
		@brief Byte type
		
		Use this type to represent byte data (8 bit length). A Byte is always unsigned.
	*/
	typedef	unsigned char Byte;

	/**	
		@brief Pointer-sized unsigned int type
		
		This type holds unsigned 64 bit integer numbers and is used to store pointers
		in a portable fashion (for both 32-bit and 64-bit systems).
	*/
	typedef OPENMS_POINTERSIZEUINT_TYPE PointerSizeUInt;

	/**	
		@brief Pointer-sized signed int type
		
		This type holds unsigned 64 bit numbers and is used to store pointers
		in a portable fashion (for both 32-bit and 64-bit systems).
	*/
	typedef OPENMS_POINTERSIZEINT_TYPE PointerSizeInt;

	/**	
		@brief A unique object ID.
		
		@see PersistentObject
	*/
	typedef OPENMS_ULONG64_TYPE UID;

/** @} */ // end of BasicTypes


	#else
   	// the ISO C99 definitions
		typedef int32_t SignedInt; 
		typedef uint32_t UnsignedInt;
		typedef time_t	Time;
   	typedef Position UnsignedInt;
		typedef float Real;
		typedef double DoubleReal;
		typedef	uint8_t Byte;
		typedef int64_t	 PointerSizeInt;
		typedef uint64_t PointerSizeUInt;
		typedef uint64_t UID;

	#endif

	//@}

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
	
} // namespace OpenMS

#endif // OPENMS_CONCEPT_TYPES_H
