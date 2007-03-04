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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_HASHFUNCTION_H
#define OPENMS_CONCEPT_HASHFUNCTION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>

namespace OpenMS
{
	/**	@brief General Hash Function Template.
	 
			This template function provides a simple wrapper
			for the specialized hash functions. It facilitates their use 
			in STL hash associative containers which expect a <b>Hasher</b>
			class as template parameter.
			
			@ingroup Concept
	*/
  template <typename T>
  class HashFunction
  {
    public:
		
		UInt operator () (const T& t) const throw()
    {
      return Hash(t);
		}
	};

	/**	@name Specialized Hash Functions.
	*/
	//@{
	
	/**
	*/
  extern UInt hashPointer(void *const ptr) throw();

	/**
	*/
  extern UInt hashString(const char* str) throw();

	/**
	*/
  extern UInt hashPJWString(const char* str) throw();

	/**
	*/
  extern UInt hashElfString(const char* str) throw();

	/** General default hash function.
			This method converts a given key to a  UInt by calling <tt>(UInt)key</tt>.
			If the key type <tt>T</tt> is not convertible to UInt by
			default, a converter should be defined (<tt>operator UInt</tt>).
			@param	key the key to be hashed
			@return	UInt the hash index
	*/
	template <typename T>
	inline UInt Hash(const T& key) throw()
	{
		return static_cast<UInt>((OPENMS_POINTERSIZEINT_TYPE)key);
	}

	/** String hash function.
			This method is optimized for the hashing of OpenMS Strings.
			In fact, it is only an inline wrapper around hashString.
	*/
	inline UInt Hash(const String& s) throw()
	{
		return hashString(s.c_str());
	}

	/** string hash function.
	  	This method is optimized for the hashing of STL strings.
	  	In fact, it is only an inline wrapper aound hashString.
	*/
	inline UInt Hash(const std::string& s) throw()
	{
		return hashString(s.c_str());
	}

	/** Pointer hash function.
			Use this function to hash pointers to objects.
	*/
	inline UInt Hash(void* const& ptr) throw()
	{
		return hashPointer(ptr);
	}

	//@}
		
	//@{
	/**	Calculate the next prime number.
			This method returns the first prime number that is 
			greater or equal to the number given as the argument.
			Only odd prime numbers are returned, the lowest number returned is 3.
	*/
	UInt getNextPrime(UInt l) throw();

	//@}

} // namespace OPENMS

#endif // OPENMS_CONCEPT_HASHFUNCTION_H
