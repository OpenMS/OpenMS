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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SERIALIZATION_H
#define OPENMS_FORMAT_SERIALIZATION_H

/**
  @file Serialization.h
  @addtogroup Serialization

	OpenMS provides (partial) support for the Boost serialization library.
	Classes which can be serialized include Serialization.h.  If you need
	more header files, add them to Serialization.h.  Currently not all classes
	in OpenMS can be serialized, but those in KERNEL can.  Serialization support
	is enabled if you "#define OPENMS_HAVE_SERIALIZATION 1" before
	Serialization.h is included.  See the documentation of
  #OPENMS_HAVE_SERIALIZATION about that.
*/

#ifndef OPENMS_HAVE_SERIALIZATION
/**@brief Enable support for boost::serialization
 
	We recognized that including the header files for boost serialization
	increases the compilation time.  Thus we provide a mechanism to compile
	things using dumb forward declarations, without including the boost header
	files, unless we are really using boost serialization.  (Their documentation
	should point to this macro.)  Serialization headers are not included by
	default.	If you want to serialize (more precisely, instantiate the "real" serialize()
	function templates), write:
	
	<code>\#define OPENMS_HAVE_SERIALIZATION 1</code>
 
	before any OpenMS header files are included.

 */
#define OPENMS_HAVE_SERIALIZATION 0
#endif

///////////////////////////////////////////////////////////////////////////////
#if OPENMS_HAVE_SERIALIZATION

// For splitting serialize() into save() and load()
#include <boost/serialization/split_member.hpp>

// Provides name-value-pairs:  make_nvp("bla",bla_);
#include <boost/serialization/nvp.hpp>

// Serialization of base classes of derived classes
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/void_cast.hpp>

// Serialization of standard types.
// All of STL can be serialized.
// Add headers as necessary.
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

///////////////////////////////////////////////////////////////////////////////
#else // OPENMS_HAVE_SERIALIZATION == false 

/// Dummy, see #OPENMS_HAVE_SERIALIZATION
namespace boost
{
	/// Dummy, see #OPENMS_HAVE_SERIALIZATION
	namespace serialization
	{

	/// Dummy, see #OPENMS_HAVE_SERIALIZATION
	class access;

#ifdef BOOST_SERIALIZATION_SPLIT_MEMBER
#error BOOST_SERIALIZATION_SPLIT_MEMBER is already #defined -- probably you are doing serialization but forgot to "#define OPENMS_HAVE_SERIALIZATION 1"
#else

/// Dummy, see #OPENMS_HAVE_SERIALIZATION
#define	BOOST_SERIALIZATION_SPLIT_MEMBER()

#endif

		/// Dummy, see #OPENMS_HAVE_SERIALIZATION
		template<class T> struct nvp;

		/// Dummy, see #OPENMS_HAVE_SERIALIZATION
		template<class T>
		inline const nvp<T> make_nvp(const char * name, T & t);

		/// Dummy, see #OPENMS_HAVE_SERIALIZATION.
		template< /* class Base, */ class Derived>
		Derived & base_object(Derived &d);
		
	}
}

#endif // OPENMS_HAVE_SERIALIZATION
///////////////////////////////////////////////////////////////////////////////




namespace OpenMS
{
	/**@brief A little helper so that we can serialize <i>mutable</i> stuff.
	 
	 Note that this will bypass a compile-time trap in Boost Serialization.  <i>Use
	 at your own risk!</i>  The dangers of serializing mutable stuff are explained in
	 "Boost documentation" > "Serialization" > "Rationale" > "Compile time trap
	 when saving a non-const value".
 
	 The hack is provided here because programmers will work around the trap
	 anyway -- so this will be something to grep for.  ;-)

	 (Note that <code>makeConstReference(x)</code> is equivalent to
	 <code>boost::cref(x).get()</code>.)
	*/
	template < typename T >
	const T & makeConstReference ( const T & t )
	{
		return t;
	}

} // namespace OpenMS

#endif // OPENMS_FORMAT_SERIALIZATION_H
