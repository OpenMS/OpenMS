// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: RTTI.h,v 1.3 2006/03/10 11:06:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_RTTI_H
#define OPENMS_FORMAT_RTTI_H

#include <OpenMS/CONCEPT/Types.h>
#include <string>
#include <typeinfo>

using std::string;

namespace OpenMS 
{

#	ifdef	__GNUC__
	// EGCS produces a nearly unreadable name mangling that requires 
	// further interpretation
	namespace GNUDemangling 
	{
		string demangle(string s);
	}
#	endif

	/**	
		@brief Returns a unique name for a class.
			
		This name contains no blanks. It is usually derived by substituting all
		blanks in the name (as returned by RTTI::getName()) with underscores ("_").
		In the case of <tt>gcc</tt>, however a name demangling decodes the string first.
		This function is needed for object persistence.
		@param	t the <tt>type_info</tt> structure as returned by <tt>typeid</tt>
		@return	string the coverted class name
		
		@ingroup Concept
		
		@todo complete tests
	*/ 
	string streamClassName(const std::type_info& t);

	/**	
		@brief Simplified RunTime Type Identification.
			
		ANSI C++ provides support for runtime type identification (RTTI). However, the support
		is very basic. The template functions of the RTTI namespace  provide a more 
		readable support for RTTI. It defines
		predicates such as  \link OpenMS::RTTI::isKindOf isKindOf \endlink that simplify tests on the hereditary relationship of
		different objects.
		
		To use the RTTI template functions, parametrize it with the type you are interested in.
		For example, to find out whether a given DRawDataPoint<1> is a DPeak<1>, the following code
		can be used:
		
		@code
			DRawDataPoint<1>& r =...;
			...
			if (RTTI::isKindOf< DPeak<1> >(r)) 
			{
				// perform some peak specific operations
			} 
			else 
			{
				// this is only a raw data point...
			}
		@endcode

		@ingroup Concept
	*/
	namespace RTTI
	{

		/**	
			@brief Return a reference to a static default instance of the corresponding class.
				
			This method is basically intended to provide a default object for certain operations
			that require an instance of a certain class without really using this instance.
			It is mainly used inside the RTTI class.
		*/
		template <typename T>
		const T& getDefault() 
		{
			static T t;
			return t;
		}

		/**	
			@brief Return a void pointer to a new instance of the class.
			
			Use this method to provide an easy factory for objects of a certain class.
			The main use of this function lies in object persistence. The PersistenceManager
			needs a function for the dynamic creation of objects.
		*/
		template <typename T>
		void* getNew()
		{
			return static_cast<void*>(new T);
		}

		/**	
			@brief Return the name of the class.
			
			This method returns the name of the class as given by <tt>typeid(\<class instance\>.name())</tt>.
			No additional name demangling and whitespace substitution are performed.
		*/
		template <typename T>
		const char* getName()
		{
			return typeid(getDefault<T>()).name();
		}

		/**	
			@brief Return a void pointer that is unique for each class.
		*/
		template <typename T>
		void* getClassID()
		{
			static char dummy;
			return (void*)&dummy;
		}

		/**	
			@brief Return the demangled class name.
			
			The class name is demangled (as far as possible) and in the
			resulting string blanks are substituted by underscores, so the
			name can be read from a stream as one string.
			The typical usage is something like
			@code
				String class_name = RTTI::getStreamName< DPeak<1> >();
				...
			@endcode
		*/
		template <typename T>
		const char* getStreamName()
		{
			// define portable names for the portable
			// types (some platforms use Size, some unsigned int, 
			// SUN CC even unsigned  for the Size type)
			if ((typeid(T) == typeid(Size)) 
					|| (typeid(T) == typeid(Position))
					|| (typeid(T) == typeid(HashIndex))
					|| (typeid(T) == typeid(Property))
					|| (typeid(T) == typeid(Handle)))
			{
				return "OpenMS::Size";
			}
			if ((typeid(T) == typeid(Index))
					|| (typeid(T) == typeid(ErrorCode))
					|| (typeid(T) == typeid(Distance)))
			{
				return "OpenMS::Index";
			}
			if (typeid(T) == typeid(::std::string))
			{
				return "::std::string";
			}
			if (typeid(T) == typeid(PointerSizeInt))
			{
				return "OpenMS::PointerSizeInt";
			}
			if (typeid(T) == typeid(bool))
			{
				return "bool";
			}
			if (typeid(T) == typeid(float))
			{
				return "float";
			}
			if (typeid(T) == typeid(char))
			{
				return "char";
			}
			if (typeid(T) == typeid(unsigned char))
			{
				return "unsigned_char";
			}
			if (typeid(T) == typeid(double))
			{
				return "double";
			}
			static string s("");
			static bool is_set = false;

			if (!is_set)
			{
				is_set = true;
				s = streamClassName(typeid(getDefault<T>()));
			}

			return s.c_str();
		}

		/**	
			@brief Return true if <tt>u</tt> is a kind of T.
			
			If <tt>u</tt> is an instance of a class derived from T,
			this predicate returns true:
			
			@code
				DPeak<1> p;

				// return true, since DPeak is derived from DRawDataPoint
				bool is_rawdatapoint = RTTI::isKindOf< DRawDataPoint<1> >(p);
			@endcode
		*/
		template <typename T, typename U>
		bool isKindOf(const U&  u)
		{
			return (0 != dynamic_cast<const T*>(&u));
		}

		/**	
			@brief Cast an object to a specialized class.
				
			<b>Example:</b>
			@code
				DRawDataPoint<1>* rawdp = ...;
				...
			
				// check whether the DRawDataPoint is also an DPeak
				if (RTTI::isKindOf< DPeak<1> >(rawdp))
				{
					// perform some peak specific actions
					DPeak<1>* peak = RTTI::castTo< DPeak<1> >(*rawdp);
					...
				}
			@endcode
		*/
		template <typename T, typename U>
		T* castTo(const U& u)
		{
			return const_cast<T*>(dynamic_cast<const T*>(&u));
		}

		/**	
			@brief Return <b>true</b> if a given object is an instance of a given class.
				
			If <tt>u</tt> is an instance of <tt>T</tt>, this predicate returns <b>true</b>.
			If <tt>u</tt> is an instance of a class that is derived from <tt>T</tt> or
			a base class of <tt>T</tt>, it returns false.
		*/
		template <typename T, typename U>
		bool isInstanceOf(const U& u)
		{
			T		t;
			return (typeid(u) == typeid(t));
		}
	} // namespace RTTI

} // namespace OpenMS

#endif // OPENMS_FORMAT_RTTI_H
