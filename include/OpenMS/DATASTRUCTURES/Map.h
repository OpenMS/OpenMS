// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MAP_H
#define OPENMS_DATASTRUCTURES_MAP_H

#include <OpenMS/CONCEPT/Exception.h>

#include <map>

namespace OpenMS
{
	/**
		@brief Map class based on the STL map (containing serveral convenience functions)

		@ingroup Datastructures
	*/
	template <class Key, class T>
	class Map
	  : public std::map<Key,T>
	{
		public:

			/**
				@brief Map illegal key exception
	  
        Thrown when trying to access an element with operator[], which is not contained in the Map
        , i.e. no default ctor is called (as done in std::map)

				@ingroup Exceptions
			*/
			class IllegalKey
				:	public Exception::BaseException
			{
				public:
				IllegalKey(const char* file, int line, const char* function)
					:	Exception::BaseException(file, line, function)
				{
				}
			};
			
			///@name OpenMS style typedefs
			//@{
			typedef std::map<Key,T> Base;
			typedef typename Base::value_type ValueType;
			typedef Key KeyType;
			typedef typename Base::value_type* PointerType;
			typedef typename Base::iterator Iterator;
			typedef typename Base::const_iterator ConstIterator;
			typedef typename Base::reverse_iterator ReverseIterator;
			typedef typename Base::const_reverse_iterator ConstReverseIterator;
			//@}

			///Test whether the map contains the given key.
			inline bool has(const Key& key) const
			{
				return Base::find(key)!=Base::end();
			}

			/**	
				@brief Return a constant reference to the element whose key is @p key.
				
				@exception IllegalKey if the given key does not exist
			*/
			const T& operator [] (const Key& key) const;

			/// Return a mutable reference to the element whose key is @p key. If an element with the key @p key does not exist, it is inserted.
			T& operator [] (const Key& key);
	};
	
	//******************************************************************************************
	// Implementations of template methods
	//******************************************************************************************
	
	template <class Key, class T>
	const T& Map<Key, T>::operator [] (const Key& key) const
	{
		ConstIterator it = this->find(key);
		if (it == Base::end())
		{
			throw IllegalKey(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		else
		{
			return it->second;
		}
	}

	template <class Key, class T>
	T& Map<Key, T>::operator [] (const Key& key)
		
	{
		Iterator it = this->find(key);
		if (it == Base::end())
		{
			it = this->insert(ValueType(key, T())).first;
		}
		return it->second;
	}

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_MAP_H
