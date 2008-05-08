// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_GROUP_H
#define OPENMS_ANALYSIS_MAPMATCHING_GROUP_H

#include <OpenMS/ANALYSIS/MAPMATCHING/IndexTuple.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>

namespace OpenMS
{
  /**
    @brief A set of IndexTuples.
    
    @todo rename (Marc, Clemens)
  */
  class Group 
    : public std::set < IndexTuple, IndexTuple::IndexLess >
  {
    	
    public:
    
      typedef std::set< IndexTuple, IndexTuple::IndexLess > Base;

      /// Default constructor
      inline Group() 
       : Base()
    	{
    	}

      /// Copy constructor
      inline Group(const Group& source)
        : Base(source)
      {
      }

      /// Assignment operator
      inline Group& operator= (const Group& source)
      {
        if (&source == this) return *this;
        Base::operator=(source);
        return *this;
      }

      /// Destructor
      virtual ~Group()
    	{
    	}

      /**
      	@brief Inserts an element into the group
      	
      	@exception InvalidValue is thrown if the element is already contained
      */
      inline Base::iterator insert(const IndexTuple& elem) throw (Exception::InvalidValue)
      {
        std::pair< Base::iterator, bool > pr = Base::insert(elem);

        if (pr.second == false)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",(String)(elem.getMapIndex())) ;
        }
        
        return pr.first;
      }

      /// Equality operator
      inline bool operator == (const Group& group) const
      {
        return std::operator==(group,*this);
      }

      /// Equality operator
      inline bool operator != (const Group& group) const
      {
        return std::operator!=(group,*this);
      }
	};
    
  ///Print the contents of a Group to a stream
  std::ostream& operator<<(std::ostream& os, const Group& cons);

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_GROUP_H
