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
    @brief A group of IndexTuple.
    
    This class stores a set of IndexTuple.
  */
  template < typename ContainerType >
  class Group : public std::set
          < IndexTuple< ContainerType >, typename IndexTuple< ContainerType >::IndexLess >
    {
    public:
      typedef std::set
        < IndexTuple< ContainerType >, typename IndexTuple< ContainerType >::IndexLess > Base;

      /** @name STL-compliance type definitions of the container interface*/
      //@{
      typedef typename Base::iterator iterator;
      typedef typename Base::const_iterator const_iterator;
      typedef typename Base::reverse_iterator reverse_iterator;
      typedef typename Base::const_reverse_iterator const_reverse_iterator;
      typedef typename Base::value_type value_type;
      typedef typename Base::reference reference;
      typedef typename Base::const_reference const_reference;
      typedef typename Base::pointer pointer;
      typedef typename Base::difference_type difference_type;
      typedef typename Base::size_type size_type;
      //@}

      typedef IndexTuple< ContainerType > Element;

      /// Default constructor
      Group() : Base()
    {}

      /// Copy constructor
      inline Group(const Group& source)
          : Base(source)
      {}

      /// Assignment operator
      Group& operator= (const Group& source)
      {
        if (&source == this)
          return *this;

        Base::operator=(source);
        return *this;
      }

      /// Destructor
      virtual ~Group()
    {}

      /// Returns the number of elements
      inline const unsigned int count() const
      {
        return this->size();
      }

      /// Returns wheter the group is empty
      inline bool isEmpty()
      {
        return (this->size() == 0);
      }

      /// Inserts an element into the group
      std::pair< typename Base::iterator, bool >
      insert(const Element& elem) throw (Exception::InvalidValue)
      {
        std::pair< typename Base::iterator, bool > pr;
        pr = Base::insert(elem);

        if (pr.second == false)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",(String)(elem.getMapIndex())) ;
        }
        else
        {
          return pr;
        }
      }

      /// Equality operator
      virtual bool operator == (const Group& group) const
      {
        return std::operator==(group,*this);
      }

      /// Equality operator
      virtual bool operator != (const Group& group) const
      {
        return !(std::operator==(group,*this));
      }

    public:
      /// Serialization interface
      //     template<class Archive>
      //     void serialize(Archive & ar, const unsigned int /* version */ )
      //     {
      //       ar & boost::serialization::make_nvp("group",boost::serialization::base_object<Group<ElementContainerType> >(*this));
      //     }
      //@}

      /// Serialization
      //     friend class boost::serialization::access;

    };

  ///Print the contents to a stream.
  template < typename ContainerT >
  std::ostream& operator << (std::ostream& os, const Group<ContainerT>& cons)
  {
    os << "---------- GROUP BEGIN -----------------\n";
    unsigned int i=0;
    for (typename Group<ContainerT>::const_iterator it = cons.begin(); it != cons.end(); ++it,i)
    {
      os  << "Element: " << i << '\n'
      << it << std::endl;
    }
    return os;
  }
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_GROUP_H
