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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H
#define OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H

#include <iostream>
#include <vector>

#include<OpenMS/KERNEL/DPosition.h>


namespace OpenMS
{
  /**
    @brief 
            
    @ingroup 
  */

  template < typename ContainerType >
  class IndexTuple
  {
    public:

      /**
        @name Type definitions
      */
      //@{
      typedef typename ContainerType::value_type ElementType;
      typedef typename ElementType::TraitsType TraitsType;
      typedef DPosition<2, TraitsType> PositionType;
      //@}

      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      IndexTuple()
          : map_index_(0),
          feature_index_(0),
          element_pointer_(0)
      {}

      /// Constructor
      inline IndexTuple(const UnsignedInt& map_index, const UnsignedInt& element_index, const ElementType& element)
      {
        map_index_ = map_index;
        feature_index_ = element_index;
        element_pointer_ = const_cast<ElementType*>(&element);
        transformed_position_ = element_pointer_->getPosition();
      }

      /// Copy constructor
      inline IndexTuple(const IndexTuple& source)
      {
        map_index_ = source.map_index_;
        feature_index_ = source.feature_index_;
        element_pointer_ = source.element_pointer_;
        transformed_position_ = source.transformed_position_;
      }

      /// Assignment operator
      IndexTuple& operator = (const IndexTuple& source)
      {
        if (&source == this)
          return *this;

        map_index_ = source.map_index_;
        feature_index_ = source.feature_index_;
        element_pointer_ = source.element_pointer_;
        transformed_position_ = source.transformed_position_;
        return *this;
      }

      /// Destructor
      virtual ~IndexTuple()
    {}
      //@}

      /** @name Accessors */
      //@{
      /// Non-mutable access to the container index
      inline const UnsignedInt& getMapIndex() const
      {
        return map_index_;
      }
      /// mutable access to the container index
      inline UnsignedInt& getMapIndex()
      {
        return map_index_;
      }
      /// Set the container index
      inline void setMapIndex(const UnsignedInt& c)
      {
        map_index_ = c;
      }

      /// Non-mutable access to the element index
      inline const UnsignedInt& getElementIndex() const
      {
        return feature_index_;
      }
      /// Set the container
      inline void setElementIndex(const UnsignedInt& e)
      {
        feature_index_= e;
      }

      /// Non-mutable access to the element
      inline const ElementType& getElement() const
      {
        return *element_pointer_;
      }
      /// Set the container
      inline void setElement(const ElementType& e)
      {
        element_pointer_ = &e;
      }

      /// Non-mutable access to the transformed position
      inline const PositionType& getTransformedPosition() const
      {
        return transformed_position_;
      }
      /// Mutable access to the transformed position
      inline PositionType& getTransformedPosition()
      {
        return transformed_position_;
      }
      /// Set the container
      inline void setTransformedPosition(const PositionType& p)
      {
        transformed_position_ = p;
      }
      //@}


      /// Equality operator
      virtual bool operator == (const IndexTuple& i) const
      {
        return ((map_index_ == i.map_index_) && (feature_index_ == i.feature_index_) && (element_pointer_ == i.element_pointer_));
      }

      /// Equality operator
      virtual bool operator != (const IndexTuple& i) const
      {
        return !((map_index_ == i.map_index_) && (feature_index_ == i.feature_index_) && (element_pointer_ == i.element_pointer_));
      }

      /// Compare by getOverallQuality()
      struct IndexLess
            : std::binary_function < IndexTuple, IndexTuple, bool >
      {
        bool operator () ( IndexTuple const & left, IndexTuple const & right ) const
        {
          return ( left.map_index_ < right.map_index_ );
        }
      };

    protected:
      PositionType transformed_position_;
      UnsignedInt map_index_;
      UnsignedInt feature_index_;
      ElementType* element_pointer_;
  };

  ///Print the contents to a stream.
  template < typename ContainerT >
  std::ostream& operator << (std::ostream& os, const IndexTuple<ContainerT>& cons)
  {
    os << "---------- IndexTuple -----------------\n"
    << "Transformed Position: " << cons.getTransformedPosition() << '\n'
    << "Original Position: " << cons.getElement() << '\n'
    << "Index: " << cons.getElementIndex() << std::endl;
    return os;
  }


} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMAPPING_INDEXTUPLE_H
