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
          : map_pointer_(0),
          element_pointer_(0)
      {}

      /// Constructor
      inline IndexTuple(const ContainerType& container, const ElementType& element)
      {
        map_pointer_ = const_cast<ContainerType*>(&container);
        element_pointer_ = const_cast<ElementType*>(&element);
        transformed_position_ = element_pointer_->getPosition();
      }

      /// Copy constructor
      inline IndexTuple(const IndexTuple& source)
      {
        map_pointer_ = source.map_pointer_;
        element_pointer_ = source.element_pointer_;
        transformed_position_ = source.transformed_position_;
      }

      /// Assignment operator
      IndexTuple& operator = (const IndexTuple& source)
      {
        if (&source == this)
          return *this;

        map_pointer_ = source.map_pointer_;
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
      /// Non-mutable access to the container
      inline const ContainerType& getContainer() const
      {
        return *map_pointer_;
      }
      /// Set the container
      inline void setContainer(const ContainerType& c)
      {
        map_pointer_ = &c;
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
        return ((map_pointer_ == i.map_pointer_) && (element_pointer_ == i.element_pointer_));
      }

      /// Equality operator
      virtual bool operator != (const IndexTuple& i) const
      {
        return !((map_pointer_ == i.map_pointer_) && (element_pointer_ == i.element_pointer_));
      }

      /// Compare by getOverallQuality()
      struct MapAdressLesser
            : std::binary_function < IndexTuple, IndexTuple, bool >
      {
        bool operator () ( IndexTuple const & left, IndexTuple const & right ) const
        {
          return ( left.map_pointer_ < right.map_pointer_ );
        }
      };

    protected:
      PositionType transformed_position_;
      ContainerType* map_pointer_;
      ElementType* element_pointer_;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMAPPING_INDEXTUPLE_H
