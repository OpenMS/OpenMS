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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H
#define OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H

#include <iostream>
#include <vector>

#include <OpenMS/KERNEL/RawDataPoint2D.h>

namespace OpenMS
{
  /**
    @brief This class stores 2D peak/feature representations.
    
    It is used to reference to peaks or features in different maps.
    
    The actual position and the intensity are stores in the base class RawDataPoint2D.
    The original datapoint is referenced by the map and element index.
    
    @todo rename (Marc, Clemens)
  */
  class IndexTuple
  	: public RawDataPoint2D
  {
  	
    public:
    
      /// Default constructor
      inline IndexTuple()
      	: RawDataPoint2D(),
        	map_index_(0),
          element_index_(0)
      {
      }

      /// Constructor with map index, element index and position
      inline IndexTuple(UInt map_index, UInt element_index, const RawDataPoint2D& point)
      {
      	this->RawDataPoint2D::operator=(point);
        map_index_ = map_index;
        element_index_ = element_index;
      }

      /// Copy constructor
      inline IndexTuple(const IndexTuple& rhs)
      	: RawDataPoint2D(rhs),
      		map_index_(rhs.map_index_),
      		element_index_(rhs.element_index_)
      {
      }

      /// Assignment operator
      inline IndexTuple& operator = (const IndexTuple& rhs)
      {
        if (&rhs == this) return *this;
        
        RawDataPoint2D::operator=(rhs);
        map_index_ = rhs.map_index_;
        element_index_ = rhs.element_index_;

        return *this;
      }

      /// Destructor
      virtual ~IndexTuple()
      {
      }
      
      /// Returns the map index
      inline UInt getMapIndex() const
      {
        return map_index_;
      }
      
      /// Set the map index
      inline void setMapIndex(UInt i)
      {
        map_index_ = i;
      }

      /// Returns the element index
      inline UInt getElementIndex() const
      {
        return element_index_;
      }
      
      /// Set the element index
      inline void setElementIndex(UInt e)
      {
        element_index_= e;
      }

      /// Equality operator
      virtual bool operator == (const IndexTuple& i) const
      {
        return ((map_index_ == i.map_index_) && (element_index_ == i.element_index_) && (intensity_ == i.intensity_));
      }

      /// Equality operator
      virtual bool operator != (const IndexTuple& i) const
      {
        return !((map_index_ == i.map_index_) && (element_index_ == i.element_index_) && (intensity_ == i.intensity_));
      }
			
			///Comparator by map and element index
      struct IndexLess
      	: std::binary_function < IndexTuple, IndexTuple, bool >
      {
        inline bool operator () ( IndexTuple const & left, IndexTuple const & right ) const
        {
        	//if map indices are equal, use element indices
        	if ( left.map_index_ == right.map_index_)
        	{
        		return left.element_index_ < right.element_index_;
          }
          //else use map indices
          return ( left.map_index_ < right.map_index_ );
        }
      };

    protected:
    	
      /// Int of the element's container
      UInt map_index_;
      /// Int of the element within element's container
      UInt element_index_;
  };

  ///Print the contents of an IndexTuple to a stream.
  std::ostream& operator << (std::ostream& os, const IndexTuple& cons);
  	
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMAPPING_INDEXTUPLE_H
