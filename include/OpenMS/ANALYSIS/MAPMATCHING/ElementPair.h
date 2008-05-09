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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_ELEMENTPAIR_H
#define OPENMS_ANALYSIS_MAPMATCHING_ELEMENTPAIR_H

#include <OpenMS/CONCEPT/Types.h>

#include <utility> // std::pair

namespace OpenMS
{

	/**
		 @brief A pair of elements that has a quality store
		 
		 @todo make this a struct or derive it non-public from the std::pair (Marc, Clemens)
	*/
	template < typename ElementType >
	class ElementPair
	  : public std::pair<ElementType,ElementType>
	{
		
		public:
	  	///Base class
	    typedef std::pair<ElementType,ElementType> Base;
	
	    /** @name Constructors and Destructor
	     */
	    //@{
	    /// Default constructor
	    ElementPair() 
	      : Base(), 
	      	quality_(0)
	    {
	    }
	
	    /// Copy constructor
	    ElementPair(const ElementPair& rhs)
	      : Base(rhs), 
	      	quality_(rhs.quality_)
	    {
	    }
			
			///Constructor with two elements and a quality
	    ElementPair(const ElementType&first, const ElementType& second, DoubleReal quality =0.0)
	      : Base(first,second),
	      	quality_(quality)
	    {
	    }
	
	    /// Destructor
	    virtual ~ElementPair()
	    {
	    }
	    //@}
	
	    /// assignment operator
	    ElementPair& operator=(const ElementPair& rhs)
	    {
	      if (&rhs==this) return *this;
	
	      Base::operator = (rhs);
	      quality_ = rhs.quality_;
	
	      return *this;
	    }
			
			///Equality operator
	    bool operator==(const ElementPair& rhs) const
	    {
	      return ( this->getFirst()   == rhs.getFirst() &&
	               this->getSecond()  == rhs.getSecond() &&
	               this->getQuality() == rhs.getQuality() );
	    }
			
			///Equality operator
	    bool operator!=(const ElementPair& rhs) const
	    {
				return ( this->getFirst()   != rhs.getFirst() ||
	               this->getSecond()  != rhs.getSecond() ||
	               this->getQuality() != rhs.getQuality() );
	    }
	
	    /**@name Accessors
	     */
	    //@{
	    /// Non-mutable access to the first element
	    const ElementType& getFirst() const
	    {
	    	return this->first;
	    }
	    
	    /// Sets the first element
	    void setFirst(const ElementType& element)
	    {
	    	this->first = element;
	    }

	    /// Non-mutable access to the second element
	    const ElementType& getSecond() const
	    {
	    	return this->second;
	    }

	    /// Sets the second element
	    void setSecond(const ElementType& element)
	    {
	    	this->second = element;
	    }
	
	    /// Non-mutable access to the quality of the pair
	    DoubleReal getQuality() const
	    {
	    	return quality_;
	    }

	    /// Mutable access to the quality of the pair
	    void setQuality(DoubleReal quality)
	    {
	    	quality_ = quality;
	    }
	    //@}
	
		protected:
	
	    /// quality of the pair (not individual elements)
	    DoubleReal quality_;
	
	}; // end of class ElementPair

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_ELEMENTPAIR_H
