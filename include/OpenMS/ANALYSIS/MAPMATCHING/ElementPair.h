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

#include<OpenMS/KERNEL/Feature.h>
// #include <OpenMS/CONCEPT/Types.h>

#include <utility> // std::pair

namespace OpenMS
{

/**
	 @brief A pair of element in two different maps.

   The purpose of the mapmatching stage is to identify pairs of features in
   different map, to estimate a transformation that maps features in a
   specified range onto each other and to execute this transform
   (dewarping).
*/
template < typename ElementType = Feature >
class ElementPair : public std::pair<ElementType,ElementType>
{
public:
    typedef std::pair<ElementType,ElementType> Base;
  	typedef DoubleReal QualityType;

    /** @name Constructors and Destructor
     */
    //@{
    /// Default constructor
    ElementPair() : Base(), quality_(0)
    {}

    /// Copy constructor
    ElementPair(const ElementPair& fp)
            : Base(fp), quality_(fp.quality_)
    {}

    ElementPair(ElementType const & first, ElementType const & second, QualityType const & quality = QualityType(0))
            : Base(first,second), quality_(quality)
    {}

    /// Destructor
    virtual ~ElementPair()
    {}
    //@}

    /// assignment operator
    ElementPair& operator = (const ElementPair& rhs)
    {
        if (&rhs==this)
            return *this;

        Base::operator = (rhs);
        quality_       = rhs.quality_;

        return *this;
    }

    bool operator == (const ElementPair& rhs) const
    {
        return ( this->getFirst()   == rhs.getFirst() &&
                 this->getSecond()  == rhs.getSecond() &&
                 this->getQuality() == rhs.getQuality() );
    }

    bool operator != (const ElementPair& rhs) const
    {
        return  !(*this == rhs);
    }

    /**@name Accessors
     */
    //@{
    /// Non-mutable access to the first feature
    const ElementType& getFirst() const
    {
        return this->first;
    }
    /// Mutable access to the first feature
    ElementType& getFirst()
    {
        return this->first;
    }
    /// Non-mutable access to the first feature
    void setFirst(const ElementType& frt)
    {
        this->first = frt;
    }

    /// Non-mutable access to the second feature
    const ElementType& getSecond() const
    {
        return this->second;
    }
    /// Mutable access to the second feature
    ElementType& getSecond()
    {
        return this->second;
    }
    /// Non-mutable access to the second feature
    void setSecond(const ElementType& sec)
    {
        this->second = sec;
    }

    /// Non-mutable access to the quality of the pair
    QualityType getQuality() const
    {
        return quality_;
    }
    /// Mutable access to the quality of the pair
    QualityType& getQuality()
    {
        return quality_;
    }
    /// Mutable access to the quality of the pair
    void setQuality(QualityType ql)
    {
        quality_ = ql;
    }
    //@}

protected:

    /// quality of the pair (not individual features)
    QualityType quality_;

}
; // end of class ElementPair

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_ELEMENTPAIR_H
