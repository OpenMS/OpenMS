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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DFEATUREPAIR_H
#define OPENMS_ANALYSIS_MAPMATCHING_DFEATUREPAIR_H

#include<OpenMS/KERNEL/DFeature.h>

#include <OpenMS/CONCEPT/Types.h>

#include <utility> // std::pair

namespace OpenMS
{

/**
	 @brief A pair of features in two different maps.

   The purpose of the mapmatching stage is to identify pairs of features in
   different map, to estimate a transformation that maps features in a
   specified range onto each other and to execute this transform
   (dewarping).
*/
template <Size D, typename FeatureT = DFeature<D> >
class DFeaturePair : public std::pair<FeatureT,FeatureT>
{
public:

    enum { DIMENSION = D };
    typedef FeatureT FeatureType;
    typedef std::pair<FeatureType,FeatureType> Base;
    // CHANGED
    typedef	 typename FeatureType::TraitsType::QualityType QualityType;

    /** @name Constructors and Destructor
     */
    //@{
    /// Default constructor
    DFeaturePair() : Base(), quality_(0)
    {}

    /// Copy constructor
    DFeaturePair(const DFeaturePair& fp)
            : Base(fp), quality_(fp.quality_)
    {}

    DFeaturePair(FeatureType const & first, FeatureType const & second, QualityType const & quality = QualityType(0))
            : Base(first,second), quality_(quality)
    {}

    /// Destructor
    virtual ~DFeaturePair()
    {}
    //@}

    /// assignment operator
    DFeaturePair& operator = (const DFeaturePair& rhs)
    {
        if (&rhs==this)
            return *this;

        Base::operator = (rhs);
        quality_       = rhs.quality_;

        return *this;
    }

    bool operator == (const DFeaturePair& rhs) const
    {
        return ( this->getFirst()   == rhs.getFirst() &&
                 this->getSecond()  == rhs.getSecond() &&
                 this->getQuality() == rhs.getQuality() );
    }

    bool operator != (const DFeaturePair& rhs) const
    {
        return  !(*this == rhs);
    }

    /**@name Accessors
     */
    //@{
    /// Non-mutable access to the first feature
    const FeatureType& getFirst() const
    {
        return this->first;
    }
    /// Mutable access to the first feature
    FeatureType& getFirst()
    {
        return this->first;
    }
    /// Non-mutable access to the first feature
    void setFirst(const FeatureType& frt)
    {
        this->first = frt;
    }

    /// Non-mutable access to the second feature
    const FeatureType& getSecond() const
    {
        return this->second;
    }
    /// Mutable access to the second feature
    FeatureType& getSecond()
    {
        return this->second;
    }
    /// Non-mutable access to the second feature
    void setSecond(const FeatureType& sec)
    {
        this->second = sec;
    }

    /// Non-mutable access to the quality of the pair
    const QualityType& getQuality() const
    {
        return quality_;
    }
    /// Mutable access to the quality of the pair
    QualityType& getQuality()
    {
        return quality_;
    }
    /// Mutable access to the quality of the pair
    void setQuality(const QualityType& ql)
    {
        quality_ = ql;
    }
    //@}

protected:

    /// quality of the pair (not individual features)
    QualityType quality_;

}
; // end of class DFeaturePair

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DFEATUREPAIR_H
