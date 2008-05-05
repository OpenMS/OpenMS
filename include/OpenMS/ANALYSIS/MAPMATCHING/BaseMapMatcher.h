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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEMAPMATCHER_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

#include <utility>

namespace OpenMS
{

  /**
    @brief The base class of the map matching algorithm.
    
    This class defines the basic interface for all mapmatching
    algorithms. It expects a list of element pairs together
    with a quality value for each pair and the coordinates
    of a grid covering the map.
    
    It then estimates the parameters of a transformation
    describing the shift in retention time and m/z between
    the two maps.  
  */
  template <typename ElementT = Feature >
  class BaseMapMatcher
  {
  public:
    /// Element type
    typedef ElementT ElementType;
    /// The element pairs are computed by the element matching class
    typedef std::vector< ElementPair < ElementType > > ElementPairVector;
    /// Quality type
    typedef DoubleReal QualityType;

    /// Constructor
    BaseMapMatcher() : min_quality_(-1)
    {}

    /// Equality operator
    bool operator == (const BaseMapMatcher& rhs)
    {
      return (grid_          == rhs.grid_ &&
              element_pairs_ == rhs.element_pairs_ &&
              min_quality_   == rhs.min_quality_);
    }

    /// Destructor
    virtual ~BaseMapMatcher()
    {}

    /// Set grid
    void setGrid(const LinearMapping& g)
    {
      grid_ = g;
    }
    /// Get grid
    LinearMapping& getGrid()
    {
      return grid_;
    }
    /// Get grid (non-mutable)
    const LinearMapping& getGrid() const
    {
      return grid_;
    }

    /// Set element pair list
    void setElementPairs(const ElementPairVector& plist)
    {
      element_pairs_ = plist;
    }
    /// Get element pair list
    ElementPairVector& getElementPairs()
    {
      return element_pairs_;
    }
    /// Get element pair list (non-mutable)
    const ElementPairVector& getElementPairs() const
    {
      return element_pairs_;
    }

    /// Set quality
    void setMinQuality(QualityType qu)
    {
      min_quality_ = qu;
    }
    /// Get quality
    QualityType& getMinQuality()
    {
      return min_quality_;
    }
    /// Get quality
    QualityType getMinQuality() const
    {
      return min_quality_;
    }

    /// Estimates the transformation for each grid cell
    virtual void estimateTransform() = 0;

  protected:
    /// Vector of DRange instances defining a grid over the map
    LinearMapping grid_;

    /// Vector of pairs of elements that have been identified by the element matcher
    ElementPairVector element_pairs_;

    /// Minimum quality that we accept for element pairs, defined in param class
    QualityType min_quality_;

  }
  ; // end of class BaseMapMatcher

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DBASEMAPMATCHER_H
