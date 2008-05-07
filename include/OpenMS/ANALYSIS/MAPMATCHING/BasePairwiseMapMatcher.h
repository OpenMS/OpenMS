// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                  OpenMS Mass Spectrometry Framework
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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
		@brief The base class of all pairwise point matching algorithms.
		
		This class defines the basic interface for all point matching
		algorithms.  
		It works on two point maps and computes a vector of corresponding points
		in both maps (given by a point pairs vector). 
		A point can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
		(wheras DFeature is the default element type).
		
		The point pairs created by the algorithm solve a
		(bipartite) matching problem between two point maps.
		Therefore first a transformation is estimated, that maps the one map 
		(the so called scene map) onto the other map (the so called model map).
		Given the transformation correspoinding elements in the two maps are determined.
		
		@note If a piecewise transformation is assumed, the user can define a grid by setting 
		the number of buckets in the RT as well as the MZ dimension. 
		Call initGridTransformation() before run()!
 */
  template < typename MapT = FeatureMap< > >
  class BasePairwiseMapMatcher 
  	: public FactoryProduct
  {
  public:
    /// Container for input elements
    typedef MapT PointMapType;

    /// Type of elements considered here
    typedef typename PointMapType::value_type ElementType;

    /// Type of element pairs
    typedef ElementPair < ElementType > ElementPairType;

    /// Container for generated element pairs
    typedef std::vector < ElementPairType > ElementPairVectorType;

    /// Position
    typedef DPosition < 2 > PositionType;

    ///
    typedef DBoundingBox< 2>  PositionBoundingBoxType;

    /// Coordinate
    typedef DoubleReal CoordinateType;


    /// Constructor
    BasePairwiseMapMatcher()
        : FactoryProduct("BasePairWiseMapMatcher")
    {
      element_map_[0] = 0;
      element_map_[1] = 0;
      
			// no need to call defaultsToParam_() as it is called in the non-abstract children 
    }

   
    /// Destructor
    virtual ~BasePairwiseMapMatcher()
  	{
  	}
		
    /// Set element map
    void setElementMap(UInt const index, const PointMapType& element_map)
    {
      element_map_[index] = &element_map;
    }

    /// Get element map
    const PointMapType& getElementMap(UInt index) const
    {
      return *element_map_[index];
    }

    /// Get element pair list
    const ElementPairVectorType& getElementPairs() const
    {
      return all_element_pairs_;
    }

    /// Get grid
    const LinearMapping& getGrid() const
    {
      return grid_;
    }

    /// Register all derived classes here
    static void registerChildren();

    /// Determine corresponding elements (element pairs)
    virtual void run() = 0;

  protected:
  	
    /// Two maps of elements to be matched
    PointMapType const * element_map_[2];

    /// Each element of the vector corresponds to all element pairs of one gridcell
    ElementPairVectorType all_element_pairs_;

    /// The estimated transformation between the two element maps
    LinearMapping grid_;


  } ; // BasePairwiseMapMatcher



} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H
