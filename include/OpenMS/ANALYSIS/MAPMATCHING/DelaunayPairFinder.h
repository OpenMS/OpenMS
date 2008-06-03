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
// $Maintainer: Clemens Groepl, Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

namespace OpenMS
{
  /**
	  @brief This class implements an element pair finding algorithm.
				
	  This class implements a point pair finding algorithm.
	  It offers a method to determine element pairs in two element maps,
	  given two point maps and a transformation defined for the second element map (if no
	  transformation is given, the pairs are found in the two original maps). 
	  The pair finder also offers a method to compute consensus elements given 
	  two element maps. This algorithm is similar to the pair finding method as mentioned above,
	  but it implies that the scene map is already dewarped.
				     
	  To speed up the search for element pairs an consensus elements, the DelaunayPairFinder
	  uses the CGAL delaunay triangulation for the nearest neighbour search.
				
	  The template parameter is the type of the consensus map.
				
	  @note The RT and the MZ dimension are not equivalent, because two elements that differ in RT by 1s (or minute) are 
	  more similar than two points that differ in MZ by 1Th. To be able to use the euclidean distance in the nearest neighbour search, 
	  we have to transform the elements MZ position m into a new MZ position m'= m * internal_mz_scaling;
	  E.g. given internal_mz_scaling=10 results in 1s difference in RT being similar to 0.1Th difference in MZ.
				 
	  @ref DelaunayPairFinder_Parameters are explained on a separate page.  
	
	  @todo work out all TODOs in the code, add offsets (Clemens)
	
	  @ingroup FeatureGrouping
  */
  class DelaunayPairFinder 
  	: public BaseGroupFinder
  {
   public:
		
		struct Point;
		class GeometricTraits;
		struct PointArray2;
		
    typedef BaseGroupFinder Base;
    
		enum { MODEL_=0, SCENE_=1 };		
		enum { RT = RawDataPoint2D::RT, MZ = RawDataPoint2D::MZ };

    /// Constructor
    DelaunayPairFinder();

    /// Destructor
    virtual ~DelaunayPairFinder()
    {
    }

    /// Returns an instance of this class
    static BaseGroupFinder* create()
    {
      return new DelaunayPairFinder();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "delaunay";
    }

    // documented in base class
    void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map);

   protected:

		inline DoubleReal squared_distance_( DoubleReal x1, DoubleReal y1, DoubleReal x2, DoubleReal y2 ) const
		{
			return pow(x1-x2,2) + pow(y1-y2,2); // TODO: check if pow(x,2) is really faster than x*x  (as claimed by AnHi)
		}
		
		virtual void updateMembers_();

		/// Maximal distance of a matched pair, in both dimension RT and MZ
    DoubleReal max_pair_distance_[2];

    /// Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant.
    DoubleReal internal_mz_scaling_;

		/// Upper bound for squared_distance_()
		DoubleReal max_squared_distance_;
		
		/// The distance of the second nearest neighbor must be this factor larger
    DoubleReal second_nearest_gap_;

  }; // class

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
