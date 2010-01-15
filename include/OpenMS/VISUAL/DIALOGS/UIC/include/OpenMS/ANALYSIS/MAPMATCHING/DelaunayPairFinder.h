// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

namespace OpenMS
{
  /**
	  @brief This class implements an element pair finding algorithm.

	  It offers a method to determine element pairs across two element maps.
	  The corresponding features must be aligned, but may have small position deviations.

	  To speed up the search for element pairs an consensus elements, the %DelaunayPairFinder
	  uses the CGAL Delaunay triangulation for the nearest neighbor search.

	  @htmlinclude OpenMS_DelaunayPairFinder.parameters

	  @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI DelaunayPairFinder
  	: public BaseGroupFinder
  {
   public:
		///Base class
    typedef BaseGroupFinder Base;

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

		/**
			@brief Run the algorithm

			@note Exactly two @em input maps must be provided.

			@exception Exception::IllegalArgument is thrown if the input data is not valid.
		*/
    void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map);

   protected:

		///@name Internal helper classes and enums
		//@{
		struct Point;
		class GeometricTraits;
		struct PointArray2;
		enum { MODEL_=0, SCENE_=1 };
		enum { RT = Peak2D::RT, MZ = Peak2D::MZ };
		//@}

		///Calculates the squared distance between two-dimensional points
		inline DoubleReal squaredDistance_( DoubleReal x1, DoubleReal y1, DoubleReal x2, DoubleReal y2 ) const
		{
			DoubleReal tmpx = x1-x2;
			tmpx *= tmpx;
			DoubleReal tmpy = y1-y2;
			tmpy *= tmpy;
			return tmpx + tmpy;
		}

		//docu in base class
		virtual void updateMembers_();

		/// Maximal distance of a matched pair, in both dimension RT and MZ
    DoubleReal max_pair_distance_[2];

    /// Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant.
    DoubleReal internal_mz_scaling_;

		/// Upper bound for squaredDistance_()
		DoubleReal max_squared_distance_;

		/// The distance of the second nearest neighbor must be this factor larger
    DoubleReal second_nearest_gap_;

  };

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
