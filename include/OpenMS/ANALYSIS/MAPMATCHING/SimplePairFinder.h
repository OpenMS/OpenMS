// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#define V_SimplePairFinder(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /**
		@brief This class implements a simple point pair finding algorithm.
		
		It offers a method to find element pairs across two element maps.
		In contrast to the @ref DelaunayPairFinder , a different distance measure is used.
		
		@todo Add formula for distance measure (Clemens)
		
		@ref SimplePairFinder_Parameters are explained on a separate page.

		@ingroup FeatureGrouping
  */
  class SimplePairFinder
  	: public BaseGroupFinder
  {
	 public:
	 	///Base class
    typedef BaseGroupFinder Base;

    /// Constructor
    SimplePairFinder();
    /// Destructor
    virtual ~SimplePairFinder()
		{
		}

    /// returns an instance of this class
    static BaseGroupFinder* create()
    {
      return new SimplePairFinder();
    }
    /// returns the name of this module
    static const String getProductName()
    {
      return "simple";
    }

		/**
			@brief Run the algorithm
			
			@note Exactly two @em input maps must be provided.
			@note All two @em input maps must be provided.
			
			@exception Exception::IllegalArgument is thrown if the input data is not valid.
		*/
		virtual void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map);

	 protected:
	 	
	 	//docu in base class
    virtual void updateMembers_();

    /// A parameter for similarity_().
    DoubleReal diff_exponent_[2];

    /// A parameter for similarity_().
    DoubleReal diff_intercept_[2];

    /// Minimal pair quality
    DoubleReal pair_min_quality_;

    /**@brief Compute the similarity for a pair of elements; larger quality
		values are better.

		The returned value should express our confidence that one element might
		possibly be matched to the other.

		The similarity is computed as follows.
		- For each dimension:
		<ul>
		<li> Take the absolute difference of the coordinates.
		<li> Add #diff_intercept_ to it.
		<li> Raise the sum to power of #diff_exponent_.
		</ul>
		- Multiply these numbers for both dimensions.
		- Take the reciprocal value of the result.
		.

		The parameter #diff_exponent_ controls the asymptotic decay rate for large
		differences.  The parameter #diff_intercept_ is important for small
		differences.

    */
    DoubleReal similarity_ ( ConsensusFeature const & left, ConsensusFeature const & right) const;

  }
  ; // SimplePairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_SimplePairFinder_H
