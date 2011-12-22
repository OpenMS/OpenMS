// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#define V_SimplePairFinder(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /** @brief This class implements a simple point pair finding algorithm.

  It offers a method to find element pairs across two element maps.

	The similarity value should express our confidence that one element might
	possibly be matched to the other.  Larger quality values are better, the
	maximal similarity is one.  Let \f$\Delta_\textit{RT}\f$ and
	\f$\Delta_\textit{MZ}\f$ be the absolute values of the RT and MZ differences
	in the data.  Then the similarity value is
	\f[
	\frac{1}{
	\big( 1 + \Delta_\textit{RT} \cdot \textit{diff\_intercept\_RT} \big)^\textit{diff\_exponent\_RT}
	\cdot
	\big( 1 + \Delta_\textit{MZ} \cdot \textit{diff\_intercept\_MZ} \big)^\textit{diff\_exponent\_MZ}
	}
	\f]
	
	Choosing diff_exponent: This parameter controls the growth rate of the
	penalty for differences.  It is for example possible to search for pairs
	using the absolute distance in RT (which should not be very susceptible to
	outliers) and the squared distance in MZ (where small difference occur
	frequently, but large differences indicate a mismatch).
	
	Choosing diff_intercept: Since we are taking the reciprocal value ("1/..."),
	we include an offset to avoid division by zero in case \f$\Delta=0\f$.  To
	set this parameter, ask yourself: How much worse is a difference of 1
	compared to no difference?

	The following image illustrates the influence of these parameters:
	\image html SimplePairFinder.png "Influence of the parameters intercept and exponent"
	
	@htmlinclude OpenMS_SimplePairFinder.parameters
	
	@ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI SimplePairFinder
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

    /**@brief Compute the similarity for a pair of elements.
    */
    DoubleReal similarity_ ( ConsensusFeature const & left, ConsensusFeature const & right) const;

  }
  ; // SimplePairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_SimplePairFinder_H

/*

gnuplot history - how the plot was created - please do not delete this receipt

f(x,intercept,exponent)=1/(1+(abs(x)*intercept)**exponent)
set terminal postscript enhanced color
set output "choosingsimplepairfinderparams.ps"
set size ratio .3
plot [-3:3] [0:1] f(x,1,1), f(x,2,1), f(x,1,2), f(x,2,2)

*/
