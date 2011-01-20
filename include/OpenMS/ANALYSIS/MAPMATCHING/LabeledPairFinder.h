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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_LABELEDPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_LABELEDPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#include <boost/math/tr1.hpp>
#include <cmath>

namespace OpenMS
{
	/**
	 	@brief The LabeledPairFinder allows the matching of labeled features (features with a fixed distance).
		
		Finds feature pairs that have a defined distance in RT and m/z in the same map.
	
		@htmlinclude OpenMS_LabeledPairFinder.parameters

		@todo Implement support for labled MRM experiments, Q1 m/z value and charges. (Andreas)
		@todo Implement support for more than one mass delta, e.g. from missed cleavages and so on (Andreas)
		
		@ingroup FeatureGrouping
	*/
	class OPENMS_DLLAPI LabeledPairFinder
		: public BaseGroupFinder
	{
		
		public:
			
			/// Default constructor
			LabeledPairFinder();
	
			/// Destructor
			inline virtual ~LabeledPairFinder()
			{
			}

	    /// Returns an instance of this class
	    static BaseGroupFinder* create()
	    {
	      return new LabeledPairFinder();
	    }
	
	    /// Returns the name of this module
	    static const String getProductName()
	    {
	      return "labeled_pair_finder";
	    }
	
			/**
				@brief Run the algorithm
				
				@note Exactly one @em input map has to be provided.
				@note The @em output map has to have two file descriptions, containing
				the same file name. The file descriptions have to be labeled 'heavy' and 'light'.
				
				@exception Exception::IllegalArgument is thrown if the input data is not valid.
			*/
			virtual void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map);
	
		protected:
			
			/// return the p-value at position x for the bi-Gaussian distribution with mean @p m and standard deviation @p sig1 (left) and @p sig2 (right)
			inline DoubleReal PValue_(DoubleReal x, DoubleReal m, DoubleReal sig1, DoubleReal sig2)
			{
				if (m<x)
				{
					return 1-boost::math::tr1::erf((x-m)/sig2/0.707106781);
				}
				else
				{
					return 1-boost::math::tr1::erf((m-x)/sig1/0.707106781);
				}
			}
		private:
			
			/// Copy constructor not implemented => private
			LabeledPairFinder(const LabeledPairFinder& source);

			/// Assignment operator not implemented => private
	    LabeledPairFinder& operator=(const LabeledPairFinder& source);

	}; // end of class LabeledPairFinder

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_LABELEDPAIRFINDER_H
