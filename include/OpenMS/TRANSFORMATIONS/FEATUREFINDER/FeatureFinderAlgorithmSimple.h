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
// $Maintainer: Marc Sturm, Marcel Grunert$
// --------------------------------------------------------------------------

// Attention: This include has to be before the include guard.
//            Otherwise the circular dependencies of the two files
//            cannot be resolved.
//            The problem is that both classes are template classes
//            and the dericed class is registered in the base class.
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H


namespace OpenMS
{
  /** 
  	@brief A simple FeatureFinderAlgorithm implementation 
   
    @ingroup FeatureFinder
  */
	template<class PeakType, class FeatureType>
  class FeatureFinderAlgorithmSimple
		: public FeatureFinderAlgorithm<PeakType, FeatureType>
  {

	  public:	  	
	    /// default constructor 
	    FeatureFinderAlgorithmSimple()
				: FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
			}
	
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				std::cout << "IT WORKED! YEAH!!!!" << std::endl;
				for (UInt i=0; i<20; ++i)
				{
					std::cout << this->getPeakIntensity_(std::make_pair(i,i)) << std::endl;
				}
			}
			
    	static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
    	{
      	return new FeatureFinderAlgorithmSimple();
    	}

    	static const String getProductName()
    	{
      	return "FeatureFinderAlgorithmSimple";
    	}
		private:
			// not implemented -> private
			FeatureFinderAlgorithmSimple& operator=(const FeatureFinderAlgorithmSimple&);
			// not implemented -> private
			FeatureFinderAlgorithmSimple(const FeatureFinderAlgorithmSimple&);
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
