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
//            and the derived class is registered in the base class.
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

namespace OpenMS
{
  /** 
  	@brief A simple FeatureFinderAlgorithm implementation 
   
    @ingroup FeatureFinder
  */
	template<class PeakType, class FeatureType>
  class FeatureFinderAlgorithmSimple
		: public FeatureFinderAlgorithm<PeakType, FeatureType>,
			public FeatureFinderDefs
  {

	  public:	  	
	    /// default constructor 
	    FeatureFinderAlgorithmSimple()
				: FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
				this->subsections_.push_back("seeder");
				this->subsections_.push_back("extender");
				this->subsections_.push_back("fitter");
				
				this->check_defaults_ =  false;
			}

			virtual Param getDefaultParameters() const
			{
				Param tmp;
				
				SimpleSeeder<PeakType,FeatureType> seeder(this->map_, this->features_, this->ff_);
	  		tmp.insert("seeder:", seeder.getParameters());
	  		tmp.setDescription("seeder", "Settings for the seeder (Determines potential feature regions)");
	  		
	  		SimpleExtender<PeakType,FeatureType> extender(this->map_, this->features_, this->ff_);
	  		tmp.insert("extender:", extender.getParameters());
	  		tmp.setDescription("extender", "Settings for the extender (Collects all peaks belonging to a feature)");
	  		
	  		SimpleModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
	  		tmp.insert("fitter:", fitter.getParameters());
	  		tmp.setDescription("fitter", "Settings for the modefitter (Fits a model to the data determinging the probapility that they represent a feature.)");
				
				return tmp;
			}


			/// Main method for actual FeatureFinder
			virtual void run()
			{
				UInt seed_nr=1;

	  		SimpleSeeder<PeakType,FeatureType> seeder(this->map_, this->features_, this->ff_);
	  		seeder.setParameters(this->getParameters().copy("seeder:",true));
	  		SimpleExtender<PeakType,FeatureType> extender(this->map_, this->features_, this->ff_);
	  		extender.setParameters(this->getParameters().copy("extender:",true));
	  		SimpleModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
	  		fitter.setParameters(this->getParameters().copy("fitter:",true));

		    try
		  	{
		      while (true)
		      {
		      	
		      	std::cout << "===============================" << std::endl;
						std::cout << "Seed # " << seed_nr++ << std::endl;
						IDX seed = seeder.nextSeed();
		        std::cout << "Extending" << std::endl;
		        	
		        ChargedIndexSet set;
		        set.insert(seed);
		        
		        ChargedIndexSet peaks = extender.extend(set);
		        std::cout << "Fitting" << std::endl;
		        
//		        try
//		        {
//		          this->features_->push_back(fitter.fit(peaks));
//		        }
//		        catch( UnableToFit ex)
//		        {
//		        	std::cout << "Failed: " << ex.what() << std::endl;
//		          // set unused flag for all data points
//		          for (IndexSet::const_iterator it=peaks.begin(); it!=peaks.end(); ++it)
//		          {
//		          	this->ff_->getPeakFlag(*it) = UNUSED;
//		          }
//		        }
		
		      } // end of while(true)
		  	}
		    catch(NoSuccessor ex)
		    {
		    }
			}
			
    	static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
    	{
      	return new FeatureFinderAlgorithmSimple();
    	}

    	static const String getProductName()
    	{
      	return "simple";
    	}
		private:
			/// Not implemented
			FeatureFinderAlgorithmSimple& operator=(const FeatureFinderAlgorithmSimple&);
			/// Not implemented
			FeatureFinderAlgorithmSimple(const FeatureFinderAlgorithmSimple&);
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
