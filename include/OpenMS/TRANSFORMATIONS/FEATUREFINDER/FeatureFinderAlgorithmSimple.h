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
// $Maintainer: Marc Sturm, Marcel Grunert $
// --------------------------------------------------------------------------

// Attention: This include has to be before the include guard.  Otherwise the
// circular dependencies of the two files cannot be resolved.  The problem is
// that both classes are template classes and the derived class is registered
// in the base class.

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
						// add subsections for DefaultParamHandler 
						this->subsections_.push_back("seeder");
						this->subsections_.push_back("extender");
						this->subsections_.push_back("fitter");

						this->check_defaults_ =  false;
					}

				// TODO(Clemens,Marc) Is this still necessary? Cf. the subsections_ mechanism of DefaultParamHandler.
				virtual Param getDefaultParameters() const
				{
					Param tmp;

					SimpleSeeder<PeakType,FeatureType> seeder(this->map_, this->features_, this->ff_);
					tmp.insert("seeder:", seeder.getParameters());
					tmp.setSectionDescription("seeder", "Settings for the seeder (Determines potential feature regions)");

					SimpleExtender<PeakType,FeatureType> extender(this->map_, this->features_, this->ff_);
					tmp.insert("extender:", extender.getParameters());
					tmp.setSectionDescription("extender", "Settings for the extender (Collects all peaks belonging to a feature)");

					SimpleModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
					tmp.insert("fitter:", fitter.getParameters());
					tmp.setSectionDescription("fitter", "Settings for the modefitter (Fits a model to the data determinging the probapility that they represent a feature.)");

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

			    //information for fitting summary
			    std::map<String,UInt> exception; //count exceptions
			    UInt no_exceptions = 0;
			    std::map<String,UInt> mz_model; //count used mz models
			    std::map<float,UInt> mz_stdev; //count used mz standard deviations
			    std::vector<UInt> charge(10); //count used charges
			    DoubleReal corr_mean=0.0, corr_max=0.0, corr_min=1.0; 	//boxplot for correlation

					try
					{
						for(;;)
						{

							std::cout << "===============================" << std::endl;
						
							std::cout << "===Seeder (seed # " << ++seed_nr << ")..." << std::endl;
							IDX seed = seeder.nextSeed();

							std::cout << "===Extender..." << std::endl;
							ChargedIndexSet index_set;
							index_set.insert(seed);
							ChargedIndexSet region;
						 	extender.extend(index_set, region);

							std::cout << "===ModelFitter..." << std::endl;
							try
							{
								this->features_->push_back(fitter.fit(region));

			         	// gather information for fitting summary
			          const Feature& f = this->features_->back();
			
			          DoubleReal corr = f.getOverallQuality();
			          corr_mean += corr;
			          if (corr<corr_min) corr_min = corr;
			          if (corr>corr_max) corr_max = corr;
			
			          // count estimated charge states
			          UInt ch = f.getCharge();
			          if (ch>= charge.size())
			          {
			          	charge.resize(ch);
			          }
			        	charge[ch]++;
			
			          const Param& p = f.getModelDescription().getParam();
			          ++mz_model[ p.getValue("MZ") ];
			
			          DataValue dp = p.getValue("MZ:isotope:stdev");
			          if (dp != DataValue::EMPTY)
			          {
			          	++mz_stdev[dp];
			          }
							}
							catch( UnableToFit ex)
							{
								std::cout << "UnableToFit: " << ex.what() << std::endl;
								
								// set unused flag for all data points
								for (IndexSet::const_iterator it=region.begin(); it!=region.end(); ++it)
								{
									this->ff_->getPeakFlag(*it) = UNUSED;
								}
			          ++no_exceptions;
			          ++exception[ex.getName()];
								
#ifdef DEBUG_FEATUREFINDER
  							writeGnuPlotFile_(peaks,false,nr_feat++);
#endif
							}
						} //for
					} // try
					catch(NoSuccessor ex)
					{
					}
			    // Print summary:
			    UInt size = this->features_->size();
			    std::cout << size << " features were found. " << std::endl;
			
			    std::cout << "FeatureFinder summary:\n"
			    << "Correlation:\n\tminimum: " << corr_min << "\n\tmean: " << corr_mean/size
			    << "\n\tmaximum: " << corr_max << std::endl;
			
			    std::cout << "Exceptions:\n";
			    for (std::map<String,UInt>::const_iterator it=exception.begin(); it!=exception.end(); ++it)
			    {
			      std::cout << "\t" << it->first << ": " << it->second*100/no_exceptions << "% (" << it->second << ")\n";
			    }
			
			    std::cout << "Chosen mz models:\n";
			    for (std::map<String,UInt>::const_iterator it=mz_model.begin(); it!=mz_model.end(); ++it)
			    {
			      std::cout << "\t" << it->first << ": " << it->second*100/size << "% (" << it->second << ")\n";
			    }
			
			    std::cout << "Chosen mz stdevs:\n";
			    for (std::map<float,UInt>::const_iterator it=mz_stdev.begin(); it!=mz_stdev.end(); ++it)
			    {
			      std::cout << "\t" << it->first << ": " << it->second*100/(size-charge[0]) << "% (" << it->second << ")\n";
			    }
			
			    std::cout << "Charges:\n";
			    for (UInt i=1; i<charge.size(); ++i)
			    {
			      if (charge[i]!=0)
			      {
			        std::cout << "\t+" << i << ": " << charge[i]*100/(size-charge[0]) << "% (" << charge[i] << ")\n";
			      }
					}
					
#ifdef DEBUG_FEATUREFINDER
			    IndexSet empty;
			    writeGnuPlotFile_(empty,true,nr_feat);
#endif

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

		}; // FeatureFinderAlgorithmSimple 

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
