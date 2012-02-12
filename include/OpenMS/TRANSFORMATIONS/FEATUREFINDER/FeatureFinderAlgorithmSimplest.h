// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLEST_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLEST_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>

namespace OpenMS
{
	/**
		@brief FeatureFinderAlgorithm implementation using the Simple* modules.

		@deprecated Deprecated in OpenMS 1.7.

    SimpleSeeder, SimpleExtender, ModelFitter (using BiGaussModel in RT dimension and IsotopeModel (charge does not equal zero) or GaussModel in dimension of mz).

    @htmlinclude OpenMS_FeatureFinderAlgorithmSimplest.parameters

		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithmSimplest :
		public FeatureFinderAlgorithm<PeakType, FeatureType>,
		public FeatureFinderDefs
	{

		public:
			/// default constructor
			FeatureFinderAlgorithmSimplest() :
				FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
				this->defaults_ = getDefaultParameters();
				this->check_defaults_ =  false;
			}

			virtual Param getDefaultParameters() const
			{
				Param tmp;

				SimpleSeeder<PeakType,FeatureType> seeder(this->map_, this->features_, this->ff_);
				tmp.insert("seeder:", seeder.getParameters());
				tmp.setSectionDescription("seeder", "Settings for the seeder (Determines potential feature regions)");

				SimpleExtender<PeakType,FeatureType> extender(this->map_, this->features_, this->ff_);
				tmp.insert("extender:", extender.getParameters());
				tmp.setSectionDescription("extender", "Settings for the extender (Collects all peaks belonging to a feature)");

				ModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
				tmp.insert("fitter:", fitter.getParameters());
				tmp.setSectionDescription("fitter", "Settings for the modefitter (Fits a model to the data determinging the probapility that they represent a feature.)");

				return tmp;
			}

			virtual void run()
			{
#ifdef DEBUG_FEATUREFINDER
				UInt seed_nr=0;
#endif
				SimpleSeeder<PeakType,FeatureType> seeder(this->map_, this->features_, this->ff_);
				seeder.setParameters(this->getParameters().copy("seeder:",true));

				SimpleExtender<PeakType,FeatureType> extender(this->map_, this->features_, this->ff_);
				extender.setParameters(this->getParameters().copy("extender:",true));

				ModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
        Param params;
        params.setDefaults(this->getParameters().copy("fitter:",true));
        params.setValue("fit_algorithm", "simplest");
        fitter.setParameters(params);

				/// Summary of fitting results
				Summary summary;

				try
				{
					for(;;)
					{
#ifdef DEBUG_FEATUREFINDER
						std::cout << "===============================" << std::endl;
						std::cout << "### Seeder (seed # " << ++seed_nr << ")..." << std::endl;
#endif
						IndexPair seed = seeder.nextSeed();

#ifdef DEBUG_FEATUREFINDER
            std::cout << "seed ... " << seed.first << " - " << seed.second << std::endl;
						std::cout << "### Extender..." << std::endl;
#endif
						ChargedIndexSet index_set;
						index_set.insert(seed);
						ChargedIndexSet region;
						extender.extend(index_set, region);

#ifdef DEBUG_FEATUREFINDER
						std::cout << "### ModelFitter..." << std::endl;
#endif
						try
						{
							this->features_->push_back(fitter.fit(region));

							// gather information for fitting summary
							{
								const Feature& f = this->features_->back();

								// quality, correlation
								DoubleReal corr = f.getOverallQuality();
								summary.corr_mean += corr;
								if (corr<summary.corr_min) summary.corr_min = corr;
								if (corr>summary.corr_max) summary.corr_max = corr;

								// charge
								UInt ch = f.getCharge();
								if (ch>= summary.charge.size())
								{
									summary.charge.resize(ch+1);
								}
								summary.charge[ch]++;

								// MZ model type
								const Param& p = f.getModelDescription().getParam();
								++summary.mz_model[ p.getValue("MZ") ];

								// standard deviation of isotopic peaks
								if (p.exists("MZ:isotope:stdev") && p.getValue("MZ:isotope:stdev")!=DataValue::EMPTY)
								{
									++summary.mz_stdev[p.getValue("MZ:isotope:stdev")];
								}
							}
						}
						catch( Exception::UnableToFit ex)
						{
							std::cout << "UnableToFit: " << ex.what() << std::endl;

							// set unused flag for all data points
							for (IndexSet::const_iterator it=region.begin(); it!=region.end(); ++it)
							{
								this->ff_->getPeakFlag(*it) = UNUSED;
							}

							// gather information for fitting summary
							{
								++summary.no_exceptions;
								++summary.exception[ex.getName()];
							}
						}
					} // for
				} // try
				catch(NoSuccessor ex)
				{
				}

				this->ff_->endProgress();
			
				// print fitting summary
				{
					Size size = this->features_->size();
					std::cout << size << " features were found. " << std::endl;

					// compute corr_mean
					summary.corr_mean /= size;

					std::cout << "FeatureFinder summary:\n"
						<< "Correlation:\n\tminimum: " << summary.corr_min << "\n\tmean: " << summary.corr_mean
						<< "\n\tmaximum: " << summary.corr_max << std::endl;

					std::cout << "Exceptions:\n";
					for (std::map<String,UInt>::const_iterator it=summary.exception.begin(); it!=summary.exception.end(); ++it)
					{
						std::cout << "\t" << it->first << ": " << it->second*100/summary.no_exceptions << "% (" << it->second << ")\n";
					}

					std::cout << "Chosen mz models:\n";
					for (std::map<String,UInt>::const_iterator it=summary.mz_model.begin(); it!=summary.mz_model.end(); ++it)
					{
						std::cout << "\t" << it->first << ": " << it->second*100/size << "% (" << it->second << ")\n";
					}

					std::cout << "Chosen mz stdevs:\n";
					for (std::map<float,UInt>::const_iterator it=summary.mz_stdev.begin(); it!=summary.mz_stdev.end(); ++it)
					{
						std::cout << "\t" << it->first << ": " << it->second*100/(size-summary.charge[0]) << "% (" << it->second << ")\n";
					}

					std::cout << "Charges:\n";
					for (Size i=1; i<summary.charge.size(); ++i)
					{
						if (summary.charge[i]!=0)
						{
							std::cout << "\t+" << i << ": " << summary.charge[i]*100/(size-summary.charge[0]) << "% (" << summary.charge[i] << ")\n";
						}
					}
				}
			} // run

			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmSimplest();
			}

			static const String getProductName()
			{
				return "simplest";
			}
		private:
			/// Not implemented
			FeatureFinderAlgorithmSimplest& operator=(const FeatureFinderAlgorithmSimplest&);
			/// Not implemented
			FeatureFinderAlgorithmSimplest(const FeatureFinderAlgorithmSimplest&);

	}; // FeatureFinderAlgorithmSimplest

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLEST_H
