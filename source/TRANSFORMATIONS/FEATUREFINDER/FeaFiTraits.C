// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: FeaFiTraits.C,v 1.9 2006/06/07 09:35:25 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/KERNEL/RangeUtils.h>

namespace OpenMS
{

	FeaFiTraits::FeaFiTraits()
		: min_intensity_(0)
	{}

	FeaFiTraits::~FeaFiTraits() {}
	
	FeaFiTraits::FeaFiTraits(const FeaFiTraits& source) 
		: peaks_(source.peaks_),
		  flags_(source.flags_),
		  scan_index_(source.scan_index_),
		  features_(source.features_),
		  min_intensity_(source.min_intensity_)
	{		
	}
	

	FeaFiTraits& FeaFiTraits::operator = (const FeaFiTraits& source)
	{
		if (&source == this) return *this;
		
		peaks_            = source.peaks_;
		flags_              = source.flags_;
		scan_index_    = source.scan_index_;
		features_          = source.features_;
		min_intensity_ = source.min_intensity_;

		return *this;
	}

	const UnsignedInt FeaFiTraits::getPeakScanNr(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "SimpleFeaFiTraits::getScanNr()", index, peaks_.size());
		CoordinateType current_rt = getPeakRt(index);
		return scan_index_.getRank(current_rt);
	}
			
	UnsignedInt FeaFiTraits::getNextMz(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index, peaks_.size());
		if (index == (peaks_.size()-1) ) throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);
				
		// check whether we walked out of the current scan i.e. the retention
		// time has changed
		if (getPeakRt(index) != getPeakRt(index+1)) throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);
		
		// since we sorted by rt and then by m/z, the peak with the same rt but 
		// the larger m/z is simply one step further in the peak vector			
		return ++index;
	}
	
	UnsignedInt FeaFiTraits::getPrevMz(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, peaks_.size());
		
		// if we are at the beginning of the peak vector, there will be no previous peak ;-)
		if (index == 0) throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
				
		// check whether we walked out of the current scan i.e. the retention
		// time has changed (same problem as above in nextMz() )
		if (getPeakRt(index) != getPeakRt(index-1)) throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
				
		// same as above
		return --index;
	}

	UnsignedInt FeaFiTraits::getNextRt(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, peaks_.size());
		
		const PeakType peak       = getPeak(index);
		PeakVector::iterator iter;
		try
		{
			iter = scan_index_.getNextRt(peak);
		} catch (Exception::Base ex)
		{
			throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
		}
		UnsignedInt peak_index = (iter - peaks_.begin());
		
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, peaks_.size());
					
		return peak_index;
	}

	UnsignedInt FeaFiTraits::getPrevRt(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevRt", index, peaks_.size());
		
		const PeakType peak       = getPeak(index);
		PeakVector::iterator iter;
		try {
			iter = scan_index_.getPrevRt(peak);
		} catch (Exception::Base ex)
		{
			throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevRt", index);
		}
		UnsignedInt peak_index = (iter - peaks_.begin());
				
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__,"FeaFiTraits::getPrevRt", index, peaks_.size());
		
		return peak_index;
	}

  const FeaFiTraits::FeatureVector& FeaFiTraits::run( const std::vector<BaseSeeder*>& seeders,
    	const std::vector<BaseExtender*>& extenders,const std::vector<BaseModelFitter*>& fitters)
	{
		// Visualize seeds and extension in TOPPView: 
		// get all Seeds and save corresponding peaks as "features"
		// get convex hull of the extension and use it for the "feature"
		
		// counts the number of features collected so far,
		// is needed for the gnuplot output.
		#ifdef DEBUG_FEATUREFINDER
		int nr_feat = 0;
		#endif
		
		if (peaks_.size() == 0)
		{
			std::cout << " No data provided! Aborting..." << std::endl;
			return features_;
		}
		
		if (seeders.size() == 0 ||
		    extenders.size() == 0 ||
			fitters.size() == 0)
		{
			std::cout << " No modules set. Aborting..." << std::endl;
			return features_;
		}
		
		// gather information for fitting summary
		std::map<String,int> exception;		//count exceptions
		int no_exceptions = 0;
		std::map<String,int> mz_model;			//count used mz models
		std::map<float,int> mz_stdev;			//count used mz standard deviations
		std::vector<int> charge(5);				//count used charges
		double corr_mean=0.0, corr_max=0.0, corr_min=1.0; //boxplot for correlation
		std::vector<int> corrs(10);	//count correlations in an interval e.g. corr[4] -> [0.4,0.5]
		
		StopWatch watch;
		try {
			while (true)
			{
				UnsignedInt seed = seeders[0]->nextSeed();
				
				watch.start();
				IndexSet peaks = extenders[0]->extend(seed);
				watch.stop();
				std::cout << "Time spent for extension: " << watch.getClockTime() << std::endl;
				watch.reset();
				try {
					
					watch.start();
					features_.push_back(fitters[0]->fit(peaks));
					watch.stop();
					std::cout << "Time spent for fitting: " << watch.getClockTime() << std::endl;
					watch.reset();
					
					#ifdef DEBUG_FEATUREFINDER
					writeGnuPlotFile_(peaks,false,nr_feat++);
					#endif
					
					// gather information for fitting summary
					const DFeature<2>& f = features_[features_.size()-1];

					float corr = f.getOverallQuality();
					corr_mean += corr;
					if (corr<corr_min) corr_min = corr;
					if (corr>corr_max) corr_max = corr;
					++corrs[ int(corr*10) ];

					charge[f.getCharge()]++;
					const Param& p = f.getModelDescription().getParam();
					++mz_model[ p.getValue("MZ") ];
					
					DataValue dp = p.getValue("MZ:isotope:stdev");
					if (!dp.isEmpty() && dp.toString() != "")
					{
						++mz_stdev[p.getValue("MZ:isotope:stdev")];
					}
					
				}catch( BaseModelFitter::UnableToFit ex)
				{
					for (IndexSet::ConstIterator it=peaks.begin(); it!=peaks.end(); ++it) 
					{
						getPeakFlag(*it) = FeaFiTraits::UNUSED;
					}
					std::cout << " " << ex.what() << std::endl;
					watch.stop();
					std::cout << "Time spent for fitting: " << watch.getClockTime() << std::endl;
					watch.reset();
					++no_exceptions;
					++exception[ex.getName()];
				}

			} // end of while(true)
		}
		catch(NoSuccessor ex) { }

		// Print summary:
		Size size = features_.size();
		
		std::cout << size << " features were found. " << std::endl;

		std::cout << "FeatureFinder summary:\n"
							<< "Correlation:\n\tminimum: " << corr_min << "\n\tmean: " << corr_mean/size
							<< "\n\tmaximum: " << corr_max << std::endl;
		for (int i=0; i<10; ++i) if (corrs[i]!=0)
		{
			std::cout << "\t[" << i/10.0 << "," << (i+1)/10.0 << "]: " << corrs[i]*100/size
							  << "% (" << corrs[i] << ")\n";
		}

		std::cout << "Exceptions:\n";
		for (std::map<String,int>::const_iterator it=exception.begin(); it!=exception.end(); ++it)
		{
			std::cout << "\t" << it->first << ": " << it->second*100/no_exceptions
								<< "% (" << it->second << ")\n";
		}

		std::cout << "Chosen mz models:\n";
		for (std::map<String,int>::const_iterator it=mz_model.begin(); it!=mz_model.end(); ++it)
		{
			std::cout << "\t" << it->first << ": " << it->second*100/size
								<< "% (" << it->second << ")\n";
		}

		std::cout << "Chosen mz stdevs:\n";
		for (std::map<float,int>::const_iterator it=mz_stdev.begin(); it!=mz_stdev.end(); ++it)
		{
			std::cout << "\t" << it->first << ": " << it->second*100/(size-charge[0])
								<< "% (" << it->second << ")\n";
		}

		std::cout << "Charges:\n";
		for (int i=1; i<5; ++i) if (charge[i]!=0)
		{
			std::cout << "\t+" << i << ": " << charge[i]*100/(size-charge[0])
							<< "% (" << charge[i] << ")\n";
		}

		#ifdef DEBUG_FEATUREFINDER
		IndexSet empty;
		writeGnuPlotFile_(empty,true,nr_feat);
		#endif
		
		return features_;
	} // end of run()

	void FeaFiTraits::setData(MSExperiment<DPeak<1> >& exp)
	{
		double it_u = std::numeric_limits<double>::max();
		double it_l = min_intensity_;
						
		// remove data points with intensity < it_l 
		for (MSExperiment< >::iterator it = exp.begin(); it!= exp.end(); ++it)
		{
			it->getContainer().erase(remove_if(it->begin(), it->end(), IntensityRange<MSExperiment< >::PeakType>(it_l, it_u, true)), it->end());
		}
		
		// resize internal data structures
		flags_.reserve(exp.getSize());
		peaks_.reserve(exp.getSize());

		exp.get2DData(peaks_);	
			
		for (unsigned int i=0; i<peaks_.size(); i++)
			flags_.push_back(FeaFiTraits::UNUSED);
				
		sortData_();
	}
			
	void FeaFiTraits::sortData_() 
	{
		
		if (peaks_.size() == 0)
		{
			std::cout << " FeatureFinder was initialized with empty peak array. Aborting... " << std::endl;
			return;
		}
		
		std::sort(peaks_.begin(),
		             peaks_.end(),
		             LexicographicComparator<RTless,MZless>());

		scan_index_.init ( peaks_.begin(), peaks_.end() );
		
		/// Temporal storage for the calculation of s/n ratios
		std::vector<PeakType> last_scan_;
		
		/// Estimates the signal to noise ratio 
		DSignalToNoiseEstimatorWindowing<2> sn_estimator_;
		sn_estimator_.setWindowSize(1000);
		
		String gp_fname("sn_ratios.txt");
		ofstream outfile( gp_fname.c_str() );
		
		// estimate s/n ratios
		for(PeakVector::const_iterator citer = peaks_.begin();
		     citer != peaks_.end();
			 ++citer)
		{
			if (last_scan_.size() > 0)
			{
				// check whether retention time has changed (i.e. a new scan has begun)
				if (citer->getPosition()[0] != last_scan_.back().getPosition()[0])
				{
					// estimate noise for last scan
					sn_estimator_.init(last_scan_.begin(),last_scan_.end());
			
					for (std::vector<PeakType>::const_iterator cit = last_scan_.begin();
			      	 	   cit != last_scan_.end(); 
				      	   ++cit)
					{
						// save s/n values
						double sn = sn_estimator_.getSignalToNoise(cit);
						if (sn < 0) 
						{
							std::cout << "Negative sn !" << std::endl;
							sn = 0;
						}  
						sn_ratios_.push_back(sn);
						outfile << cit->getPosition()[0] << " " << cit->getPosition()[1] << " "  << sn_estimator_.getSignalToNoise(cit) << std::endl;
					}
					// empty container
					last_scan_.clear();
				}
			}
						
			// store new peak
			last_scan_.push_back(*citer);		
		}
		
		// estimate noise for last scan
		sn_estimator_.init(last_scan_.begin(),last_scan_.end());
			
		for (std::vector<PeakType>::const_iterator cit = last_scan_.begin();
		 	   cit != last_scan_.end(); 
		   	   ++cit)
		{
			// save s/n values
			sn_ratios_.push_back(sn_estimator_.getSignalToNoise(cit));
			outfile << cit->getPosition()[0] << " " << cit->getPosition()[1] << " "  << sn_estimator_.getSignalToNoise(cit) << std::endl;
		}
		// empty container
		last_scan_.clear();
		
		outfile.close();
		std::cout << sn_ratios_.size() << " vs " << peaks_.size() << std::endl;
			
	}	// end sortData()
		
	void FeaFiTraits::writeGnuPlotFile_(IndexSet peaks, bool last,int nr_feat)
	{
		
		// ONLY FOR DEBUGGING PURPOSES:
		// write feature + surrounding region to file
		if (!last)
		{	
				String gp_fname("plot.gp");
				ofstream gpfile( gp_fname.c_str(), std::ios_base::app  );
				  			
  			String file_pl = "feature" + String(nr_feat);
  			String file    = "feature" + String(nr_feat);
  		
  			// write feature to output stream
  			gpfile << "replot \'" << file_pl << "\' w i title \"\" " << std::endl;
  			 			
  			ofstream myfile(file.c_str()); // data file
  			IndexSet::const_iterator citer = peaks.begin();
  			
  			while (citer != peaks.end())
  			{
  				myfile << getPeakRt(*citer) << " " << getPeakMz(*citer) << " " << getPeakIntensity(*citer) << std::endl;
  				citer++;
  			}
  			myfile.close();
  	  	gpfile.close();	
  	  	
		} else {

			String gp_fname("plot.gp");

			ofstream gpfile( gp_fname.c_str() , std::ios_base::app );
			gpfile << "pause -1 \'Please hit enter to continue....\' " << std::endl;
			gpfile.close();
  	} 			  	  		
  	  	
	}

	const FeaFiTraits::ConvexHullType FeaFiTraits::calculateConvexHull(const IndexSet& set)
	{
		ConvexHullType convex_hull;
		const double PRECISION = 0.0001;
		if (set.size()<3)
			return convex_hull;

		// keep track of already in hull included peaks to avoid unnecessary computations of triangle area
		std::map<UnsignedInt, bool> isIncluded;

		CoordinateType min_mz = std::numeric_limits<CoordinateType>::max();
		IndexSet::const_iterator min = set.begin();

		// Find peak with minimal mz to start wrapping
		for (IndexSet::const_iterator it = set.begin(); it!=set.end(); ++it)
		{
			if (getPeakMz(*it) < min_mz)
			{
				min_mz = getPeakMz(*it);
				min = it;
			}
			isIncluded[*it] = false;
		}
		convex_hull.push_back( getPeak(*min).getPosition() );

		// Hull peaks denoting current hull line
		IndexSet::const_iterator hull_peak1 = min;
		IndexSet::const_iterator start = set.begin();
		if (start==min) ++start;  // don't start at "min" because of while-condition
		IndexSet::const_iterator hull_peak2 = start;

		while (hull_peak2!=min)
		{
			bool found_any = false;
			for (IndexSet::const_iterator it = set.begin(); it!=set.end(); ++it)
			{
				// skip if already used
				if (isIncluded[*it] || it==hull_peak1 || it==hull_peak2) continue;

				found_any = true;
				// "it" lies to the right of the line [hull_peak1,hull_peak2]
				double area = triangleArea_(hull_peak1,hull_peak2,it);
				if (area>-PRECISION)
				{
					// area almost 0 -> collinear points
					// -> avoid consecutive peaks with equal mz or rt coordinate
					if (fabs(area)<PRECISION)
					{
						double mz1 = getPeakMz(*hull_peak1);
						double mz2 = getPeakMz(*hull_peak2);
						double mz3 = getPeakMz(*it);
						double rt1 = getPeakRt(*hull_peak1);
						double rt2 = getPeakRt(*hull_peak2);
						double rt3 = getPeakRt(*it);
			      if ( 	( fabs(mz2-mz3)<PRECISION && fabs(rt2-rt1) > fabs(rt3-rt1) )
								||( fabs(rt2-rt3)<PRECISION && fabs(mz2-mz1) > fabs(mz3-mz1) ))
						{
							isIncluded[*it] = true;
							continue;
						}
					}
					hull_peak2 = it;  // "it" becomes new hull peak
				}
			}

			if (!found_any){
				hull_peak2 = min; // no available peaks anymore
				continue;
			}

			if (hull_peak2 == min) continue;  // finish loop
			isIncluded[*hull_peak2] = true;

			// continue wrapping
			hull_peak1 = hull_peak2;
			// hull_peak2 satisfies the contition: all peaks lie to the left of [hull_peak1,hull_peak2]
			convex_hull.push_back( getPeak(*hull_peak2).getPosition() );


			start = set.begin();
			if (start==min) ++start;  // don't start at "min" because of while-condition
			hull_peak2 = start;
		}

		return convex_hull;
	}
	
	
} // end of namespace OpenMS	


