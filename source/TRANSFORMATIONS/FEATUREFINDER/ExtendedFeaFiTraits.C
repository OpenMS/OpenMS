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
// $Id: ExtendedFeaFiTraits.C,v 1.36 2006/05/30 15:46:44 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedFeaFiTraits.h>

namespace OpenMS
{

	ExtendedFeaFiTraits::ExtendedFeaFiTraits() : BaseFeaFiTraits()
	{
		name_ = ExtendedFeaFiTraits::getName();
		defaults_.setValue("min_intensity",0.0f);

		param_ = defaults_;
	}

	ExtendedFeaFiTraits::~ExtendedFeaFiTraits(){}

	const ExtendedFeaFiTraits::Flag& ExtendedFeaFiTraits::getPeakFlag(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

 	ExtendedFeaFiTraits::Flag& ExtendedFeaFiTraits::getPeakFlag(const UnsignedInt index)
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

	const ExtendedFeaFiTraits::FlagRefVector& ExtendedFeaFiTraits::getFlags(const IndexSet& set)
		throw (Exception::IndexOverflow)
	{
		IndexSet::ConstIterator it = --set.end();
		if ( *it >= flags_.size()) // last Index too big
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, *it, flags_.size());

		selected_flags_.clear();
		for (it=set.begin(); it!=set.end(); it++)
			selected_flags_.push_back( &flags_[*it] );

		return selected_flags_;
	}
  
	const ExtendedFeaFiTraits::FlagVector& ExtendedFeaFiTraits::getAllFlags() const
	{
		return flags_;
	}
  
	ExtendedFeaFiTraits::FlagVector& ExtendedFeaFiTraits::getAllFlags()
	{
		return flags_;
	}

	const ExtendedFeaFiTraits::PeakType& ExtendedFeaFiTraits::getPeak(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		return peaks_[index];
	}

  const ExtendedFeaFiTraits::PeakRefVector& ExtendedFeaFiTraits::getPeaks(const IndexSet& set)
		throw (Exception::IndexOverflow)
	{
		IndexSet::ConstIterator it = --set.end();
		if ( *it >= peaks_.size()) // last Index too big
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, *it, peaks_.size());

		selected_peaks_.clear();
		for (it=set.begin(); it!=set.end(); it++)
			selected_peaks_.push_back( &peaks_[*it] );

		return selected_peaks_; 
	}
  
	const ExtendedFeaFiTraits::PeakVector& ExtendedFeaFiTraits::getAllPeaks()
	{
		return peaks_;
	}
	
	const UnsignedInt ExtendedFeaFiTraits::getNumberOfPeaks()
	{
		return peaks_.size();
	}

   
	const ExtendedFeaFiTraits::IntensityType& ExtendedFeaFiTraits::getPeakIntensity(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());

		return peaks_[index].getIntensity();
	}
   
	const ExtendedFeaFiTraits::CoordinateType& ExtendedFeaFiTraits::getPeakMz(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peaks_[index].getPosition()[MZ];
	}

	const ExtendedFeaFiTraits::CoordinateType& ExtendedFeaFiTraits::getPeakRt(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peaks_[index].getPosition()[RT];
	}
	
	const UnsignedInt ExtendedFeaFiTraits::getPeakScanNr(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "SimpleFeaFiTraits::getScanNr()", index, peaks_.size());
		CoordinateType current_rt = getPeakRt(index);
		return scan_index_.getRank(current_rt);
	}
			
	UnsignedInt ExtendedFeaFiTraits::getNextMz(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		if (index == (peaks_.size()-1) ) throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
				
		// check whether we walked out of the current scan i.e. the retention
		// time has changed
		if (getPeakRt(index) != getPeakRt(index+1)) throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
		
		// since we sorted by rt and then by m/z, the peak with the same rt but 
		// the larger m/z is simply one step further in the peak vector			
		return ++index;
	}
	
	UnsignedInt ExtendedFeaFiTraits::getPrevMz(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		// if we are at the beginning of the peak vector, there will be no previous peak ;-)
		if (index == 0) throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
				
		// check whether we walked out of the current scan i.e. the retention
		// time has changed (same problem as above in nextMz() )
		if (getPeakRt(index) != getPeakRt(index-1)) throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
				
		// same as above
		return --index;
	}

	UnsignedInt ExtendedFeaFiTraits::getNextRt(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		const PeakType peak       = getPeak(index);
		PeakVector::iterator iter;
		try
		{
			iter = scan_index_.getNextRt(peak);
		} catch (Exception::Base ex)
		{
			throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
		}
		UnsignedInt peak_index = (iter - peaks_.begin());
		
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
					
		return peak_index;
	}

	UnsignedInt ExtendedFeaFiTraits::getPrevRt(UnsignedInt index) const
		throw (Exception::IndexOverflow, NoSuccessor)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		const PeakType peak       = getPeak(index);
		PeakVector::iterator iter;
		try {
			iter = scan_index_.getPrevRt(peak);
		} catch (Exception::Base ex)
		{
			throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, index);
		}
		UnsignedInt peak_index = (iter - peaks_.begin());
				
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peak_index;
	}

	const ExtendedFeaFiTraits::FeatureVector& ExtendedFeaFiTraits::run()
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
			*debug_stream_ << Date::now() << " " << instance_ << " No data provided! Aborting..." << std::endl;
			return features_;
		}

		// gather information for fitting summary
		std::map<String,int> exception;		//count exceptions
		int no_exceptions = 0;
		std::map<String,int> mz_model;		//count used mz models
		std::map<float,int> mz_stdev;			//count used mz standard deviations
		std::vector<int> charge(5);				//count used charges
		double corr_mean=0.0, corr_max=0.0, corr_min=1.0; //boxplot for correlation
		std::vector<int> corrs(10);	//count correlations in an interval e.g. corr[4] -> [0.4,0.5]

		try {
			while (true)
			{
				UnsignedInt seed = seeders_[0]->nextSeed();
				IndexSet peaks = extenders_[0]->extend(seed);
				try {
					features_.push_back(fitters_[0]->fit(peaks));
					
					#ifdef DEBUG_FEATUREFINDER
					writeGnuPlotFile(peaks,false,nr_feat++);
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
					if (f.getCharge()!=0)	++mz_stdev[p.getValue("MZ:isotope:stdev")];
					
				}catch( BaseModelFitter::UnableToFit ex)
				{
						for (IndexSet::ConstIterator it=peaks.begin(); it!=peaks.end(); ++it) {
							getPeakFlag(*it) = BaseFeaFiTraits::UNUSED;
						}
					*debug_stream_ << Date::now() << " " << instance_ << " " << ex.what() << std::endl;
					++no_exceptions;
					++exception[ex.getName()];
				}

			} // end of while(true)
		}
		catch(NoSuccessor ex) { }

		// Print summary:
		Size size = features_.size();
		
		if (debug_ > 0) 
			*debug_stream_ << Date::now() << " " << instance_ << " " << size << " features were found. " << std::endl;

		*debug_stream_ << Date::now() << " " << instance_ << " FeatureFinder summary:\n"
							<< "Correlation:\n\tminimum: " << corr_min << "\n\tmean: " << corr_mean/size
							<< "\n\tmaximum: " << corr_max << std::endl;
		for (int i=0; i<10; ++i) if (corrs[i]!=0)
		{
			*debug_stream_ << "\t[" << i/10.0 << "," << (i+1)/10.0 << "]: " << corrs[i]*100/size
							<< "% (" << corrs[i] << ")\n";
		}

		*debug_stream_ << "Exceptions:\n";
		for (std::map<String,int>::const_iterator it=exception.begin(); it!=exception.end(); ++it)
		{
			*debug_stream_ << "\t" << it->first << ": " << it->second*100/no_exceptions
								<< "% (" << it->second << ")\n";
		}

		*debug_stream_ << "Chosen mz models:\n";
		for (std::map<String,int>::const_iterator it=mz_model.begin(); it!=mz_model.end(); ++it)
		{
			*debug_stream_ << "\t" << it->first << ": " << it->second*100/size
								<< "% (" << it->second << ")\n";
		}

		*debug_stream_ << "Chosen mz stdevs:\n";
		for (std::map<float,int>::const_iterator it=mz_stdev.begin(); it!=mz_stdev.end(); ++it)
		{
			*debug_stream_ << "\t" << it->first << ": " << it->second*100/(size-charge[0])
								<< "% (" << it->second << ")\n";
		}

		*debug_stream_ << "Charges:\n";
		for (int i=1; i<5; ++i) if (charge[i]!=0)
		{
			*debug_stream_ << "\t+" << i << ": " << charge[i]*100/(size-charge[0])
							<< "% (" << charge[i] << ")\n";
		}

		#ifdef DEBUG_FEATUREFINDER
		IndexSet empty;
		writeGnuPlotFile(empty,true,nr_feat);
		#endif
		
		return features_;
	}

	void ExtendedFeaFiTraits::addSinglePeak(const DRawDataPoint<2>& peak) 
	{
		min_intensity_ = param_.getValue("min_intensity");

		if (peak.getIntensity() > min_intensity_) 
		{
			peaks_.push_back(peak);
			flags_.push_back(BaseFeaFiTraits::UNUSED);
		}
	}

	void ExtendedFeaFiTraits::setData(MSExperiment<DPeak<1> >& exp)
	{
		double it_u = std::numeric_limits<double>::max();
		double it_l = param_.getValue("min_intensity");
		
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
			flags_.push_back(BaseFeaFiTraits::UNUSED);
				
		sortData_();
	}
	
	void ExtendedFeaFiTraits::sortData_() 
	{
		
		if (peaks_.size() == 0)
		{
			*debug_stream_ << Date::now() << " " << instance_ << " FeatureFinder was initialized with empty peak array. Aborting... " << std::endl;
			return;
		}
		
		std::sort(peaks_.begin(),
		          peaks_.end(),
		          LexicographicComparator<RTless,MZless>());

		scan_index_.init ( peaks_.begin(), peaks_.end() );
			
	}	// end sortData()
		
	void ExtendedFeaFiTraits::writeGnuPlotFile(IndexSet peaks, bool last,int nr_feat)
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
	
} // end of namespace OpenMS	


