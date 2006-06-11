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
// $Id: FastFeaFiTraits.C,v 1.7 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FastFeaFiTraits.h>

namespace OpenMS
{

	FastFeaFiTraits::FastFeaFiTraits() : BaseFeaFiTraits()
	{
		name_ = FastFeaFiTraits::getName();
	}

	FastFeaFiTraits::~FastFeaFiTraits(){}

	const FastFeaFiTraits::Flag& FastFeaFiTraits::getPeakFlag(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

 	FastFeaFiTraits::Flag& FastFeaFiTraits::getPeakFlag(const UnsignedInt index)
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

	const FastFeaFiTraits::FlagRefVector& FastFeaFiTraits::getFlags(const IndexSet& set)
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

	const FastFeaFiTraits::FlagVector& FastFeaFiTraits::getAllFlags() const
	{
		return flags_;
	}

	FastFeaFiTraits::FlagVector& FastFeaFiTraits::getAllFlags()
	{
		return flags_;
	}

	const FastFeaFiTraits::PeakType& FastFeaFiTraits::getPeak(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		return peaks_[index];
	}

  const FastFeaFiTraits::PeakRefVector& FastFeaFiTraits::getPeaks(const IndexSet& set)
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
  
	const FastFeaFiTraits::PeakVector& FastFeaFiTraits::getAllPeaks()
	{
		return peaks_;
	}
	
	const UnsignedInt FastFeaFiTraits::getNumberOfPeaks()
	{
		return peaks_.size();
	}

	const UnsignedInt FastFeaFiTraits::getPeakScanNr(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "FastFeaFiTraits::getScanNr()", index, peaks_.size());
		CoordinateType current_rt = getPeakRt(index);
		return scan_index_.getRank(current_rt);
	}

   
	const FastFeaFiTraits::IntensityType& FastFeaFiTraits::getPeakIntensity(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());

		return peaks_.at(index).getIntensity();
	}
   
	const FastFeaFiTraits::CoordinateType& FastFeaFiTraits::getPeakMz(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peaks_.at(index).getPosition()[MZ];
	}

	const FastFeaFiTraits::CoordinateType& FastFeaFiTraits::getPeakRt(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peaks_.at(index).getPosition()[RT];
	}

	UnsignedInt FastFeaFiTraits::getNextMz(UnsignedInt index) const
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
	
	UnsignedInt FastFeaFiTraits::getPrevMz(UnsignedInt index) const
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

	UnsignedInt FastFeaFiTraits::getNextRt(UnsignedInt index) const
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
		UnsignedInt peak_index    = (iter - peaks_.begin());
		
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
					
		return peak_index;
	}

	UnsignedInt FastFeaFiTraits::getPrevRt(UnsignedInt index) const
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
		
		UnsignedInt peak_index    = (iter - peaks_.begin());
				
		if (peak_index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peak_index;
	}
	
	const FastFeaFiTraits::FeatureVector& FastFeaFiTraits::run()
	{
		// Visualize seeds and extension in TOPPView: 
		// get all Seeds and save corresponding peaks as "features"
		// get convex hull of the extension and use it for the "feature"
		try {
			while (true) 
			{
				UnsignedInt seed = seeders_[0]->nextSeed();
				IndexSet peaks = extenders_[0]->extend(seed);
				
				try {
					features_.push_back(fitters_[0]->fit(peaks));
				}catch( BaseModelFitter::UnableToFit ex){}

				/*
				// ONLY FOR DEBUGGING Seeder with TOPPView: Use seed as feature center
				const DPeak<2>&  p = getPeak(seed);
				DFeature<2>& f = features_[features_.size()-1];
				f.getIntensity() = p.getIntensity();
				f.getPosition()[RT] = p.getPosition()[RT];
				f.getPosition()[MZ] = p.getPosition()[MZ];

				// ONLY FOR DEBUGGING Extender with TOPPView: Use convex hull of the extended region
				f.getConvexHulls().push_back(calculateConvexHull_(peaks));
				*/
			} // end of while(true)
		} 
		catch(NoSuccessor ex) { }
		
		if (debug_ > 0) 
			*debug_stream_ << instance_ << " " << features_.size() << " features were found. " << std::endl;
			
		return features_;
	}

	void FastFeaFiTraits::addSinglePeak(const DRawDataPoint<2>& peak) 
	{
			peaks_.push_back(peak);
			flags_.push_back(BaseFeaFiTraits::UNUSED);
	}
	
	void FastFeaFiTraits::setData(const MSExperiment<DPeak<1> >& exp)
	{
		exp.get2DData(peaks_);	
		
		for (unsigned int i=0; i<peaks_.size(); i++)
			flags_.push_back(BaseFeaFiTraits::UNUSED);
				
		sortData_();
	}
	
	void FastFeaFiTraits::sortData_() 
	{
		
		std::sort(peaks_.begin(),
		          peaks_.end(),
		          LexicographicComparator<RTless,MZless>());
		
		scan_index_.init ( peaks_.begin(), peaks_.end() );
	}
	
} // end of namespace OpenMS	
