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
// $Id: SimpleFeaFiTraits.C,v 1.39 2006/05/30 15:46:44 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleFeaFiTraits.h>

namespace OpenMS
{

	SimpleFeaFiTraits::SimpleFeaFiTraits() : BaseFeaFiTraits()
	{
		name_ = SimpleFeaFiTraits::getName();
		defaults_.setValue("min_intensity",0.0f);

		param_ = defaults_;
	}

	SimpleFeaFiTraits::~SimpleFeaFiTraits(){}

	const SimpleFeaFiTraits::Flag& SimpleFeaFiTraits::getPeakFlag(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

 	SimpleFeaFiTraits::Flag& SimpleFeaFiTraits::getPeakFlag(const UnsignedInt index)
		throw (Exception::IndexOverflow)
	{
		if (index>=flags_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, flags_.size());
		return flags_[index];
	}

	const SimpleFeaFiTraits::FlagRefVector& SimpleFeaFiTraits::getFlags(const IndexSet& set)
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

	const SimpleFeaFiTraits::FlagVector& SimpleFeaFiTraits::getAllFlags() const
	{
		return flags_;
	}

	SimpleFeaFiTraits::FlagVector& SimpleFeaFiTraits::getAllFlags()
	{
		return flags_;
	}

	const SimpleFeaFiTraits::PeakType& SimpleFeaFiTraits::getPeak(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		return peaks_[index];
	}

  const SimpleFeaFiTraits::PeakRefVector& SimpleFeaFiTraits::getPeaks(const IndexSet& set)
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
  
	const SimpleFeaFiTraits::PeakVector& SimpleFeaFiTraits::getAllPeaks()
	{
		return peaks_;
	}

	const UnsignedInt SimpleFeaFiTraits::getNumberOfPeaks()
	{
		return peaks_.size();
	}

	const UnsignedInt SimpleFeaFiTraits::getPeakScanNr(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, "SimpleFeaFiTraits::getScanNr()", index, peaks_.size());
		CoordinateType current_rt = getPeakRt(index);
		return scan_index_.getRank(current_rt);
	}

   
	const SimpleFeaFiTraits::IntensityType& SimpleFeaFiTraits::getPeakIntensity(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());

		return peaks_[index].getIntensity();
	}
   
	const SimpleFeaFiTraits::CoordinateType& SimpleFeaFiTraits::getPeakMz(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());

		return peaks_[index].getPosition()[MZ];
	}

	const SimpleFeaFiTraits::CoordinateType& SimpleFeaFiTraits::getPeakRt(const UnsignedInt index) const
		throw (Exception::IndexOverflow)
	{
		if (index>=peaks_.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peaks_.size());
		
		return peaks_[index].getPosition()[RT];
	}
	
	UnsignedInt SimpleFeaFiTraits::getNextMz(UnsignedInt index) const
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
	
	UnsignedInt SimpleFeaFiTraits::getPrevMz(UnsignedInt index) const
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

	UnsignedInt SimpleFeaFiTraits::getNextRt(UnsignedInt index) const
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

	UnsignedInt SimpleFeaFiTraits::getPrevRt(UnsignedInt index) const
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
	
	const SimpleFeaFiTraits::FeatureVector& SimpleFeaFiTraits::run()
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

			} // end of while(true)
		} 
		catch(NoSuccessor ex) { }
		
		if (debug_ > 0) 
			*debug_stream_ << instance_ << " " << features_.size() << " features were found. " << std::endl;
			
		return features_;
	}

	void SimpleFeaFiTraits::addSinglePeak(const DRawDataPoint<2>& peak)
	{
			min_intensity_ = param_.getValue("min_intensity");
			if (peak.getIntensity() > min_intensity_)
			{
				peaks_.push_back(peak);
				flags_.push_back(BaseFeaFiTraits::UNUSED);
			}
	}
	
	void SimpleFeaFiTraits::setData(MSExperiment<DPeak<1> >& exp)
	{
		exp.get2DData(peaks_);	
		
		for (unsigned int i=0; i<peaks_.size(); i++)
			flags_.push_back(BaseFeaFiTraits::UNUSED);
				
		sortData_();
	}
		
	void SimpleFeaFiTraits::sortData_() 
	{
		
		std::sort(peaks_.begin(),
		          peaks_.end(),
		          LexicographicComparator<RTless,MZless>());
		
		scan_index_.init ( peaks_.begin(), peaks_.end() );
				
	}
	
} // end of namespace OpenMS	


