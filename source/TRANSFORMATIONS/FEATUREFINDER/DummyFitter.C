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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummyFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>

#include <iostream>
#include <fstream>
#include <numeric>


using namespace std;

namespace OpenMS
{
	using namespace Internal;

	DummyFitter::DummyFitter()
	: BaseModelFitter(),
		counter_(0)
	{
		setName(getProductName());
		
		defaults_.setValue("min_num_peaks:final",5);
		defaults_.setValue("min_num_peaks:extended",10);
		defaults_.setValue("use_fwhm_intensity",0);
		
		defaultsToParam_();
	}

	DummyFitter::~DummyFitter()  { }


  DummyFitter::DummyFitter(const DummyFitter& rhs)
    : BaseModelFitter(rhs)
  {
  }
  
  DummyFitter& DummyFitter::operator= (const DummyFitter& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseModelFitter::operator=(rhs);
    
    return *this;
  }

  Feature DummyFitter::fit(const ChargedIndexSet& set) throw (UnableToFit)
	{		
		// not enough peaks to fit
		if (set.size() < static_cast<UInt>(param_.getValue("min_num_peaks:extended")))
		{
			for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
			{
				traits_->getPeakFlag(*it) = FeaFiTraits::UNUSED;
			}
			
			String mess = String("Skipping feature, IndexSet size too small: ") + set.size();
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-IndexSet", mess.c_str());
		}
				
		// Build Feature
		// The feature coordinate in rt dimension is given
		// by the centroid of the rt model whereas the coordinate
		// in mz dimension is equal to the monoisotopic peak.
		Feature f;
		f.setOverallQuality(1.0);
		f.setCharge(0);		
		
		// set feature coordinates and intensity
		IntensityType intensity_sum  = 0.0;		
		IDX max_intensity_index;
		IntensityType max_intensity = 0.0;
		
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
		{
			intensity_sum += traits_->getPeakIntensity(*it);
			traits_->getPeakFlag(*it) = FeaFiTraits::INSIDE_FEATURE;
			if (traits_->getPeakIntensity(*it) > max_intensity)
			{
				max_intensity_index = *it;
				max_intensity = traits_->getPeakIntensity(*it);
			}		
		}		
		UInt use_fwhm_intensity = param_.getValue("use_fwhm_intensity");
		
		if (use_fwhm_intensity == 0)
		{
			f.setIntensity(intensity_sum);
		}
		else
		{
			f.setIntensity(max_intensity);
// 			IntensityType intensity_avg = intensity_sum /= set.size();
// 			IntensityType intensity_std = 0.0;
// 		
// 			// compute standard deviation of intensities
// 			for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
// 			{
// 				intensity_std += ( (traits_->getPeakIntensity(*it) - intensity_avg) * (traits_->getPeakIntensity(*it) - intensity_avg) );
// 			}		
// 			intensity_std /= set.size();
// 			intensity_std = sqrt(intensity_std);			
// 			
// 			f.setIntensity(2.345 * intensity_std);
		}
		
		f.setRT(traits_->getPeakRt(max_intensity_index));
		f.setMZ(traits_->getPeakMz(max_intensity_index));
		
		traits_->addConvexHull(set, f);
		
		std::cout << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << " Feature " << counter_ << ": (" << f.getRT();
		std::cout	<< "," << f.getMZ() << ") Qual.:" << f.getOverallQuality() << "\n";
		
		// save meta data in feature for TOPPView
		stringstream s;
		s <<  "Feature #" << counter_ << ", +" << f.getCharge() << ", " << set.size() << "->" << set.size() 
			<< ", Corr: (" << f.getOverallQuality() << "," << f.getQuality(1) << "," << f.getQuality(0) << ")";
		f.setMetaValue(3,s.str());
		
		#ifdef DEBUG_FEATUREFINDER
		// write debug output
		CoordinateType rt = f.getRT();
		CoordinateType mz = f.getMZ();
				
		// wrote peaks remaining after model fit
		String fname = String("feature") + counter_ + "_" + rt + "_" + mz;
		ofstream file2(fname.c_str()); 
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
		{
			FeaFiTraits::PositionType2D p = traits_->getPeakPos(*it);
			file2 << p[1] << " " << p[0] << " " << traits_->getPeakIntensity(*it) << "\n";						
		}
		file2.close();
		#endif
		++counter_;

	  return f;
	}

}



