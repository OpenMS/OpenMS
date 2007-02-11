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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakFitter.h>
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

	PeakFitter::PeakFitter()
	: BaseModelFitter(),
		counter_(0)
	{
		setName(getProductName());
		
		defaults_.setValue("min_num_peaks:final",5);
		defaults_.setValue("min_num_peaks:extended",10);
		defaults_.setValue("use_max_intensity",0);
		
	}

	PeakFitter::~PeakFitter()  { }


  PeakFitter::PeakFitter(const PeakFitter& rhs)
    : BaseModelFitter(rhs)
  {
  }
  
  PeakFitter& PeakFitter::operator= (const PeakFitter& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseModelFitter::operator=(rhs);
    
    return *this;
  }

  DFeature<2> PeakFitter::fit(const IndexSet& set) throw (UnableToFit)
	{		
		// not enough peaks to fit
		if (set.size() < static_cast<Size>(param_.getValue("min_num_peaks:extended")))
		{
			String mess = String("Skipping feature, IndexSet size too small: ") + set.size();
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-IndexSet", mess.c_str());
		}
				
		// Build Feature
		// The feature coordinate in rt dimension is given
		// by the centroid of the rt model whereas the coordinate
		// in mz dimension is equal to the monoisotopic peak.
		DFeature<2> f;
		f.setOverallQuality(1.0);
		f.setCharge(0);		
		
		// set feature coordinates and intensity
		IntensityType intensity_sum  = 0.0;		
		IntensityType max_intensity  = 0.0;
		IDX max_intensity_index;
		
		// intensity of the feature is either the sum of all included data points 
		// or the maximum ion count in the feature region
		// coordinates are given by point with max intensity
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
		{
			intensity_sum += traits_->getPeakIntensity(*it);
			if (traits_->getPeakIntensity(*it) > max_intensity)
			{
				max_intensity          = traits_->getPeakIntensity(*it);
				max_intensity_index = *it;
			} 
		}		
		
		UnsignedInt use_max_intensity = param_.getValue("use_max_intensity");
		
		if (use_max_intensity == 0)
		{
			f.setIntensity(intensity_sum);
		}
		else
		{
			f.setIntensity(max_intensity);
		}
		
		f.getPosition()[RT] = traits_->getPeakRt(max_intensity_index);
		f.getPosition()[MZ] = traits_->getPeakMz(max_intensity_index);
		
		traits_->addConvexHull(set, f);
		
		std::cout << Date::now() << " Feature " << counter_ << ": (" << f.getPosition()[RT];
		std::cout	<< "," << f.getPosition()[MZ] << ") Qual.:" << f.getOverallQuality() << "\n";
		
		// save meta data in feature for TOPPView
		stringstream s;
		s <<  "Feature #" << counter_ << ", +" << f.getCharge() << ", " << set.size() << "->" << set.size() 
			<< ", Corr: (" << f.getOverallQuality() << "," << f.getQuality(RT) << "," << f.getQuality(MZ) << ")";
		f.setMetaValue(3,s.str());
		
		#ifdef DEBUG_FEATUREFINDER
		// write debug output
		CoordinateType rt = f.getPosition()[RT];
		CoordinateType mz = f.getPosition()[MZ];
				
		// wrote peaks remaining after model fit
		String fname = String("feature") + counter_ + "_" + rt + "_" + mz;
		ofstream file2(fname.c_str()); 
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
		{
			FeaFiTraits::PositionType2D p = traits_->getPeakPos(*it);
			file2 << p[RT] << " " << p[MZ] << " " << traits_->getPeakIntensity(*it) << "\n";						
		}
		file2.close();
		#endif
		++counter_;

	  return f;
	}

}



