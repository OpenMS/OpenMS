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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PickedPeakSeeder.h>

using namespace std;

namespace OpenMS
{

	PickedPeakSeeder::PickedPeakSeeder()
		: BaseSweepSeeder(),
			isomodel_()
	{
    setName(getProductName());

    // lower and upper bounds for distances between isotopic peaks (defaults)
    // charge 1
    defaults_.setValue("charge1_ub",1.6f);
    defaults_.setValue("charge1_lb",0.7f);
    // charge 2
    defaults_.setValue("charge2_ub",0.7f);
    defaults_.setValue("charge2_lb",0.4f);
    // charge 3
    defaults_.setValue("charge3_ub",0.4f);
    defaults_.setValue("charge3_lb",0.27f);
    // charge 4
    defaults_.setValue("charge4_ub",0.27f);
    defaults_.setValue("charge4_lb",0.24f);
    // charge 5
    defaults_.setValue("charge5_ub",0.24f);
    defaults_.setValue("charge5_lb",0.15f);

    // minimum number of peaks per pattern (in one scan)
    defaults_.setValue("min_peaks_per_scan",3);

    defaultsToParam_();
	}
	
	PickedPeakSeeder::~PickedPeakSeeder()
	{
	}

  PickedPeakSeeder::PickedPeakSeeder(const PickedPeakSeeder& rhs)
    : BaseSweepSeeder(rhs)
  {
    updateMembers_();
  }
  
  PickedPeakSeeder& PickedPeakSeeder::operator= (const PickedPeakSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSweepSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
	
	void PickedPeakSeeder::updateMembers_()
	{	
		// update member of base class first
		BaseSweepSeeder::updateMembers_();
	
		// retrieve values for accepted peaks distances
		charge1_ub_	= param_.getValue("charge1_ub");
		charge1_lb_	 = param_.getValue("charge1_lb");
	
		charge2_ub_	= param_.getValue("charge2_ub");
		charge2_lb_	 = param_.getValue("charge2_lb");
	
		charge3_ub_	= param_.getValue("charge3_ub");
		charge3_lb_	 = param_.getValue("charge3_lb");
	
		charge4_ub_	= param_.getValue("charge4_ub");
		charge4_lb_	 = param_.getValue("charge4_lb");
	
		charge5_ub_	= param_.getValue("charge5_ub");
		charge5_lb_	 = param_.getValue("charge5_lb");		
		
		min_peaks_ = param_.getValue("min_peaks_per_scan");
	}
	
	PickedPeakSeeder::ScoredMZVector PickedPeakSeeder::detectIsotopicPattern_(SpectrumType& scan )
	{
		UInt current_charge	 = 0;			// charge state of the current isotopic cluster
		
		ScoredMZVector scored_positions;
		
		vector<IntensityType> data_intensities;
		vector<IntensityType> model_intensities;
	
		// it would be better to look at more peaks not only the next neighbour.
		for (UInt j=0; j < (scan.size()-1); ++j)
		{
				CoordinateType dist 	= scan[j+1].getMZ() - scan[j].getMZ();
						
				// test for different charge states
				current_charge = distanceToCharge_(dist);

				// remove false positives by looking at the intensity ratio of the first two peaks peaks
				// their intensities should not be equal.
				if ( fabs( scan[j].getIntensity()/scan[j+1].getIntensity() - 1.0) < 0.001)
				{					
					current_charge = 0;	// reset charge
				}
								
				if (current_charge > 0) 	// 0 => no pattern
				{
					ScoredChargeType sc_charge;
					sc_charge.first = current_charge;
					
					// initialize averagine model
					Param p;
					p.setValue("statistics:mean", (scan[j].getMZ() + 1) );
					p.setValue("isotope:stdev",0.1);	
					isomodel_.setParameters(p);
					isomodel_.setSamples();
									
					data_intensities.push_back( scan[j].getIntensity() );	
					model_intensities.push_back( isomodel_.getIntensity( scan[j].getMZ() ) );
					
					// count number of peaks supporting this charge
					UInt count = 1; 
					for (UInt c = j + 1; c < scan.size(); ++c)
					{	
						CoordinateType this_mz = 	scan[c].getMZ();
						CoordinateType prev_mz = scan[ (c-1) ].getMZ();
					
						UInt next_charge = 	distanceToCharge_( (this_mz - prev_mz) );
						
						if (next_charge != current_charge) 	break;
					
						data_intensities.push_back( scan[c].getIntensity() );	
						model_intensities.push_back( isomodel_.getIntensity( scan[c].getMZ() ) );
						++count;	
						++j;
					}
													
					#ifdef DEBUG_FEATUREFINDER
					std::cout	<< "There are " << count << " peaks supporting this charge. " << std::endl;
					#endif
					
					if (count >= min_peaks_)  										
					{
						// use number of peaks as score
						sc_charge.second = scorePattern_(data_intensities, model_intensities);
						scored_positions.push_back( make_pair(j,sc_charge) );				
					}
				}
		
		} // end for (all peaks in scan)

		return scored_positions;
		
} // end of detectIsotopicPattern_(SpectrumType& scan )
	

PickedPeakSeeder::ProbabilityType PickedPeakSeeder::scorePattern_(std::vector<IntensityType>& data, std::vector<IntensityType>& model)
{
	IntensityType data_sum   = 0.0;
	IntensityType model_sum = 0.0;
	
	// normalize...
	for (UInt i=0;i<data.size(); ++i)
	{
		data_sum += data[i];	
		model_sum += model[i];	
	}
	
	// compute chi^2 statistic
	ProbabilityType chi_stat = 0.0;
	IntensityType temp = 0.0;
	for (UInt j=0;j<data.size();++j)
	{
			if (model[j] <= 0) continue; // skip zeros in model
		
			cout << "data: " << ( data[j] / data_sum) << " model: " << (model[j] / model_sum) << endl;
				
			temp = (data[j] / data_sum) - ( model[j] / model_sum);
			chi_stat += (temp * temp) / (model[j] /  model_sum) ;	
	}

	cout << "test statistic " << chi_stat << endl;
	cout << "p-value is " << (1 - gsl_cdf_chisq_P(chi_stat, (data.size() - 1 ))) << endl;
	
	return chi_stat;
}	
	
UInt PickedPeakSeeder::distanceToCharge_(CoordinateType dist)
	{
	  if (dist <= charge1_ub_ && dist >= charge1_lb_)
    {
    	return 1;
    }
    else if (dist <= charge2_ub_ && dist >= charge2_lb_)
    {
    	return 2;
    }
    else if (dist <= charge3_ub_ && dist >= charge3_lb_)
    {
    	return 3;
    }
    else if (dist <= charge4_ub_ && dist >= charge4_lb_)
    {
    	return 4;
    }
    else if (dist <= charge5_ub_ && dist >= charge5_lb_)
    {
    	return 5;
    }
    else
    {
    	return 0;
    }
	}

} // end of namespace OpenMS
