// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H
#define OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <limits>
#include <cmath>

namespace OpenMS
{
	/**
		@brief Linear Resampling of raw data.
		
		This class can be used to generate uniform data from non-uniform raw data (e.g. ESI-TOF or MALDI-TOF experiments).
		Therefore the intensity at every position x in the input raw data is spread to the two
		adjacent resampling points.
		This method preserves the area of the input signal and also the centroid position of a peak.
		Therefore it is recommended for quantitation as well as for ProteinIdentification experiments.
		
		@note Use this method only for high resoluted data (< 0.1 Th between two adjacent raw data points).
		     The resampling rate should be >= the precision.
		 
		@htmlinclude OpenMS_LinearResampler.parameters
	*/
	class OPENMS_DLLAPI LinearResampler 
		: public DefaultParamHandler, 
			public ProgressLogger
	{
	
	public:
	
	    /// Constructor
	    LinearResampler()
      	: DefaultParamHandler("LinearResampler")
      {
        defaults_.setValue("spacing",0.05,"Spacing of the resampled output peaks.");
        defaultsToParam_();
	    }
		
	    /// Destructor.
	    ~LinearResampler()
	    {
	    }
	
	    /** 
	    	@brief Applies the resampling algorithm to to an MSSpectrum.
	    */
	    template < typename PeakType >
	    void raster(MSSpectrum<PeakType>& spectrum)
      {
      	//return if nothing to do
      	if (spectrum.size()==0) return;
      	
      	typename MSSpectrum<PeakType>::iterator first = spectrum.begin();
      	typename MSSpectrum<PeakType>::iterator last = spectrum.end();
      	
        double end_pos = (last-1)->getMZ();
        double start_pos = first->getMZ();
        int number_raw_points = (int)spectrum.size();
        int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing_ + 1));

        typename std::vector<PeakType> resampled_peak_container;
        resampled_peak_container.resize(number_resampled_points);

        // generate the resampled peaks at positions origin+i*spacing_
        typename std::vector<PeakType>::iterator it = resampled_peak_container.begin();
        for (int i=0; i < number_resampled_points; ++i)
        {
            it->setMZ( start_pos + i*spacing_);
            ++it;
        }

        // spread the intensity h of the data point at position x to the left and right
        // adjacent resampled peaks
        double distance_left = 0.;
        double distance_right = 0.;
        int left_index = 0;
        int right_index = 0;

        it = resampled_peak_container.begin();
        for (int i=0; i < number_raw_points ; ++i)
        {
	          int help = (int)floor(((first+i)->getMZ() - start_pos) / spacing_);
            left_index = (help < 0) ? 0 : help;
            help = distance(first,last) - 1;
            right_index = (left_index >= help) ? help : left_index + 1;

            // compute the distance between x and the left adjacent resampled peak
            distance_left = fabs((first+i)->getMZ() - (it + left_index)->getMZ()) / spacing_;
            //std::cout << "Distance left " << distance_left << std::endl;
            // compute the distance between x and the right adjacent resampled peak
            distance_right = fabs((first+i)->getMZ() - (it + right_index)->getMZ());
            //std::cout << "Distance right " << distance_right << std::endl;
            

            // add the distance_right*h to the left resampled peak and distance_left*h to the right resampled peak
            DoubleReal intensity = (it + left_index)->getIntensity();
            intensity += (first+i)->getIntensity()*distance_right / spacing_;
            (it + left_index)->setIntensity(intensity);
            intensity = (it + right_index)->getIntensity();
            intensity += (first+i)->getIntensity()*distance_left;
            (it + right_index)->setIntensity(intensity);
        }
        
        resampled_peak_container.swap(spectrum);
	    }
	
	    /**
	    	@brief Resamples the data in an MSExperiment.
	    */
	    template <typename PeakType >
	    void rasterExperiment(MSExperiment<PeakType>& exp)
			{
				startProgress(0,exp.size(),"resampling of data");
				for (Size i = 0; i < exp.size(); ++i)
				{
					raster(exp[i]);
					setProgress(i);
				}
				endProgress();
			}
	
		protected:
		
	    /// Spacing of the resampled data
	    double spacing_;
      
      virtual void updateMembers_()
      {
        spacing_ =  param_.getValue("spacing");
      }
	};


} // namespace OpenMS

#endif // OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H
