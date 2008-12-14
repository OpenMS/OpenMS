// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
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
	    	@brief Applies the resampling algorithm to to an given iterator range.
	
	      Creates uniform data from the raw data given iterator intervall [first,last) and writes the 
	      resulting data to the resampled_peak_container.
	      
	      @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<Peak1D >::const_iterator)
	            points to a data point of type Peak1D or any other class derived from Peak1D.
	      
	            The resulting raw data in the resampled_peak_container (e.g. of type MSSpectrum<Peak1D >)
	            can be of type Peak1D or any other class derived from Peak1D.
	    */
	    template < typename InputPeakIterator, typename OutputPeakContainer >
	    void raster(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& resampled_peak_container)
      {
	        double end_pos = (last-1)->getMZ();
	        double start_pos = first->getMZ();
	        int number_raw_points = distance(first,last);
	        int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing_ + 1));
	
	        resampled_peak_container.resize(number_resampled_points);
	
	        // generate the resampled peaks at positions origin+i*spacing_
	        typename OutputPeakContainer::iterator it = resampled_peak_container.begin();
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
//								
//							std::cout << "Number of points " << number_raw_points << '\n'
//												<< "Act pos " << (first+i)->getPosition()[0] 
//												<< " Start pos " << start_pos << " End position " << end_pos << std::endl;

//	
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
	    }
	
	
	    /** 
	    	@brief Applies the resampling algorithm to a raw data point container.
	
	      Creates uniform data from the raw data in the input container (e.g. of type MSSpectrum<Peak1D >) and writes the 
	      resulting data to the resampled_peak_container.
	      
	      @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<Peak1D >::const_iterator)
	            points to a data point of type Peak1D or any other class derived from Peak1D.
	      
	            The resulting raw data in the resampled_peak_container (e.g. of type MSSpectrum<Peak1D >)
	            can be of type Peak1D or any other class derived from Peak1D.
	    */
	    template <typename InputPeakContainer, typename OutputPeakContainer >
	    void raster(const InputPeakContainer& input_peak_container, OutputPeakContainer& baseline_filtered_container)
	    {
        // copy the spectrum settings
        static_cast<SpectrumSettings&>(baseline_filtered_container) = input_peak_container;

	    	raster(input_peak_container.begin(), input_peak_container.end(), baseline_filtered_container);
	    }
	
	
	    /**
	    	@brief Resamples the data in a range of MSSpectra.
	
	    	Rasters the raw data successive in every scan in the intervall [first,last).
	    	The resampled data are stored in a MSExperiment.
	    		
	    	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectra should be of type Peak1D 
	            or any other derived class of Peak1D.
	
	     	@note You have to copy the ExperimentalSettings of the raw data on your own. 	
	    */
	    template <typename InputSpectrumIterator, typename OutputPeakType >
	    void rasterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_filtered)
	    {
	        UInt n = distance(first,last);
          ms_exp_filtered.reserve(n);
          startProgress(0,n,"resampling of data");
	        // pick peaks on each scan
	        for (UInt i = 0; i < n; ++i)
	        {
	            MSSpectrum< OutputPeakType > spectrum;
	            InputSpectrumIterator input_it(first+i);
	
	            // pick the peaks in scan i
	            raster(*input_it,spectrum);
              setProgress(i);
	
	            // if any peaks are found copy the spectrum settings
	            if (spectrum.size() > 0)
	            {
	                // copy the spectrum settings
	                static_cast<SpectrumSettings&>(spectrum) = *input_it;
	                spectrum.setType(SpectrumSettings::RAWDATA);
	
	                // copy the spectrum information
	                spectrum.getPrecursorPeak() = input_it->getPrecursorPeak();
	                spectrum.setRT(input_it->getRT());
	                spectrum.setMSLevel(input_it->getMSLevel());
	                spectrum.getName() = input_it->getName();
	
	                ms_exp_filtered.push_back(spectrum);
	            }
	        }
          endProgress();
	    }
	
	
	
	    /** 
	    	@brief Resamples the data in an MSExperiment.
	
	    	Rasters the raw data of every scan in the MSExperiment.
	    	The resampled data are stored in a MSExperiment.	
	    				
	    	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectra should be of type Peak1D 
	            or any other derived class of Peak1D.
	    */
	    template <typename InputPeakType, typename OutputPeakType >
	    void rasterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)
	    {
        // copy the experimental settings
        static_cast<ExperimentalSettings&>(ms_exp_filtered) = ms_exp_raw;

        rasterExperiment(ms_exp_raw.begin(), ms_exp_raw.end(), ms_exp_filtered);
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
