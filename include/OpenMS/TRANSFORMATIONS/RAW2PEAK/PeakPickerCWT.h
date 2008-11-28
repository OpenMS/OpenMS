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
//
#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <math.h>
#include <vector>
#include <algorithm>


#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
		@brief This class implements a peak picking algorithm using wavelet techniques
			
		The algorithm is descripted in detail in Lange et al. (2006) Proc. PSB-06.
			
		This peak picking algorithm uses the continuous wavelet transform of a raw data signal to detect mass peaks.
		Afterwards a given asymmetric peak function is fitted to the raw data and important peak parameters (e.g. fwhm)
		are extracted.
		In an optional step these parameters can be optimized using a non-linear opimization method.
			
		The peak parameters are stored in the meta data arrays of the spectra (see DSpectrum) in this order:
		- rValue
		- area
		- fwhm
		- leftWidth
		- rightWidth
		- peakShape
		- SignalToNoise
		.
	
		@htmlinclude OpenMS_PeakPickerCWT.parameters
	  
		@ingroup PeakPicking
  */
  class PeakPickerCWT
  	: public DefaultParamHandler, 
  		public ProgressLogger
  {
	 public:

    /// Raw data point type
    typedef Peak1D PeakType;
    /// Raw data container type using for the temporary storage of the input data
    typedef std::vector<PeakType> RawDataArrayType;
    /// Raw data iterator type
    typedef RawDataArrayType::iterator PeakIterator;
    /// Position type
    typedef DPosition<1> PositionType;

    /// Constructor
    PeakPickerCWT();

    /// Destructor
    virtual ~PeakPickerCWT();

    /** 
			@brief Applies the peak picking algorithm to an given iterator range.
	        
			Picks the peaks in the given iterator intervall [first,last) and writes the
			resulting peaks to the picked_peak_container.
			The ms_level should be one if the spectrum is a normal mass spectrum, or two if it is a tandem mass spectrum.
	        
			@note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<Peak1D>::const_iterator)
			points to a data point of type Peak1D or any other class derived from Peak1D.
	        
			@note The resulting peaks in the picked_peak_container (e.g. of type MSSpectrum<>)
			can be of type Peak1D or any other class derived from Peak1D.
    */
    template <typename InputPeakIterator, typename OutputPeakContainer  >
    void pick(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& picked_peak_container, int ms_level = 1)
    {

      if (peak_bound_cwt_==0.0 || peak_bound_ms2_level_cwt_==0.0)
      {
        initializeWT_();
      }

      // empty spectra shouldn't be picked
      if(first == last)
			{
        return;
			}

      typedef typename OutputPeakContainer::value_type OutputPeakType;

			//prepare the container
			picked_peak_container.clear();
      //clear the peak shapes vector
      peak_shapes_.clear();
      //set up meta data arrays
      picked_peak_container.getMetaDataArrays().clear();
      picked_peak_container.getMetaDataArrays().resize(7);
      picked_peak_container.getMetaDataArrays()[0].setName("rValue");
      picked_peak_container.getMetaDataArrays()[1].setName("maximumIntensity");
      picked_peak_container.getMetaDataArrays()[2].setName("fwhm");
      picked_peak_container.getMetaDataArrays()[3].setName("leftWidth");
      picked_peak_container.getMetaDataArrays()[4].setName("rightWidth");
      picked_peak_container.getMetaDataArrays()[5].setName("peakShape");
      picked_peak_container.getMetaDataArrays()[6].setName("SignalToNoise");

#ifdef DEBUG_PEAK_PICKING
      std::cout << "****************** PICK ******************" << std::endl;
#endif

      // vector of peak endpoint positions
      std::vector<double> peak_endpoints;

      // copy the raw data into a std::vector<Peak1D>
			//  raw_peak_array_original is needed for the separation of overlapping peaks
      RawDataArrayType raw_peak_array, raw_peak_array_original;
      // signal to noise estimator
			SignalToNoiseEstimatorMeanIterative< RawDataArrayType > sne;      
      Param sne_param(param_.copy("SignalToNoiseEstimationParameter:",true));
      sne.setParameters(sne_param);
      
			UInt n = distance(first, last);
      raw_peak_array.resize(n);
			raw_peak_array_original.resize(n);
      
      for (UInt i = 0; i < n; ++i)
      {
        PeakType raw_data_point;
        raw_data_point.setIntensity((first + i)->getIntensity());
        raw_data_point.setPosition((first + i)->getPosition());
        raw_peak_array[i] = raw_data_point;
        raw_peak_array_original[i] = raw_data_point;
      }


      PeakIterator it_pick_begin = raw_peak_array.begin();
      PeakIterator it_pick_end   = raw_peak_array.end();

			std::set<double> ints;
			for (PeakIterator it = raw_peak_array.begin(); it != raw_peak_array.end(); ++it)
			{
				ints.insert(it->getIntensity());
			}
			if (ints.size() < 2)
			{
				return;
			}
			
      if(it_pick_begin == it_pick_end)
        return;
      sne.init(it_pick_begin,it_pick_end);

      // thresholds for deconvolution
      double fwhm_threshold = (float)param_.getValue("deconvolution:fitting:fwhm_threshold");
      double symm_threshold = (float)param_.getValue("deconvolution:asym_threshold");


      // Points to the actual maximum position in the raw data
      PeakIterator it_max_pos;

      // start the peak picking until no more maxima can be found in the wavelet transform
      UInt number_of_peaks = 0;

      do
      {
        number_of_peaks = 0;
        int peak_left_index, peak_right_index;

        // compute the continious wavelet transform with resolution 1
        DoubleReal resolution = 1;
        wt_.transform(it_pick_begin, it_pick_end,resolution);
        PeakArea_ area;
        bool centroid_fit=false;
        bool regular_endpoints=true;

        // search for maximum positions in the cwt and extract potential peaks
        int direction=1;
        int distance_from_scan_border = 0;
        while ((distance(it_pick_begin, it_pick_end) > 3)
               && getMaxPosition_(it_pick_begin,
                                  it_pick_end,
                                  wt_,
                                  area,
                                  distance_from_scan_border,
                                  ms_level,
                                  direction))
        {
          // if the signal to noise ratio at the max position is too small
          // the peak isn't considered

          if((area.max  != it_pick_end) && (sne.getSignalToNoise(area.max) < signal_to_noise_) )
          {
            it_pick_begin = area.max;
            distance_from_scan_border = distance(raw_peak_array.begin(),it_pick_begin);
            
            continue;
          }
          else if(area.max >= it_pick_end) break;

          //search for the endpoints of the peak
          regular_endpoints = getPeakEndPoints_(it_pick_begin,
                                                it_pick_end,
                                                area,
                                                distance_from_scan_border,
                                                peak_left_index,
                                                peak_right_index);

          // compute the centroid position
          getPeakCentroid_(area);

          // if the peak achieves a minimal width, start the peak fitting
          if (regular_endpoints)
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "The endpoints are "
											<< area.left->getPosition()
											<< " and "
											<< area.right->getPosition()
											<< std::endl;
#endif
            // determine the best fitting lorezian or sech2 function
            PeakShape shape = fitPeakShape_(area,centroid_fit);
            shape.setLeftEndpoint( (raw_peak_array_original.begin() + distance(raw_peak_array.begin(), area.left)));
            shape.setRightEndpoint ( (raw_peak_array_original.begin() + distance(raw_peak_array.begin(), area.right)));
						if(shape.getRightEndpoint() == raw_peak_array_original.end()) shape.setRightEndpoint(raw_peak_array_original.end()-1); 
            // Use the centroid for Optimization
            shape.mz_position=area.centroid_position[0];
            if ( (shape.r_value > peak_corr_bound_)
                 && (shape.getFWHM() >= fwhm_bound_))
            {
              shape.signal_to_noise = sne.getSignalToNoise(area.max);
							peak_shapes_.push_back(shape);
              ++number_of_peaks;
            }

            else
            {
#ifdef DEBUG_PEAK_PICKING
              std::cout << "Corr: " << shape.r_value << " SN: " << sne.getSignalToNoise(area.max) << " FWHM: " << shape.getFWHM() << std::endl;
              std::cout << "Bad fitting peak "<< std::endl;
#endif
            }
          }

          // remove the peak from the signal
          // TODO: does this work as expected???
          for (PeakIterator pi=area.left; pi!=area.right+1; pi++)
          {
            pi->setIntensity(0.);
          }

          // search for the next peak
          it_pick_begin = area.right;
          distance_from_scan_border = distance(raw_peak_array.begin(),it_pick_begin);

        } //end while (getMaxPosition_(it_pick_begin, it_pick_end, wt_, area, distance_from_scan_border, ms_level, direction))
        it_pick_begin = raw_peak_array.begin();
      }
      while (number_of_peaks != 0);

      // start the nonlinear optimization for all peaks in split
#ifdef DEBUG_PEAK_PICKING
      std::cout << "Try the optimization run... with " << peak_shapes_.size() << std::endl;
#endif

      if (peak_shapes_.size() > 0)
      {
				// overlapping peaks are mostly broad or asymmetric
        // we distinguish them from broad or asymmetric isotopic peaks 
        // (e.g. charge one peaks, or peaks in the high mass range)
        // by a simple heuristic: if the distances to adjacent peaks
				// are dissimilar, the fhwm is much broader than the fhwm of
				// adjacent peaks or if the peak has no near neighbors
				// we assume a convolved peak pattern and start the deconvolution.
        // sort the peaks according to their positions
        sort(peak_shapes_.begin(), peak_shapes_.end(), PeakShape::PositionLess());
				std::vector<UInt> peaks_to_skip;
        // search for broad or asymmetric peaks
        UInt n = peak_shapes_.size();
				if( deconvolution_)
				{
					for (UInt i = 0; i < n; ++i)
					{
						if ((peak_shapes_[i].getFWHM() > fwhm_threshold) 
								|| (peak_shapes_[i].getSymmetricMeasure() < symm_threshold))
						{
#ifdef DEBUG_DECONV
							std::cout << "check " << peak_shapes_[i].mz_position 
												<< " with fwhm: " << peak_shapes_[i].getFWHM() 
												<< " and " << peak_shapes_[i].left_width 
												<< ' ' << peak_shapes_[i].right_width 
												<< ' ' << peak_shapes_[i].getLeftEndpoint()->getMZ()
												<< ' ' << peak_shapes_[i].getRightEndpoint()->getMZ()
												<< std::endl;
#endif
							// this might be a convolved peak pattern
							// and we check the distance to the neighboring peaks as well as their fwhm values:
							//float max_distance = 1.1;
							float dist_left = ((i > 0) && (fabs(peak_shapes_[i].mz_position-peak_shapes_[i-1].mz_position) < 1.2)) ? fabs(peak_shapes_[i].mz_position-peak_shapes_[i-1].mz_position) : -1;
							float dist_right = ((i < (n-1)) && (fabs(peak_shapes_[i].mz_position-peak_shapes_[i+1].mz_position) < 1.2)) ? fabs(peak_shapes_[i].mz_position-peak_shapes_[i+1].mz_position) : -1;
          
							// left and right neighbor
							if ((dist_left > 0) && (dist_right > 0))
							{
								// if distances to left and right adjacent peaks is dissimilar deconvolute
								DoubleReal ratio = (dist_left > dist_right) ? dist_right/dist_left : dist_left/dist_right;
#ifdef DEBUG_DECONV
								std::cout << "Ratio " << ratio << std::endl;
#endif
								if (ratio < 0.6)
								{
#ifdef DEBUG_DECONV
									std::cout << "deconvolute: dissimilar left and right neighbor "  << peak_shapes_[i-1].mz_position << ' ' << peak_shapes_[i+1].mz_position << std::endl;
#endif
									if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
								}
							}
							// has only one or no neighbor peak
							else
							{
								// only left neighbor
								if (dist_left > 0)
								{
									// check distance and compare fwhm
									DoubleReal dist = 1.00235;
									//check charge 1 or 2
									bool dist_ok = ((fabs(dist-dist_left) < 0.21) || (fabs(dist/2.-dist_left) < 0.11)) ? true : false ;  
									// distance complies peptide mass rule
									if (dist_ok)
									{
#ifdef DEBUG_DECONV
										std::cout << "left neighbor " << peak_shapes_[i-1].mz_position << ' ' << peak_shapes_[i-1].getFWHM() << std::endl;
#endif
										// if the left peak has a fwhm which is smaller than 60% of the fwhm of the broad peak deconvolute
										if ((peak_shapes_[i-1].getFWHM()/peak_shapes_[i].getFWHM()) < 0.6)
										{
#ifdef DEBUG_DECONV
											std::cout << " too small fwhm" << std::endl;
#endif
											if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
										}
									}
									else
									{
#ifdef DEBUG_DECONV
										std::cout << "distance not ok" << dist_left << ' ' << peak_shapes_[i-1].mz_position << std::endl;
#endif
										if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
									}
								}
								else
								{ 
									// only right neighbor 
									if (dist_right > 0)
									{
										// check distance and compare fwhm
										DoubleReal dist = 1.00235;
										//check charge 1 or 2
										bool dist_ok = ((fabs(dist-dist_right) < 0.21) || (fabs(dist/2.-dist_right) < 0.11)) ? true : false ;  
										// distance complies peptide mass rule
										if (dist_ok)
										{
#ifdef DEBUG_DECONV
											std::cout << "right neighbor " << peak_shapes_[i+1].mz_position << ' ' << peak_shapes_[i+1].getFWHM() << std::endl;
#endif
											// if the left peak has a fwhm which is smaller than 60% of the fwhm of the broad peak deconvolute
											if ((peak_shapes_[i+1].getFWHM()/peak_shapes_[i].getFWHM()) < 0.6)
											{
#ifdef DEBUG_DECONV
												std::cout << "too small fwhm"  << std::endl;
#endif
												if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
											}
										}
										else
										{
#ifdef DEBUG_DECONV
											std::cout << "distance not ok" << dist_right << ' ' << peak_shapes_[i+1].mz_position << std::endl;
#endif
											if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
										}
									}
									// no neighbor
									else
									{
#ifdef DEBUG_DECONV
										std::cout << "no neighbor" << std::endl;
#endif
										if(deconvolutePeak_(peak_shapes_[i])) peaks_to_skip.push_back(i);
									} 
								}
							}
						}
					}
				}
				
        // write the picked peaks to the outputcontainer
        for (UInt i = 0; i < peak_shapes_.size(); ++i)
        {
					// put it out only if the peak was not deconvoluted
					if(find(peaks_to_skip.begin(),peaks_to_skip.end(),i) == peaks_to_skip.end() )
					{
						//store output peak
						OutputPeakType picked_peak;						
						picked_peak.setIntensity(peak_shapes_[i].height);
						picked_peak.setMZ(peak_shapes_[i].mz_position);
						picked_peak_container.push_back(picked_peak);
						//store meta data
						picked_peak_container.getMetaDataArrays()[0].push_back(peak_shapes_[i].r_value);
				    picked_peak_container.getMetaDataArrays()[1].push_back(peak_shapes_[i].area);
				    picked_peak_container.getMetaDataArrays()[2].push_back(peak_shapes_[i].getFWHM());
				    picked_peak_container.getMetaDataArrays()[3].push_back(peak_shapes_[i].left_width);
				    picked_peak_container.getMetaDataArrays()[4].push_back(peak_shapes_[i].right_width);
				    picked_peak_container.getMetaDataArrays()[5].push_back(peak_shapes_[i].type);
				    picked_peak_container.getMetaDataArrays()[6].push_back(peak_shapes_[i].signal_to_noise);
					}
        }
      } // if (peak_shapes_.size() > 0)
			
			// set MS level of output container to match the input container
			picked_peak_container.setMSLevel(ms_level);
			
    }

    /** 
			@brief Applies the peak picking algorithm to a raw data point container.
	        
			Picks the peaks in the input container (e.g. of type MSSpectrum<Peak1D >) 
			and writes the resulting peaks to the picked_peak_container (e.g. MSSpectrum<>).
	
			The ms_level should be one if the spectrum is a normal mass spectrum, or two if it is a tandem mass spectrum.
	        
			@note This method assumes that the input_peak_container contains data points of type 
			Peak1D or any other class derived from Peak1D. 
	              
			@note The resulting peaks in the picked_peak_container (e.g. of type MSSpectrum<>)
			can be of type Peak1D or any other class derived from Peak1D.
    */
    template <typename InputPeakContainer, typename OutputPeakContainer >
    void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container, int ms_level = 1)
    {
      // copy the spectrum settings
      static_cast<SpectrumSettings&>(picked_peaks_container) = input_peak_container;
      
      pick(input_peak_container.begin(), input_peak_container.end(), picked_peaks_container, ms_level);
    }


    /** 
			@brief Picks the peaks in a range of MSSpectra.
	        
			Picks the peaks successive in every scan in the intervall [first,last).
			The detected peaks are stored in a MSExperiment.
	              
			@note The InputSpectrumIterator should point to a MSSpectrum.
			Elements of the input spectra should be of type Peak1D 
			or any other derived class of Peak1D.
	
			@note You have to copy the ExperimentalSettings of the raw data by your own.  
    */
    template <typename InputSpectrumIterator, typename OutputPeakType >
    void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_peaks)
    {
      UInt n = distance(first,last);
      ms_exp_peaks.reserve(n);
      startProgress(0,n,"picking peaks");
      // pick peaks on each scan
      for (UInt i = 0; i < n; ++i)
      {
				setProgress(i);
        MSSpectrum< OutputPeakType > spectrum;
        InputSpectrumIterator input_it(first+i);
#ifdef DEBUG_PEAK_PICKING
        std::cout << "PeakPicker: Picking Scan " << input_it->getRT()<< std::endl;
#endif
        // pick the peaks in scan i
        pick(*input_it,spectrum,input_it->getMSLevel());
        setProgress(i);

        // if any peaks are found copy the spectrum settings
        if (spectrum.size() > 0)
        {
          // copy the spectrum settings
          static_cast<SpectrumSettings&>(spectrum) = *input_it;
          spectrum.setType(SpectrumSettings::PEAKS);

          // copy the spectrum information
          spectrum.setPrecursorPeak(input_it->getPrecursorPeak());
          spectrum.setRT(input_it->getRT());
          spectrum.setMSLevel(input_it->getMSLevel());
          spectrum.getName() = input_it->getName();

          ms_exp_peaks.push_back(spectrum);
        }
      }
      // sort spectra
      ms_exp_peaks.sortSpectra(true);

      if(two_d_optimization_ || optimization_)
      {
				Param two_d_param(param_.copy("optimization:",true));
			
				TwoDOptimization my_2d;
				my_2d.setParameters(two_d_param);
        my_2d.optimize(first,last,ms_exp_peaks,two_d_optimization_);
      }
      // sort spectra
      ms_exp_peaks.sortSpectra(true);
      endProgress();
    }

    /** 
			@brief Picks the peaks in a MSExperiment.
	        
			Picks the peaks on every scan in the MSExperiment.
			The detected peaks are stored in a MSExperiment.
	              
			@note The input peaks should be of type Peak1D or any other derived class of Peak1D.
    */
    template <typename InputPeakType, typename OutputPeakType >
    void pickExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_peaks)
    {
      // copy the experimental settings
      static_cast<ExperimentalSettings&>(ms_exp_peaks) = ms_exp_raw;

      pickExperiment(ms_exp_raw.begin(),ms_exp_raw.end(),ms_exp_peaks);
    }

	 protected:

    /// Threshold for the peak height in the MS 1 level
    float peak_bound_;

    /// Threshold for the peak height in the MS 2 level
    float peak_bound_ms2_level_;

    /// Signal to noise threshold
    float signal_to_noise_;

    /// The minimal full width at half maximum
    float fwhm_bound_;

    /// Container the determined peak shapes
    std::vector<PeakShape> peak_shapes_;

    /// The continuous wavelet "transformer"
    ContinuousWaveletTransformNumIntegration wt_;

    /// The continuous wavelet "transformer" for the deconvolution
    ContinuousWaveletTransformNumIntegration wtDC_;

    /// The search radius for the determination of a peak's maximum position
    UInt radius_;

    /// The dilation of the wavelet
    float scale_;

    /// The minimal height which defines a peak in the CWT (MS 1 level)
    float peak_bound_cwt_;

    /// The minimal height which defines a peak in the CWT (MS 2 level)
    float peak_bound_ms2_level_cwt_;

    /// The threshold for correlation
    float peak_corr_bound_;

    /// The threshold for the noise level (TODO: Use the information of the signal to noise estimator)
    float noise_level_;

    /// Switch for the optimization of peak parameters
    bool optimization_;

    /// Switch for the deconvolution of peak parameters
    bool deconvolution_;

    /// Switch for the 2D optimization of peak parameters
    bool two_d_optimization_;


		void updateMembers_();

    /// Initializes the members and parses the parameter object
    void init_();


    /**
			@brief Class for the internal peak representation
	        
			A regularData-Object which contents some additional useful informations
			for analysing peaks and their properties
    */
    class PeakArea_
    {
      typedef std::vector<PeakType>::iterator PeakIterator;

		 public:
      PeakArea_() : left(), max(), right(), left_behind_centroid()
      {
      }

      /** 
			@brief Iterator defining a raw data peak.
         
			The left and right iterators delimit a range in the raw data which represents a raw peak.
			They define the raw peak endpoints. Max points to the raw data point in [left, right] with the highest intensity, the 
			maximum of the raw peak. 
				
			Left_behind_centroid points to the raw data point next to the estimates centroid position.
      */
      PeakIterator left;
      PeakIterator max;
      PeakIterator right;
      PeakIterator left_behind_centroid;
      /// The estimated centroid position.
      DPosition<1> centroid_position;
    };

    /// Computes the peak's left and right area
    void getPeakArea_(const PeakArea_& area, double &area_left, double &area_right);

    /// Returns the best fitting peakshape
    PeakShape fitPeakShape_(const PeakArea_& area, bool enable_centroid_fit);

    /** 
			@brief Returns the squared pearson coefficient.
	
			Computes the correlation of the peak and the original data given by the peak enpoints area.left and area.right.
			If the value is near 1, the fitted peakshape and the raw data are expected to be very similar. 
    */
    double correlate_(const PeakShape& peak, const PeakArea_& area, int direction=0) const;


    /** 
			@brief Finds the next maximum position in the wavelet transform wt.
	        
			If the maximum is greater than peak_bound_cwt we search for the corresponding maximum in the raw data interval [first,last)
			given a predefined search radius radius. Only peaks with intensities greater than peak_bound_ 
			are relevant. If no peak is detected the method return false.
			For direction=1, the method runs from first to last given direction=-1 it runs the other way around.
    */
    bool getMaxPosition_(PeakIterator first, PeakIterator last, const ContinuousWaveletTransform& wt, PeakArea_& area, int distance_from_scan_border, int ms_level, int direction=1);


    /** 
			@brief Determines a peaks's endpoints.
	      
			The algorithm does the following:
			- let x_m be the position of the maximum in the data and let (x_l, x_r) be
				the left and right neighbours
	    -	(1) starting from x_l', walk left until one of the following happens
				- the new point is lower than the original bound => we found our left endpoint
				- the new point is larger than the last, but the point left from
					the new point is smaller. In that case, we either ran into another
					peak, or we encounter some noise. Therefore we now look in the cwt
					at the position corresponding to this value. If the cwt here is
					monotonous, we consider the point as noise and continue further to the
					left. Otherwise, we probably found the beginning of a new peak and
					therefore stop here.
				.
			-	(2) analogous procedure to the right of x_r
			.
    */
    bool getPeakEndPoints_(PeakIterator first, PeakIterator last,  PeakArea_ &area, int distance_from_scan_border, int& peak_left_index, int& peak_right_index);


    /** 
			@brief Estimates a peak's centroid position.
	
			Computes the centroid position of the peak using all raw data points which are greater than 
			60% of the most intensive raw data point.
    */
    void getPeakCentroid_(PeakArea_& area);

    /// Computes the value of a theroretical lorentz peak at position x
    double lorentz_(double height, double lambda, double pos, double x);

    /** 
			@brief Computes the threshold for the peak height in the wavelet transform and initializes the wavelet transform.
	
			Given the threshold for the peak height a corresponding value peak_bound_cwt can be computed
			for the continious wavelet transform. 
			Therefore we compute a theoretical lorentzian peakshape with height=peak_bound_ and a width which 
			is similar to the width of the wavelet. Taking the maximum in the wavelet transform of the
			lorentzian peak we have a peak bound in the wavelet transform. 
    */
    void initializeWT_();

		/** @name Methods needed for separation of overlapping peaks
    */
    //@{
		
		/** 
			@brief Separates overlapping peaks.
	
			It determines the number of peaks lying underneath the initial peak using the cwt with different scales.
			Then a nonlinear optimzation procedure is applied to optimize the peak parameters.
		*/
    bool deconvolutePeak_(PeakShape& shape);

		/// Determines the number of peaks in the given mass range using the cwt
    int getNumberOfPeaks_(PeakIterator first,PeakIterator last, std::vector<double>& peak_values,
													int direction,DoubleReal resolution, ContinuousWaveletTransformNumIntegration& wt);

		/// Estimate the charge state of the peaks
    int determineChargeState_(std::vector<double>& peak_values);

		/// Add a peak
    void addPeak_(std::vector<PeakShape>& peaks_DC,PeakArea_& area,double left_width,double right_width);
		//@}
  }
		; // end PeakPickerCWT


}// namespace OpenMS

#endif
