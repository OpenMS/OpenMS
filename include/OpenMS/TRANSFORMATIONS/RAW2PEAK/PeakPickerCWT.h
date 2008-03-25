// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
#include <OpenMS/KERNEL/PickedPeak1D.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>

#include <math.h>
#include <vector>
#include <algorithm>


#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
		@brief This class implements a peak picking algorithm using wavelet techniques (as described by Lange et al. (2006) Proc. PSB-06).
		
		This peak picking algorithm uses the continuous wavelet transform of a raw data signal to detect mass peaks.
		Afterwards a given asymmetric peak function is fitted to the raw data and important peak parameters (e.g. fwhm)
		are extracted.
		In an optional step these parameters can be optimized using a non-linear opimization method.
		 
		@ref PeakPickerCWT_Parameters are explained on a separate page.

    @ingroup PeakPicking
  */
  class PeakPickerCWT : public PeakPicker
  {
  public:

    /// Raw data point type
    typedef RawDataPoint1D RawDataPointType;
    /// Raw data container type using for the temporary storage of the input data
    typedef DPeakArray<RawDataPointType > RawDataArrayType;
    /// Raw data iterator type
    typedef RawDataArrayType::iterator RawDataPointIterator;
    /// Position type
    typedef DPosition<1> PositionType;


    using PeakPicker::peak_bound_;
    using PeakPicker::peak_bound_ms2_level_;
    using PeakPicker::signal_to_noise_;
    using PeakPicker::fwhm_bound_;

    /// Constructor
    PeakPickerCWT();

    /// Destructor
    virtual ~PeakPickerCWT();

		/** @name Accesssor methods
		 */
    //@{
    /// Non-mutable access to the vector of peak shapes
    inline const std::vector<PeakShape>& getPeakShapes() const
    {
      return peak_shapes_;
    }

    /// Non-mutable access to the wavelet transform
    inline const ContinuousWaveletTransformNumIntegration& getWaveletTransform() const
    {
      return wt_;
    }

    /// Non-mutable access to the search radius for the peak maximum
    inline UInt getSearchRadius() const
    {
      return radius_;
    }
    /// Mutable access to the search radius for the peak maximum
    inline void setSearchRadius(UInt radius)
    {
      radius_ = radius;
      param_.setValue("thresholds:search_radius",radius);
    }

    /// Non-mutable access to the scale of the wavelet transform
    inline Real getWaveletScale() const
    {
      return scale_;
    }
    /// Mutable access to the scale of the wavelet transform
    inline void setWaveletScale(Real scale)
    {
      scale_ = scale;
      param_.setValue("wavelet_transform:scale",scale);
    }

    /// Non-mutable access to the threshold of the height
    inline Real getPeakBound() const
    {
      return peak_bound_;
    }

    /// Non-mutable access to the peak bound in the wavelet transform for the MS 1 level
    inline Real getPeakBoundCWT() const
    {
      return peak_bound_cwt_;
    }
    /// Non-mutable access to the peak bound in the wavelet transform for the MS 2 level
    inline Real getPeakBoundMs2LevelCWT() const
    {
      return peak_bound_ms2_level_cwt_;
    }

    /// Non-mutable access to the minimum peak correlation coefficient
    inline Real getPeakCorrBound() const
    {
      return peak_corr_bound_;
    }
    /// Mutable access to the minimum peak correlation coefficient
    inline void setPeakCorrBound(Real peak_corr_bound)
    {
      peak_corr_bound_ = peak_corr_bound;
      param_.setValue("thresholds:correlation",peak_corr_bound);
    }

    /// Non-mutable access to the noise level
    inline Real getNoiseLevel() const
    {
      return noise_level_;
    }
    /// Mutable access to the noise level
    inline void setNoiseLevel(Real noise_level)
    {
      noise_level_ = noise_level;
      param_.setValue("thresholds:noise_level",noise_level);
    }

    /// Non-mutable access to the optimization switch
    inline const bool& getOptimizationFlag() const
    {
      return optimization_;
    }
    /// Mutable access to the optimization switch
    inline void setOptimizationFlag(const bool& optimization)
    {
      optimization_ = optimization;
      if (optimization)
      {
        param_.setValue("optimization:optimization","one_dimensional");
      }
      else
      {
        param_.setValue("optimization:optimization","no");
      }
    }

    /// Non-mutable access to the deconvolution switch
    inline const bool& getDeconvolutionFlag() const
    {
      return deconvolution_;
    }
    /// Mutable access to the deconvolution switch
    inline void setDeconvolutionFlag(const bool& deconvolution)
    {
      deconvolution_ = deconvolution;
      if (deconvolution)
      {
        param_.setValue("deconvolution:skip_deconvolution","no");
      }
      else
      {
        param_.setValue("deconvolution:skip_deconvolution","yes");
      }
    }
    /// Non-mutable access to the optimization switch
    inline const bool& get2DOptimizationFlag() const
    {
      return two_d_optimization_;
    }
    /// Mutable access to the optimization switch
    inline void set2DOptimizationFlag(const bool& two_d_optimization)
    {
      two_d_optimization_ = two_d_optimization;
      if (two_d_optimization)
      {
        param_.setValue("optimization:optimization","two_dimensional");
      }

    }
		//@}
    //@}


    /** @brief Applies the peak picking algorithm to an given iterator range.
        
      Picks the peaks in the given iterator intervall [first,last) and writes the
      resulting peaks to the picked_peak_container.
        The ms_level should be one if the spectrum is a normal mass spectrum, or two if it is a tandem mass spectrum.
        
      @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<RawDataPoint1D >::const_iterator)
            points to a data point of type RawDataPoint1D or any other class derived from RawDataPoint1D.
        
           The resulting peaks in the picked_peak_container (e.g. of type MSSpectrum<PickedPeak1D>)
           can be of type RawDataPoint1D or any other class derived from DRawDataPoint. 
           We recommend to use the PickedPeak1Dbecause it stores important information gained during
           the peak picking algorithm.
              
           If you use MSSpectrum iterators you have to set the SpectrumSettings on your own.              
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
			else
			{
				// only one "peak"
				if (last - first == 1)
				{
					// @todo handle spectra with only one raw data point (seg. faulted) (Chris Bielow)
					return;
				}
			}
      typedef typename OutputPeakContainer::value_type OutputPeakType;

      //clear the peak shapes vector
      peak_shapes_.clear();

#ifdef DEBUG_PEAK_PICKING

      std::cout << "****************** PICK ******************" << std::endl;
#endif

      // vector of peak endpoint positions
      std::vector<double> peak_endpoints;

      // copy the raw data into a DPeakArray<DRawDataPoint<D> >
			//  raw_peak_array_original is needed for the separation of overlapping peaks
      RawDataArrayType raw_peak_array, raw_peak_array_original;
      // signal to noise estimator
			SignalToNoiseEstimatorMeanIterative< RawDataArrayType > sne;      
      Param sne_param(param_.copy("SignalToNoiseEstimationParameter:",true));
      if(sne_param.empty()) sne.setParameters(Param());
      else sne.setParameters(sne_param);
      unsigned int n = distance(first, last);
      raw_peak_array.resize(n);
			raw_peak_array_original.resize(n);
      
      for (unsigned int i = 0; i < n; ++i)
      {
        RawDataPointType raw_data_point;
        raw_data_point.setIntensity((first + i)->getIntensity());
        raw_data_point.setPosition((first + i)->getPosition());
        raw_peak_array[i] = raw_data_point;
        raw_peak_array_original[i] = raw_data_point;
      }


      RawDataPointIterator it_pick_begin = raw_peak_array.begin();
      RawDataPointIterator it_pick_end   = raw_peak_array.end();

      if(it_pick_begin == it_pick_end)
        return;
      StopWatch timer;
      timer.start();
      sne.init(it_pick_begin,it_pick_end);
      timer.stop();
#ifdef DEBUG_PEAK_PICKING
      std::cout << "SNE init " << timer.getCPUTime() << std::endl;
#endif

      // thresholds for deconvolution
      double fwhm_threshold = (float)param_.getValue("deconvolution:fitting:fwhm_threshold");
      double symm_threshold = (float)param_.getValue("deconvolution:asym_threshold");


      // Points to the actual maximum position in the raw data
      RawDataPointIterator it_max_pos;

      // start the peak picking until no more maxima can be found in the wavelet transform
      unsigned int number_of_peaks = 0;

      do
      {
        number_of_peaks = 0;
        int peak_left_index, peak_right_index;

        // compute the continious wavelet transform with resolution 1
        timer.reset();
        timer.start();
        DoubleReal resolution = 1;
        wt_.transform(it_pick_begin, it_pick_end,resolution);
        timer.stop();
#ifdef DEBUG_PEAK_PICKING
        std::cout << "TRANSFORM " << timer.getCPUTime() << std::endl;
#endif
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
            shape.left_endpoint = (raw_peak_array_original.begin() + distance(raw_peak_array.begin(), area.left));
            shape.right_endpoint = (raw_peak_array_original.begin() + distance(raw_peak_array.begin(), area.right));
						if(shape.right_endpoint == raw_peak_array_original.end()) --shape.right_endpoint; 
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
          for (RawDataPointIterator pi=area.left; pi!=area.right+1; pi++)
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
															<< ' ' << peak_shapes_[i].left_endpoint->getMZ()
															<< ' ' << peak_shapes_[i].right_endpoint->getMZ()
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
        for (unsigned int i = 0; i < peak_shapes_.size(); ++i)
        {
					// put it out only if the peak was not deconvoluted
					if(find(peaks_to_skip.begin(),peaks_to_skip.end(),i) == peaks_to_skip.end() )
						{
							

							OutputPeakType picked_peak;
							
							picked_peak.setIntensity(peak_shapes_[i].height);
							picked_peak.setMZ(peak_shapes_[i].mz_position);
							
							
							
							fillPeak(peak_shapes_[i],picked_peak);
							picked_peak_container.push_back(picked_peak);
						}
        }

      } // if (peak_shapes_.size() > 0)
    }

    /** @brief Applies the peak picking algorithm to a raw data point container.
        
      Picks the peaks in the input container (e.g. of type MSSpectrum<RawDataPoint1D >) 
      and writes the resulting peaks to the picked_peak_container (e.g. MSSpectrum<PickedPeak1D>).

        The ms_level should be one if the spectrum is a normal mass spectrum, or two if it is a tandem mass spectrum.
        
      @note This method assumes that the input_peak_container contains data points of type 
           RawDataPoint1D or any other class derived from DRawDataPoint. 
              
           The resulting peaks in the picked_peak_container (e.g. of type MSSpectrum<PickedPeak1D>)
           can be of type RawDataPoint1D or any other class derived from DRawDataPoint. 
           We recommend to use the PickedPeak1Dbecause it stores important information gained during
           the peak picking algorithm.
            
           If you use MSSpectrum you have to set the SpectrumSettings by your own.
    */
    template <typename InputPeakContainer, typename OutputPeakContainer >
    void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container, int ms_level = 1)
    {
      pick(input_peak_container.begin(), input_peak_container.end(), picked_peaks_container, ms_level);
    }


    /** @brief Picks the peaks in a range of MSSpectren.
        
      Picks the peaks successive in every scan in the intervall [first,last).
      The detected peaks are stored in a MSExperiment.
              
      @note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type RawDataPoint1D 
              or any other derived class of DRawDataPoint.
            For the resulting peaks we recommend to use the PickedPeak1Dbecause it stores important information gained during
            the peak picking algorithm.  

        @note You have to copy the ExperimentalSettings of the raw data by your own.  
    */
    template <typename InputSpectrumIterator, typename OutputPeakType >
    void pickExperiment(InputSpectrumIterator first,
                        InputSpectrumIterator last,
                        MSExperiment<OutputPeakType>& ms_exp_peaks)
    {
      unsigned int n = distance(first,last);
      startProgress(0,n,"pick peaks in mzData file");
      // pick peaks on each scan
      for (unsigned int i = 0; i < n; ++i)
      {
				setProgress(i);
        MSSpectrum< OutputPeakType > spectrum;
        InputSpectrumIterator input_it(first+i);
#ifdef DEBUG_PEAK_PICKING
        std::cout << "PeakPicker: Picking Scan " << input_it->getRT()<< std::endl;
#endif
        StopWatch timer;
        timer.start();

        // pick the peaks in scan i
        pick(*input_it,spectrum,input_it->getMSLevel());
        setProgress(i);
        timer.stop();
#ifdef DEBUG_PEAK_PICKING
        std::cout << "Picking took " << timer.getClockTime()  << std::endl;
        std::cout << "Found " << spectrum.size() << " peaks. " << std::endl;
#endif

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

    /** @brief Picks the peaks in a range of MSSpectren (and output data structure MSExperimentExtern).
          
        Picks the peaks successive in every scan in the intervall [first,last).
        The detected peaks are stored in a MSExperiment.
                
        @note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type RawDataPoint1D 
                 or any other derived class of DRawDataPoint.
              For the resulting peaks we recommend to use the PickedPeak1Dbecause it stores important information gained during
              the peak picking algorithm.  

           @note You have to copy the ExperimentalSettings of the raw data on your own.   
       */
    template <typename InputSpectrumIterator, typename OutputPeakType >
    void pickExperiment(InputSpectrumIterator first,
                        InputSpectrumIterator last,
                        MSExperimentExtern<OutputPeakType>& ms_exp_peaks)
    {
      unsigned int n = distance(first,last);
      startProgress(0,n,"pick peaks in mzData file");
      // pick peaks on each scan
      for (unsigned int i = 0; i < n; ++i)
      {
        MSSpectrum< OutputPeakType > spectrum;
        InputSpectrumIterator input_it(first+i);

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
          spectrum.getPrecursorPeak() = input_it->getPrecursorPeak();
          spectrum.setRT(input_it->getRT());
          spectrum.setMSLevel(input_it->getMSLevel());
          spectrum.getName() = input_it->getName();

          ms_exp_peaks.push_back(spectrum);
        }
      }
      endProgress();
    }

    /** @brief Picks the peaks in a MSExperiment.
        
      Picks the peaks on every scan in the MSExperiment.
      The detected peaks are stored in a MSExperiment.
              
      @note The input peaks should be of type RawDataPoint1D or any other derived class of DRawDataPoint.
            For the resulting peaks we recommend to use the PickedPeak1Dbecause it stores important information gained during
            the peak picking algorithm.   
    */
    template <typename InputPeakType, typename OutputPeakType >
    void pickExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_peaks)
    {
      // copy the experimental settings
      static_cast<ExperimentalSettings&>(ms_exp_peaks) = ms_exp_raw;

      pickExperiment(ms_exp_raw.begin(),ms_exp_raw.end(),ms_exp_peaks);
    }

    /** @brief Picks the peaks in a MSExperimentExtern.
          
        Picks the peaks on every scan in the MSExperiment.
        The detected peaks are stored in a MSExperiment.
                
        @note The input peaks should be of type RawDataPoint1D or any other derived class of DRawDataPoint.
              For the resulting peaks we recommend to use the PickedPeak1Dbecause it stores important information gained during
              the peak picking algorithm.   
       */
    template <typename InputPeakType, typename OutputPeakType >
    void pickExperiment(const MSExperimentExtern< InputPeakType >& ms_exp_raw, MSExperimentExtern<OutputPeakType>& ms_exp_peaks)
    {
      for (unsigned int i=0; i<ms_exp_raw.size();++i)
      {
        MSSpectrum< OutputPeakType > out_spec;

#ifdef DEBUG_PEAK_PICKING
        std::cout << "Picking scan " << i << std::endl;
        std::cout << "Size of input: " << ms_exp_raw[i].size() << std::endl;
#endif
        StopWatch watch;
        watch.start();
        // pick the peaks in scan i
        pick(ms_exp_raw[i],out_spec,ms_exp_raw[i].getMSLevel());
        watch.stop();

#ifdef DEBUG_PEAK_PICKING
        std::cout << "Picking this scan took " << watch.getClockTime() << " seconds. " << std::endl;
#endif
        // copy spectrum settings
        out_spec.setType(SpectrumSettings::PEAKS);

        // copy the spectrum information
        out_spec.getPrecursorPeak() = ms_exp_raw[i].getPrecursorPeak();
        out_spec.setRT(ms_exp_raw[i].getRT());
        out_spec.setMSLevel(ms_exp_raw[i].getMSLevel());
        out_spec.getName() = ms_exp_raw[i].getName();

        ms_exp_peaks.push_back(out_spec);
      }

      //pickExperiment(ms_exp_raw.begin(),ms_exp_raw.end(),ms_exp_peaks);
    }


    /// This function fills the members of a picked peak of type OutputPeakType.
    template <typename OutputPeakType>
    void fillPeak(const PeakShape& /* peak_shape */, OutputPeakType& /* picked_peak */)
    {}

   

  protected:
    /// Container the determined peak shapes
    std::vector<PeakShape> peak_shapes_;

    /// The continuous wavelet "transformer"
    ContinuousWaveletTransformNumIntegration wt_;

    /// The continuous wavelet "transformer" for the deconvolution
    ContinuousWaveletTransformNumIntegration wtDC_;

    /// The search radius for the determination of a peak's maximum position
    unsigned int radius_;

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


    /** @brief Class for the internal peak representation
        
        A regularData-Object which contents some additional useful informations
        for analysing peaks and their properties
    */
    class PeakArea_
    {
      typedef std::vector<RawDataPointType>::iterator RawDataPointIterator;

    public:
      PeakArea_() : left(0), max(0), right(0), left_behind_centroid(0)
      {}

      /** @brief Iterator defining a raw data peak.
         
         The left and right iterators delimit a range in the raw data which represents a raw peak.
         They define the raw peak endpoints. Max points to the raw data point in [left, right] with the highest intensity, the 
         maximum of the raw peak. 
         
         Left_behind_centroid points to the raw data point next to the estimates centroid position.
      */
      RawDataPointIterator left;
      RawDataPointIterator max;
      RawDataPointIterator right;
      RawDataPointIterator left_behind_centroid;
      /// The estimated centroid position.
      DPosition<1> centroid_position;
    };

    /// Computes the peak's left and right area
    void getPeakArea_(const PeakArea_& area, double &area_left, double &area_right);

    /// Returns the best fitting peakshape
    PeakShape fitPeakShape_(const PeakArea_& area, bool enable_centroid_fit);

    /** @brief Returns the squared pearson coefficient.

        Computes the correlation of the peak and the original data given by the peak enpoints area.left and area.right.
        If the value is near 1, the fitted peakshape and the raw data are expected to be very similar. 
    */
    double correlate_(const PeakShape& peak, const PeakArea_& area, int direction=0) const;


    /** @brief Finds the next maximum position in the wavelet transform wt.
        
        If the maximum is greater than peak_bound_cwt we search for the corresponding maximum in the raw data interval [first,last)
      given a predefined search radius radius. Only peaks with intensities greater than peak_bound_ 
      are relevant. If no peak is detected the method return false.
      For direction=1, the method runs from first to last given direction=-1 it runs the other way around.
    */
    bool getMaxPosition_(RawDataPointIterator first, RawDataPointIterator last, const ContinuousWaveletTransform& wt, PeakArea_& area, int distance_from_scan_border, int ms_level, int direction=1);


    /** @brief Determines a peaks's endpoints.
      
      The algorithm does the following:
          - let x_m be the position of the maximum in the data and let (x_l, x_r) be
            the left and right neighbours
      
      
        (1) starting from x_l', walk left until one of the following happens
               - the new point is lower than the original bound
                    => we found our left endpoint
      
               - the new point is larger than the last, but the point left from
                 the new point is smaller. In that case, we either ran into another
                 peak, or we encounter some noise. Therefore we now look in the cwt
                 at the position corresponding to this value. If the cwt here is
                 monotonous, we consider the point as noise and continue further to the
                 left. Otherwise, we probably found the beginning of a new peak and
                 therefore stop here.
      
        (2) analogous procedure to the right of x_r
    */
    bool getPeakEndPoints_(RawDataPointIterator first, RawDataPointIterator last,  PeakArea_ &area, int distance_from_scan_border, int& peak_left_index, int& peak_right_index);


    /** @brief Estimates a peak's centroid position.

       Computes the centroid position of the peak using all raw data points which are greater than 
       60% of the most intensive raw data point.
    */
    void getPeakCentroid_(PeakArea_& area);



    /// Computes the value of a theroretical lorentz peak at position x
    double lorentz_(double height, double lambda, double pos, double x);

    /** @brief Computes the threshold for the peak height in the wavelet transform and initializes the wavelet transform.

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
		
		/** @brief Separates overlapping peaks.

		    It determines the number of peaks lying underneath the initial peak using the cwt with different scales.
				Then a nonlinear optimzation procedure is applied to optimize the peak parameters.

		*/
    bool deconvolutePeak_(PeakShape& shape);

		/// Determines the number of peaks in the given mass range using the cwt
    int getNumberOfPeaks_(RawDataPointIterator& first,RawDataPointIterator& last, std::vector<double>& peak_values,
													int direction,DoubleReal resolution, ContinuousWaveletTransformNumIntegration& wt);

		/// Estimate the charge state of the peaks
    int determineChargeState_(std::vector<double>& peak_values);

		/// Add a peak
    void addPeak_(std::vector<PeakShape>& peaks_DC,PeakArea_& area,double left_width,double right_width);
		//@}
  }
  ; // end PeakPickerCWT


  /// Fills the members of a PickedPeak1D given an PeakShape
  template <>
  void PeakPickerCWT::fillPeak< PickedPeak1D >(const PeakShape& peak_shape, PickedPeak1D& picked_peak);


}// namespace OpenMS

#endif
