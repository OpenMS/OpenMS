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
// $Authors: Eva Lange, Alexandra Zerck $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

//#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#define DEBUG_DECONV
namespace OpenMS
{
  /**
		 @brief This class implements a peak picking algorithm using wavelet techniques
			
		 The algorithm is descripted in detail in Lange et al. (2006) Proc. PSB-06.
			
		 This peak picking algorithm uses the continuous wavelet transform of a raw data signal to detect mass peaks.
		 Afterwards a given asymmetric peak function is fitted to the raw data and important peak parameters (e.g. fwhm)
		 are extracted.
		 In an optional step these parameters can be optimized using a non-linear opimization method.
			
		 The peak parameters are stored in the meta data arrays of the spectra (see MSSpectrum) in this order:
		 - rValue
		 - area
		 - fwhm
		 - leftWidth
		 - rightWidth
		 - peakShape
		 - SignalToNoise
		 .

		 @note The peaks must be sorted according to ascending m/z!
		 
		 @htmlinclude OpenMS_PeakPickerCWT.parameters
	  
		 @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerCWT
		: public DefaultParamHandler,
			public ProgressLogger
  {
	public:
    /// Raw data iterator type
    typedef MSSpectrum<>::iterator PeakIterator;
    /// Const raw data iterator type
    typedef MSSpectrum<>::const_iterator ConstPeakIterator;

    /// Constructor
    PeakPickerCWT();

    /// Destructor
    virtual ~PeakPickerCWT();

    /** 
				@brief Applies the peak picking algorithm to a single spectrum.
	        
				Picks the peaks in the input spectrum and writes the resulting peaks to the output container.
    */
    void pick(const MSSpectrum<>& input, MSSpectrum<>& output);

    /** 
				@brief Picks the peaks in an MSExperiment.
			
				Picks the peaks successive in every scan in the spectrum range. The detected peaks are stored in the output MSExperiment.
    */
    void pickExperiment(const MSExperiment<>& input, MSExperiment<>& output);

		/**
			 @brief Estimates average peak width that can then be used for peak picking.

			 The spectra with the highest TICs are used to estimate an average peak width that
			 can be used as the peak_width parameter for picking the complete data set.
			 Typically, the number of peaks increases with decreasing peak width until a plateau
			 is reached. The beginning of this plateau is our estimate for the peak width.
			 This estimate is averaged over several spectra.

		*/
		DoubleReal estimatePeakWidth(const MSExperiment<>& input);
	protected:

    /// Threshold for the peak height in the MS 1 level
    float peak_bound_;

    /// Threshold for the peak height in the MS 2 level
    float peak_bound_ms2_level_;

    /// Signal to noise threshold
    float signal_to_noise_;

    /// The minimal full width at half maximum
    float fwhm_bound_;

    /// The search radius for the determination of a peak's maximum position
    UInt radius_;

    /// The dilation of the wavelet
    float scale_;

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
    struct OPENMS_DLLAPI PeakArea_
    {
      typedef MSSpectrum<>::iterator PeakIterator;

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
    double correlate_(const PeakShape& peak, const PeakArea_& area, Int direction=0) const;


    /** 
				@brief Finds the next maximum position in the wavelet transform wt.
	        
				If the maximum is greater than peak_bound_cwt we search for the corresponding maximum in the raw data interval [first,last)
				given a predefined search radius radius. Only peaks with intensities greater than peak_bound_ 
				are relevant. If no peak is detected the method return false.
				For direction=1, the method runs from first to last given direction=-1 it runs the other way around.
    */
    bool getMaxPosition_(PeakIterator first, PeakIterator last, const ContinuousWaveletTransform& wt, PeakArea_& area, Int distance_from_scan_border, Int ms_level, DoubleReal peak_bound_cwt,DoubleReal peak_bound_ms2_level_cwt,Int direction=1);


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
    bool getPeakEndPoints_(PeakIterator first, PeakIterator last,  PeakArea_ &area, Int distance_from_scan_border, Int& peak_left_index, Int& peak_right_index,ContinuousWaveletTransformNumIntegration& wt);


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
    void initializeWT_(ContinuousWaveletTransformNumIntegration& wt,DoubleReal& peak_bound_cwt,DoubleReal& peak_bound_ms2_level_cwt);

		/** @name Methods needed for separation of overlapping peaks
		 */
    //@{
		
		/** 
				@brief Separates overlapping peaks.
	
				It determines the number of peaks lying underneath the initial peak using the cwt with different scales.
				Then a nonlinear optimzation procedure is applied to optimize the peak parameters.
		*/
    bool deconvolutePeak_(PeakShape& shape,std::vector<PeakShape>& peak_shapes,DoubleReal peak_bound_cwt);

		/// Determines the number of peaks in the given mass range using the cwt
    Int getNumberOfPeaks_(ConstPeakIterator first,ConstPeakIterator last, std::vector<double>& peak_values,
													Int direction,DoubleReal resolution, ContinuousWaveletTransformNumIntegration& wt,DoubleReal peak_bound_cwt);

		/// Estimate the charge state of the peaks
    Int determineChargeState_(std::vector<double>& peak_values);

		/// Add a peak
    void addPeak_(std::vector<PeakShape>& peaks_DC,PeakArea_& area,double left_width,double right_width,OptimizePeakDeconvolution::Data& data);
		//@}
  }
		; // end PeakPickerCWT


}// namespace OpenMS

#endif
