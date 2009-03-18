// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Rene Hussong$
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETTRANSFORM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETTRANSFORM_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgMzFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletConstants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/CoupledMarrWavelet.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics_double.h>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>

#ifdef OPENMS_HAS_CUDA
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletCudaKernel.h>
#endif


// we are not yet sure if we really want to drag in cutil.h and the CUDA_SAFE_CALL definitions...
#ifndef CUDA_SAFE_CALL
#define CUDA_SAFE_CALL(call) call;
#endif

namespace OpenMS
{

	class cudaHelp
	{
		public:
			float getMZ()
			{ return mz; }
			float getIntensity()
			{ return intens; }
		
			float mz;
			float intens;
			float score;
	};	


	/** @brief A class implementing the isotope wavelet transform.
 		* If you just want to find features using the isotope wavelet, take a look at the IsotopeWaveletFF class. Usually, you only
 		* have to consider the class at hand if you plan to change the basic implementation of the transform.
 	*/
	template <typename PeakType>
	class IsotopeWaveletTransform
	{
		public:

			/** @brief Internally used data structure. */
			struct BoxElement
			{
				DoubleReal mz;
				UInt c; //<Note, this is not the charge (it is charge-1!!!)
				DoubleReal score;
				DoubleReal intens;
				DoubleReal max_intens;
				DoubleReal RT; //<The elution time (not the scan index)
				UInt RT_index;
				UInt MZ_begin; //<Index
				UInt MZ_end; //<Index
			};

			typedef std::multimap<UInt, BoxElement> Box; ///<Key: RT index, value: BoxElement


			/** @brief Constructor.
 				*
 				* @param min_mz The smallest m/z value occurring in your map.
 				* @param max_mz The largest m/z value occurring in your map.
 				* @param max_charge The highest charge state you would like to consider. */
			IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const DoubleReal sigma=0.2, const UInt max_scan_size=0) ;

			/** @brief Destructor. */
			virtual ~IsotopeWaveletTransform () ;
	

			virtual bool estimateCMarrWidth (const MSSpectrum<PeakType>& scan);
	
			/** @brief Computes the discrete-time continuous wavelet transform simultaneously for several charges.
 				*
 				* The function computes the isotope wavelet transformed versions of @p scan.
 				* The transform is determined for several charge states (up to charge @p max_charge) at the same time.
 				* Hence, the user has to ensure that the size of @p transforms equals to @p max_charge and that each spectrum in
 				* @p transforms has the same length as @p scan.
 				*
 				* @param scan The MS scan you wish to transform.
 				* @param transforms A vector (with indices running from 0 to @p max_charge-1) of MS spectra (each of the size of @p scan).
 				* The code will NOT check the allocated memory (the sizes) for @p transforms and its entries.
 				* @param max_charge The maximal charge state that is considered. */
			virtual void getTransforms (const MSSpectrum<PeakType>& scan,
				std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge);

	
			#ifdef OPENMS_HAS_CUDA
				virtual void getCudaTransforms (MSSpectrum<PeakType> &c_trans, const UInt c);
				
				virtual int initializeCudaScan (const MSSpectrum<PeakType>& scan, const UInt cudaDevice); 
				
				virtual void finalizeCudaScan ();
			
				virtual void identifyCudaCharges (const MSSpectrum<PeakType>& candidates,
					const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt c, const DoubleReal ampl_cutoff, const bool use_cmarr=false);
	
				virtual int cudaSort (MSSpectrum<PeakType>& sorted);
			#endif	


			virtual void getCMarrTransforms (const MSSpectrum<PeakType>& scan, 
				std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge);


			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
 				* process is stored internally. Important: Before calling this function, apply updateRanges() to the original map.
 				*
 				* @param candidates A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
 				* charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function @see getTransforms.)
 				* @param ref The reference scan (the untransformed raw data) corresponding to @p candidates.
				* @param scan_index The index of the scan (w.r.t. to some map) currently under consideration.
 				* @param ampl_cutoff The thresholding parameter. This parameter is the only (and hence a really important)
 				* parameter of the isotope wavelet transform. On the basis of @p ampl_cutoff the program tries to distinguish between
 				* noise and signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
 				* line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
 				* get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
 				* makes sense to start @p ampl_cutoff=0 or even @p ampl_cutoff=-1,
 				* indicating no thresholding at all. Note that also ampl_cutoff=0 triggers (a moderate) thresholding based on the
 				* average intensity in the wavelet transform. */
			virtual void identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates,
				const MSSpectrum<PeakType>& ref, const UInt scan_index, const DoubleReal ampl_cutoff, const bool use_cmarr=false) ;

			/** @brief A function keeping track of currently open and closed sweep line boxes.
 				* This function is used by the isotope wavelet feature finder and must be called for each processed scan.
 				* @param map The original map containing the data set to be analyzed.
				* @param scan_index The index of the scan currently under consideration w.r.t. its MS map.
 				* This information is necessary to sweep across the map after each scan has been evaluated.
 				* @param RT_interleave See the IsotopeWaveletFF class.
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class. */
			void updateBoxStates (const MSExperiment<PeakType>& map, const Size scan_index, const UInt RT_interleave,
				const UInt RT_votes_cutoff) ;


			/** @brief Filters the candidates further more and maps the internally used data structures to the OpenMS framework.
 				* @param map The original map containing the data set to be analyzed.
 				* @param max_charge The maximal charge state under consideration.
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class.*/
			FeatureMap<Feature> mapSeeds2Features (const MSExperiment<PeakType>& map, const UInt max_charge, const UInt RT_votes_cutoff) ;


			/** @brief Returns the closed boxes. */
			virtual std::multimap<DoubleReal, Box> getClosedBoxes ()
				{ return (closed_boxes_); };


			inline DoubleReal getLinearInterpolation (const typename MSSpectrum<PeakType>::const_iterator& left_iter, const DoubleReal mz_pos, const typename MSSpectrum<PeakType>::const_iterator& right_iter)
			{
				return (left_iter->getIntensity() + (right_iter->getIntensity() - left_iter->getIntensity())/(right_iter->getMZ() - left_iter->getMZ()) * (mz_pos-left_iter->getMZ())); 
			};


			/** @brief Estimates the number of peaks of an isotopic pattern at mass @p mass and charge state @p z.
 				*
 				* @param mass The mass.
 				* @param z The charge. */
			inline UInt getPeakCutOff (const DoubleReal mass, const UInt z)
			{
				return ((UInt) ceil(peak_cutoff_intercept_+peak_cutoff_slope_*mass*z));
			};

			inline UInt getNumPeakCutOff (const DoubleReal mass, const UInt z)
			{ 
				DoubleReal mz=mass*z;
				return ((int) ceil(Constants::CUTOFF_FIT99_0+Constants::CUTOFF_FIT99_1*mz+Constants::CUTOFF_FIT99_2*mz*mz-Constants::IW_QUARTER_NEUTRON_MASS)); 
			};	

			inline UInt getMzPeakCutOffAtMonoPos (const DoubleReal mass, const UInt z)
			{ 
				DoubleReal mz=mass*z;
				//in principle we need here (-0.25+0.25)*IW_NEUTRON_MASS to include the peaks as a whole (draw a picture to see this)
				return ((int) ceil(Constants::CUTOFF_FIT99_0+Constants::CUTOFF_FIT99_1*mz+Constants::CUTOFF_FIT99_2*mz*mz)); 
			};	

			/*std::vector<DoubleReal> getErrorProneScans () const
			{
				return (error_prone_scans_);
			};	

			void clearErrorProneScans () 
			{
				error_prone_scans_.clear();
			};*/  

			inline DoubleReal getSigma () const
			{
				return (sigma_);
			};

			inline void setSigma (const DoubleReal sigma)
			{
				sigma_ = sigma;
			};


		protected:


			/** @brief Default Constructor.
 				* @note Provided just for inheritance reasons. You should always use the other constructor. */
			IsotopeWaveletTransform () ;


			void estimatePeakCutOffs_ (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge) ;


			/** @brief Samples the wavelet at discrete time points, s.t. they match automatically the m/z positions provided
 				* in @p scan. The discrete values of psi are stored in the member variable psi_.
 				*
 				*	@param scan Provides the sampling positions.
 				*	@param wavelet_length The number of sampling points for the wavelet.
 				*	@param mz_index The start index of @p scan for which the wavelet should be adapted.
 				* @param offset The offset the wavelet function needs to be aligned with a signal point.
 				* @param charge The charge (not the index c!) the wavelet function should adapt (corresponds to z in the paper).
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern. */
			inline void sampleTheIsotopeWavelet (const MSSpectrum<PeakType>& scan, const UInt wavelet_length, 
				const UInt mz_index, const DoubleReal offset, const UInt charge);
			
			inline void sampleTheCMarrWavelet_ (const MSSpectrum<PeakType>& scan, const Int wavelet_length, const Int mz_index, const UInt charge);	
			

			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern.
				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined.
 				* @param intens The intensity of the transform at @p seed_mz.
 				* @param ampl_cutoff The threshold. */
			virtual DoubleReal scoreThis_ (const MSSpectrum<PeakType>& candidate, const UInt peak_cutoff,
				const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) ;
			virtual DoubleReal scoreThisAlternative_ (const MSSpectrum<PeakType>& candidate, const UInt peak_cutoff, 
				const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) ;

			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score 
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern.
				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined. 
 				* @param intens The intensity of the transform at @p seed_mz.
 				* @param ampl_cutoff The threshold. */
			//virtual DoubleReal cudaScoreThis_ (const MSSpectrum<PeakType>& candidate, const UInt peak_cutoff, 
				//const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) ;


			/** @brief A ugly but necessary function to handle "off-by-1-Dalton predictions" due to idiosyncrasies of the data set
 				* (in comparison to the averagine model)
 				* @param candidate The wavelet transformed spectrum containing the candidate.
 				* @param ref The original spectrum containing the candidate.
 				* @param seed_mz The m/z position of the candidate pattern.
 				* @param c The predicted charge state of the candidate.
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). */
			virtual bool checkPosition_ (const MSSpectrum<PeakType>& candidate, const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, 
				const UInt c, const UInt scan_index, const bool use_cmarr) ;
			

			virtual std::pair<DoubleReal, DoubleReal> correctMZ_ (const MSSpectrum<PeakType>& ref, const DoubleReal c_mz, const DoubleReal c) ;


			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getAvIntens_ (const MSSpectrum<PeakType>& scan) ;

			/** @brief Computes the standard deviation (neglecting negative values) of the intensity of @p scan. */
			inline DoubleReal getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) ;

			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getPointerAvIntens_ (const MSSpectrum<PeakType*>& scan) ; 		
			
			/** @brief Computes the standard deviation (neglecting negative values) of the intensity of @p scan. */
			inline DoubleReal getPointerSdIntens_ (const MSSpectrum<PeakType*>& scan, const DoubleReal mean) ;


			/** @brief A wrapper function to the GSL interpolation routine. */
			DoubleReal getCubicInterpolatedValue_ (const std::vector<DoubleReal>& x, const DoubleReal xi, const std::vector<DoubleReal>& y) ;

			/** @brief A function to map m/z values to m/z indices. In particular useful if you know already the
 				* approximate position of the corresponding entry which can be indicated by @p start. */
			inline std::pair<Int, Int> getNearBys_ (const MSSpectrum<PeakType>& signal, const DoubleReal mz,
				const UInt start=0) const ;

			/** @brief Inserts a potential isotopic pattern into an open box or - if no such box exists - creates a new one.
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the post-processing (sweep lining).
 				* @param charge The estimated charge state of the pattern.
 				* @param score The pattern's score.
 				* @param intens The intensity at the monoisotopic peak.
 				* @param rt The retention time of the scan (similar to @p scan, but here: no index, but the real value).
 				* @param MZ_begin The starting index of the pattern (m/z) w.r.t. the current scan.
 				* @param MZ_end The end index (w.r.t. the monoisotopic position!) of the pattern (m/z) w.r.t. the current scan. */
			virtual void push2Box_ (const DoubleReal mz, const UInt scan, UInt c, const DoubleReal score, 
				const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) ;

			/** @brief Essentially the same function as @see push2Box_.
 				* In contrast to @see push2Box this function stores its candidates only temporarily. In particular, this
 				* function is only used within a single scan transform. After the wavelet transform is computed on
 				* that scan, all candidates are pushed by this function and finally clustered together by @see clusterSeeds_.
 				* Afterwards, a final push by @see push2Box_ is performed storing the clustered candidates.
 				*
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the post-processing (sweep lining).
 				* @param charge The estimated charge state of the pattern.
 				* @param score The pattern's score.
 				* @param intens The intensity at the monoisotopic peak.
 				* @param rt The retention time of the scan (similar to @p scan, but here: no index, but the real value).
 				* @param MZ_begin The starting index of the pattern (m/z) w.r.t. the current scan.
 				* @param MZ_end The end index (w.r.t. the monoisotopic position!) of the pattern (m/z) w.r.t. the current scan.*/
			virtual void push2TmpBox_ (const DoubleReal mz, const UInt scan, UInt charge, const DoubleReal score,
				const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) ;


			/** @brief Computes the average MZ spacing of @p scan in the range @p start_index to @p end_index.
 				*
 				* @param scan	The scan we are interested in.
 				* @param start_index An optional starting position (index) w.r.t. @p scan.
 				* @param end_index An optional final position (index) w.r.t. @p scan.*/
			inline DoubleReal getAvMZSpacing_ (const MSSpectrum<PeakType>& scan);//, Int start_index=0, Int end_index=-1) ;


			/** @brief The trapezoid rule for integration.
 				* @param a first x coordinate.
 				* @param b second x coordinate.
 				* @param fa a's corresponding function value.
 				* @param fb b's corresponding function value. */
			inline DoubleReal chordTrapezoidRule_ (const DoubleReal a, const DoubleReal b, const DoubleReal fa, const DoubleReal fb)
			{
				return ((fb+fa)*0.5*(b-a));
			};

			/** @brief The trapezoid rule for integration.
 				*	@param x The x coordinates.
 				*	@param y The function values. */
			inline DoubleReal chordTrapezoidRule_ (const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y)
			{
				DoubleReal res=(x[1]-x[0])*(y[0]);
				for (Size i=0; i<x.size()-1; ++i)
				{
					res += (x[i+1]-x[i])*(y[i+1]+y[i]);
				};
				return (0.5*res);
			};


			/** @brief Clusters the seeds stored by push2TmpBox_.
 				* @param candidates A isotope wavelet transformed spectrum.
 				* @param ref The corresponding original spectrum (w.r.t. @p candidates).
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map).
 				* @param max_charge The maximal charge state we will consider. */
			void clusterSeeds_ (const std::vector<MSSpectrum<PeakType> >& candidates,
				const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt max_charge, const bool use_cmarr) ;

			void clusterCudaSeeds_ (const MSSpectrum<PeakType>& candidates, 
				const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt c, const bool use_cmarr) ;



			void extendBox_ (const MSExperiment<PeakType>& map, const Box box);

			inline DoubleReal peptideMassRule (const DoubleReal c_mass) const
			{
				DoubleReal correction_fac = c_mass / Constants::PEPTIDE_MASS_RULE_BOUND;
				DoubleReal old_frac_mass = c_mass - trunc(c_mass);
				DoubleReal new_mass = trunc(c_mass)* (1.+Constants::PEPTIDE_MASS_RULE_FACTOR)-trunc(correction_fac);
				DoubleReal new_frac_mass = new_mass - trunc(new_mass);
				
				if (new_frac_mass-old_frac_mass > 0.5)
				{
					new_mass -= 1.;
				}			

				if (new_frac_mass-old_frac_mass < -0.5)
				{
					new_mass += 1.;
				}			

				return (new_mass);
			};

			inline DoubleReal getPPMs (const DoubleReal mass_a, const DoubleReal mass_b) const
			{
				return (fabs(mass_a-mass_b)/(0.5*(mass_a+mass_b))*1e6);
			};


			//internally used data structures for the sweep line algorithm
			std::multimap<DoubleReal, Box> open_boxes_, closed_boxes_;	//DoubleReal = average m/z position
			std::vector<std::multimap<DoubleReal, Box> >* tmp_boxes_; //for each charge we need a separate container

			gsl_interp_accel* acc_;
			gsl_spline* spline_;
			DoubleReal av_MZ_spacing_, peak_cutoff_intercept_, peak_cutoff_slope_, sigma_;  
			std::vector<DoubleReal> c_mzs_, c_spacings_, psi_, prod_, xs_;
			//std::vector<DoubleReal> error_prone_scans_;
			std::vector<DoubleReal> interpol_xs_, interpol_ys_;

			Int num_elements_;
			UInt max_scan_size_, largest_array_size_, max_num_peaks_per_pattern_;
			std::vector<int> indices_;
			float *h_data_;		
			int *h_pos_;	
			
			UInt overall_size_, from_max_to_left_, from_max_to_right_, block_size_, data_length_, to_load_, to_compute_;
			MSSpectrum<PeakType> c_sorted_candidate_;
			DoubleReal min_spacing_, max_mz_cutoff_;
			std::vector<float> scores_, zeros_;	

			#ifdef OPENMS_HAS_CUDA
				void* cuda_device_intens_;
				void* cuda_device_pos_;
				void* cuda_device_trans_intens_;
				void* cuda_device_fwd2_;
				void* cuda_device_posindices_sorted_;
				void* cuda_device_trans_intens_sorted_;
				void* cuda_device_scores_;	
				std::vector<float> cuda_positions_, cuda_intensities_;
				dim3 dimGrid_, dimBlock_;
			#endif
	};

	bool myCudaComparator (const cudaHelp& a, const cudaHelp& b);

	template <typename PeakType>	
	bool intensityComparator (const PeakType& a, const PeakType& b)
	{
		return (a.getIntensity() > b.getIntensity());
	}		

	template <typename PeakType>	
	bool intensityAscendingComparator (const PeakType& a, const PeakType& b)
	{
		return (a.getIntensity() < b.getIntensity());
	}			

	template <typename PeakType>
	bool intensityPointerComparator (PeakType* a, PeakType* b)
	{
		return (a->getIntensity() > b->getIntensity());
	}			

	template <typename PeakType>	
	bool positionComparator (const PeakType& a, const PeakType& b)
	{
		return (a.getMZ() < b.getMZ());
	}		

	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform ()
	{
		acc_ = gsl_interp_accel_alloc ();
		spline_ = gsl_spline_alloc (gsl_interp_cspline, Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS); 
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box> > (1);
		av_MZ_spacing_=1;
		num_elements_ = 0;
		max_scan_size_ = 0;
		largest_array_size_ = 0;
		max_mz_cutoff_ = 3;
		max_num_peaks_per_pattern_ = 3;
		h_data_ = NULL;
		h_pos_ = NULL;
	}

	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, 
		const DoubleReal sigma, const UInt max_scan_size) 
	{
		acc_ = gsl_interp_accel_alloc ();
		spline_ = gsl_spline_alloc (gsl_interp_cspline, Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS); 
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box> > (max_charge);
		IsotopeWavelet::init (max_mz, max_charge);				
		CoupledMarrWavelet::init (max_mz, max_charge, sigma);
		av_MZ_spacing_=1;
		estimatePeakCutOffs_ (min_mz, max_mz, max_charge);
		max_mz_cutoff_ =  getMzPeakCutOffAtMonoPos(max_mz, max_charge);
		max_num_peaks_per_pattern_=  getNumPeakCutOff(max_mz, max_charge);
		
		#ifdef OPENMS_HAS_CUDA //if (max_scan_size > 0) //only important for the GPU
			max_scan_size_ = max_scan_size;
			largest_array_size_ =  pow(2, ceil(log(max_scan_size +  Constants::CUDA_EXTENDED_BLOCK_SIZE_MAX)/log(2.0)));
			//std::cout << "largest_array_size: " << max_scan_size << "\t" << largest_array_size_ << std::endl;
		
			cuda_positions_.reserve(largest_array_size_);
		  cuda_intensities_.reserve(largest_array_size_);	
			indices_.resize (largest_array_size_);
			for (UInt q=0; q<largest_array_size_; ++q)
				indices_[q] = q; 

			h_data_ = (float*) malloc (largest_array_size_*sizeof(float));		
			h_pos_ = (int*) malloc (largest_array_size_*sizeof(int)); 	
		#else
			largest_array_size_ = 0;		
			h_data_ = NULL;
			h_pos_ = NULL;
			Int size_estimate =  (Int)ceil(max_scan_size/(max_mz-min_mz));
			psi_.reserve (size_estimate*max_num_peaks_per_pattern_*Constants::IW_NEUTRON_MASS); //The wavelet
			prod_.reserve (size_estimate*max_num_peaks_per_pattern_*Constants::IW_NEUTRON_MASS); 
			xs_.reserve (size_estimate*max_num_peaks_per_pattern_*Constants::IW_NEUTRON_MASS);
			interpol_xs_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
			interpol_ys_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
		#endif
	}


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::~IsotopeWaveletTransform ()
	{
		gsl_interp_accel_free (acc_);
		gsl_spline_free (spline_);
		if (h_data_ != NULL) free (h_data_);
		if (h_pos_ != NULL) free (h_pos_);
		delete (tmp_boxes_);
	}

			
	template <typename PeakType>
 	bool IsotopeWaveletTransform<PeakType>::estimateCMarrWidth (const MSSpectrum<PeakType>& scan)
	{
		typename MSSpectrum<PeakType>::const_iterator max_iter, l_bound, r_bound;
		MSSpectrum<PeakType> help;
		help.assign (scan.begin(), scan.end());
		sort (help.begin(), help.end(), intensityComparator<PeakType>);	

		UInt trials=0; sigma_ = -1;	
		DoubleReal c_sigma=0, success=0;
		for (max_iter=help.begin(); max_iter != help.end(); ++max_iter, ++trials)
		{
			if (trials >= 10)
			{
				if (success < 5)
				{
					av_MZ_spacing_ = getAvMZSpacing_(scan);
					sigma_ = 3*av_MZ_spacing_;
					return (false);
				}
				else
				{
					sigma_ = c_sigma/success;
					return (true);
				};
			};
		
			InterpolationModel* model = NULL;
			Fitter1D* fitter = EmgMzFitter1D::create();
			Param params;
			DPeakArray<PeakType> range;	MSSpectrum<PeakType> s_range;
			l_bound = scan.MZBegin(max_iter->getMZ()-Constants::IW_QUARTER_NEUTRON_MASS);
			r_bound = scan.MZBegin(max_iter->getMZ()+Constants::IW_QUARTER_NEUTRON_MASS);
			range.assign (l_bound, r_bound); s_range.assign (l_bound, r_bound);
			av_MZ_spacing_ = getAvMZSpacing_(scan);
			params.setValue ("interpolation_step", av_MZ_spacing_);
			params.setValue ("max_iteration", 30);
			params.setValue ("tolerance_stdev_bounding_box", 0.);
			
			fitter->setParameters (params);
			if (fitter->fit1d (range, model) > 0.8)
			{
				//std::cout << "success" << std::endl;
				model = reinterpret_cast<EmgModel*> (model);
				params = model->getParameters();

				DoubleReal width = params.getValue ("emg:width");
				DoubleReal alpha = params.getValue ("emg:symmetry");
		
				//This due to Cai et al.: "Statistical Moment Analysis and Deconvolution of Overlapping Chromatographic Peaks". 
				//Chromatographia Vol. 31, No. 11/12, 1991.
				c_sigma += sqrt(width*width+alpha*alpha);
				//std::cout << max_iter->getMZ() << "\t" << sqrt(width*width+alpha*alpha) << std::endl;
				++success;
			}
			/*else
			{
				std::cout << "failed" << std::endl;
			};*/
		};

		return (false);
	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getTransforms (const MSSpectrum<PeakType>& scan,
		std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge) 
	{
		UInt scan_size = scan.size(), wavelet_length=0, old_length=0, peak_cutoff=0;
		DoubleReal mz_cutoff=0;
		av_MZ_spacing_ = getAvMZSpacing_(scan); 

		DoubleReal cum_spacing, c_spacing, //Helping variables
			max_w_monoi_intens=Constants::IW_QUARTER_NEUTRON_MASS, //The position of the monoisotopic peak within the coordinate sys. of the wavelet 
			sums=0, //Helping variables
			max_position_scan=0, //The position of the data point (within the scan) we want to align with
			align_offset, //Correction term; shifts the wavelet to get the desired alignment
			last;
		UInt c=0, k=0, j=0;
		DoubleReal c_charge; //DoubleReal, since we will oven divide by c_charge
		typename MSSpectrum<PeakType>::const_iterator wave_start, wave_end;

		//The upcoming variable is necessary to capture strange effects in special types of unequally spaced data sets.
		//Imagine some wholes in the m/z range (points the mass spectrometer did not sample). If they become larger than
		//Constants::IW_QUARTER_NEUTRON_MASS (considering the case of charge 1), several data points will share the same max_position, 
		//causing the upcoming code to crash since suddenly some m/z positions will occur twice. The interval of multiple
		//occurring points is stored by multiple_s and implicitly by i.
		std::vector<int> multiple_s (max_charge,-1);
		std::vector<DoubleReal> last_max_position_scan (max_charge, -1);
		bool repair=false;

		//Starting convolution
		for (Size i=0; i<scan_size; ++i)
		{
			//Now, let's sample the wavelets
			for (c=0; c<max_charge; ++c)
			{
				c_charge=c+1;
				cum_spacing=0;
				max_w_monoi_intens=Constants::IW_QUARTER_NEUTRON_MASS/c_charge; //This is the position of the monoisotopic peak (centered)

				//Align the maximum monoisotopic peak of the wavelet with some scan point. This is step is critical, since
				//otherwise we might - especially in the case of badly resolved data - miss patterns, since scan maxima and
				//wavelet maxima are "anticorrelated"
				j=0; last=0;
				while (cum_spacing < max_w_monoi_intens)
				{
					c_spacing = scan[(i+j+1)%scan_size].getMZ() - scan[(i+j)%scan_size].getMZ();
					last=cum_spacing;
					if (c_spacing < 0) //I.e. we are at the end of the scan
					{
						cum_spacing += av_MZ_spacing_;
					}
					else //The "normal" case
					{
						cum_spacing += c_spacing;
					};
					++j;
				};

				align_offset = max_w_monoi_intens - last; //I.e. we have to shift the wavelet by this amount to align the data
				--j;

				//The upcoming variable holds the position of the spectrum that is aligned with the monoisotopic
				//maximum of the wavelet. We do not add the overall correction term for the left shift at this point,
				//since we will get trouble by the Constants::IW_NEUTRON_MASS and the resulting numerical instabilities. 
				//We will add this correcting term at the end of the whole processing.
				if (i+j >= scan_size)
				{
					max_position_scan =  last_max_position_scan[c]+av_MZ_spacing_;
				}
				else
				{
					max_position_scan = scan[i+j].getMZ();
				};

				if (max_position_scan == last_max_position_scan[c]) //Uuups, multiple times the same m/z coordinate
				{
					if (multiple_s[c] < 0) //This is the first entry where this artifact occurred
					{
						multiple_s[c] = i-1;
					};
					//Notice that the problematic case of multiple_s being at the end of the spectrum (this might happen for
					//the overlapping part of the transform) can be ignored.
					//The special case if we are the boundary (exactly the last point in the spectrum).
					if (i == scan_size-1)
					{
						repair = true;
					};
				}
				else //Denotes the end of the multiple pos interval and triggers a repair.
				{
					if (multiple_s[c] >= 0)
					{
						repair=true; //We cannot do this now. Just after the transform at the actual point is completed.
					};
				};
	
				last_max_position_scan[c] = max_position_scan;
				cum_spacing = align_offset;

				//peak_cutoff = getPeakCutOff (scan[i].getMZ(), (UInt) c_charge);
				mz_cutoff = getMzPeakCutOffAtMonoPos(scan[i].getMZ(), (UInt) c_charge);//(peak_cutoff-1+0.75)*Constants::IW_NEUTRON_MASS;
				wave_start = scan.begin()+i;
				wave_end = scan.MZBegin(scan[i].getMZ()+mz_cutoff/c_charge);
					
				wavelet_length = distance (wave_start, wave_end);
				
				if (wavelet_length >= scan_size || wavelet_length <=1 || (scan[i+wavelet_length-1].getMZ() - scan[i].getMZ() > mz_cutoff+Constants::IW_NEUTRON_MASS/c_charge))
				{			
					sums=-1;
					/*if (error_prone_scans_.empty())
					{
						std::cout << "pushing to error_prone_scans_ a " << i << std::endl;
						std::cout << "wavelet length: " << wavelet_length << "\t scan_size: " << scan_size << std::endl; 
						error_prone_scans_.push_back(i);
					}
					else
					{
						if (*(--error_prone_scans_.end()) != i)
						{
							std::cout << "pushing to error_prone_scans_ b " << i << std::endl;
							error_prone_scans_.push_back(i);
						}
					};*/
				}
				else
				{
					if (wavelet_length != old_length)
					{
						psi_.resize (wavelet_length, 0);
						prod_.resize (wavelet_length, 0);
						xs_.resize (wavelet_length, 0);
						c_mzs_.resize (wavelet_length+1, 0);
						c_spacings_.resize (wavelet_length, 0);
						old_length = wavelet_length;
					};

					memset(&(psi_[0]), 0, sizeof(DoubleReal)*psi_.size());
					memset(&(prod_[0]), 0, sizeof(DoubleReal)*prod_.size());
					memset(&(xs_[0]), 0, sizeof(DoubleReal)*xs_.size());
					memset(&(c_mzs_[0]), 0, sizeof(DoubleReal)*c_mzs_.size());
					memset(&(c_spacings_[0]), 0, sizeof(DoubleReal)*c_spacings_.size());

					//Sampling the wavelet
					sampleTheIsotopeWavelet (scan, wavelet_length, i, cum_spacing, (UInt) c_charge);

					k=0; sums=0;
					for (UInt j=i; j<scan_size && k<wavelet_length; ++j, ++k)
					{
						prod_[k] = scan[j].getIntensity()*psi_[k];
						xs_[k] = scan[j].getMZ();
					};

					if (k< wavelet_length) // I.e. we have an overlapping wavelet
					{
						sums = 0;
						max_position_scan = transforms[c][i-1].getMZ()+av_MZ_spacing_;
					}
					else
					{
						sums = chordTrapezoidRule_ (xs_, prod_);
							sums += prod_[k];
						};*/

					};
				}
				//Store the current convolution result
				PeakType c_peak1 (transforms[c][i]);
				c_peak1.setIntensity(sums);
				c_peak1.setMZ(max_position_scan);
				transforms[c][i] = c_peak1;

				//std::cout << "Setting" << i << "\t\t" << transforms[c][i] << std::endl;

				if (repair)
				{	
					//std::cout << "repair" << std::endl;
	
					UInt noi2interpol = i - multiple_s[c]; //NOT +1

					//The special case if we are the boundary (exactly the last point in the spectrum)
					if (i == scan_size-1)
					{
						//We do not care about the intensities, since we will set them to zero anyway.
						//We would just like to avoid multiple positions to occur in the transform
						for (Size ii=0; ii<=noi2interpol; ++ii)
						//it must be "<=noi..." !!! not "<", since in this case we do not want to keep the last multiple
						//the same holds for "ii=0"
						{
							transforms[c][multiple_s[c]+ii].setMZ(transforms[c][multiple_s[c]-1].getMZ() + (ii+1)*av_MZ_spacing_);
						};

						last_max_position_scan[c] = max_position_scan; //Reset
						multiple_s[c]=-1; //Reset
						repair=false;
						continue;
					}

					PeakType& c_peak2 (transforms[c][multiple_s[c]]);
					DoubleReal x1 = c_peak2.getMZ();
					DoubleReal y1 = c_peak2.getIntensity();
					DoubleReal x2 = c_peak1.getMZ();
					DoubleReal y2 = c_peak1.getIntensity();
					DoubleReal dx = (x2-x1)/(DoubleReal)(noi2interpol);
					for (Size ii=1; ii<noi2interpol; ++ii) //ii=1, not 0, since we want to keep the first of the multiples
					{
						transforms[c][multiple_s[c]+ii].setMZ(c_peak2.getMZ()+ii*dx);
						transforms[c][multiple_s[c]+ii].setIntensity(y1 + (y2-y1)/(x2-x1)*(c_peak2.getMZ()+ii*dx-x1));
					};

					last_max_position_scan[c] = max_position_scan; //Reset
					multiple_s[c]=-1; //Reset
					repair=false;
				}					
			}
		}

		#ifdef DEBUG_FEATUREFINDER
			for (c=0; c<max_charge; ++c)
			{
				std::stringstream name; name << "trans_isowave_" << scan.getRT() << "_" << c+1 << ".dat\0"; 
				std::ofstream ofile (name.str().c_str());
				for (UInt i=0; i<transforms[c].size(); ++i)
					ofile << transforms[c][i].getMZ() << "\t" <<  transforms[c][i].getIntensity() << std::endl;
				ofile.close();
			};
		#endif

		return;
	}


	#ifdef OPENMS_HAS_CUDA
		template <typename PeakType>
		void IsotopeWaveletTransform<PeakType>::finalizeCudaScan ()
		{			
			(cudaFree(cuda_device_pos_));	
			(cudaFree(cuda_device_intens_));		
			(cudaFree(cuda_device_trans_intens_));
			(cudaFree(cuda_device_fwd2_));
			(cudaFree(cuda_device_trans_intens_sorted_));
			(cudaFree(cuda_device_posindices_sorted_));
			(cudaFree(cuda_device_scores_));
		}

		template <typename PeakType>
		int IsotopeWaveletTransform<PeakType>::initializeCudaScan (const MSSpectrum<PeakType>& scan, const UInt cudaDevice) 
		{
			(cudaSetDevice (cudaDevice));
			cudaDeviceProp prop;
			(cudaGetDeviceProperties(&prop, cudaDevice));
			
			UInt scan_size = scan.size();		

			std::vector<float> pre_positions (scan_size), pre_intensities (scan_size);
			for (UInt i=0; i<scan_size; ++i)
			{
				pre_positions[i] = scan[i].getMZ();
				pre_intensities[i] = scan[i].getIntensity();
			};
				
			/*(cudaMalloc(&cuda_device_pos_, scan_size*sizeof(float)));
			(cudaMemcpy(cuda_device_pos_, &(pre_positions[0]), scan_size*sizeof(float), cudaMemcpyHostToDevice));

			int max_wavelet_length = getMaxWaveletLength ((float*) cuda_device_pos_, scan_size, max_peaks_per_pattern_*Constants::IW_NEUTRON_MASS);
			//std::cout << "max_wavelet_length: " << max_wavelet_length << std::endl;
			
			cudaFree (cuda_device_pos_);*/


			//Determining the minimal spacing
			float c_spacing;
			min_spacing_=INT_MAX;
			for (UInt i=0; i<scan_size-1; ++i)
			{
				c_spacing = pre_positions[i+1]-pre_positions[i];
				if (c_spacing < min_spacing_)
				{
					min_spacing_ = c_spacing;
				};	
			};
			if (min_spacing_ == INT_MAX) //spectrum consists of a single data point
			{
				return (Constants::CUDA_INIT_FAIL);
			};
			

			//std::cout << "max_mz_cutoff: " << max_mz_cutoff_ << std::endl;
			//std::cout << "min_spacing_: " << min_spacing_ << "\t max_peaks_per_pattern: " << max_peaks_per_pattern_ << std::endl;	
			UInt wavelet_length = (UInt) ceil(max_mz_cutoff_/(min_spacing_));

			if (wavelet_length > scan_size || wavelet_length == 1) //==1, because of 'ceil'
			{
				return (Constants::CUDA_INIT_FAIL);
			}; 

			UInt max_index = (UInt) (Constants::IW_QUARTER_NEUTRON_MASS/min_spacing_);
			from_max_to_left_ = max_index;
			from_max_to_right_ = wavelet_length-1-from_max_to_left_;

			//Int problem_size = Constants::CUDA_BLOCK_SIZE_MAX - (wavelet_length-1);
		
			//Warum setzen wir die Problemgroesse nicht einfach auf 512?????????????????????????????????????????????????????
			//if (problem_size <= 0)
			//{
				Int problem_size = Constants::CUDA_BLOCK_SIZE_MAX;
			//};
			to_load_ = problem_size + from_max_to_left_ + from_max_to_right_;	
			
			data_length_ = pre_intensities.size(); // number of original data points
			//std::cout << "problem_size: " << problem_size << "\t" << (data_length_ % problem_size) << std::endl; 
			UInt missing_points = problem_size - (data_length_ % problem_size);
			overall_size_ = wavelet_length-1+data_length_+missing_points;
			//std::cout << "overall_size_: " << overall_size_ << "\t" << wavelet_length << "\t" << data_length_ << "\t" << missing_points << std::endl;

			//to_load_ += missing_points;

			num_elements_ = overall_size_;
			Int dev_num_elements = 1, tmp = overall_size_ >> 1;

			//Get power of 2 elements (necessary for the sorting algorithm)
			while (tmp) 
			{
					dev_num_elements <<= 1;
					tmp >>= 1;
			};

			if (num_elements_ > dev_num_elements)
			{
				dev_num_elements <<= 1;
			};

			if (dev_num_elements < Constants::CUDA_MIN_SORT_SIZE)
			{
				dev_num_elements = Constants::CUDA_MIN_SORT_SIZE;
			};

			overall_size_ = dev_num_elements;

			//std::cout << "overall_size: " << overall_size_ << "\t from_max_to_left_: " << from_max_to_left_ << "\t from_max_to_right_: " << from_max_to_right_ << "\t data_length_: " << data_length_ 
			//	<< "\t dev_num_elements-num_elements_: " << dev_num_elements-num_elements_ << std::endl;


			cuda_intensities_.resize (overall_size_, 0); cuda_positions_.resize (overall_size_, 0);
			//Pad the values to the left; the positions should not matter if the values are zero
			float first_pos = pre_positions[0];
			for (UInt i=0; i<from_max_to_left_; ++i)
			{
				cuda_positions_[i] = first_pos-(from_max_to_left_-i)*min_spacing_;
			};

			for (UInt i=0; i<data_length_; ++i)
			{
				cuda_positions_[from_max_to_left_+i] = pre_positions[i];
				cuda_intensities_[from_max_to_left_+i] = pre_intensities[i];
			};
			
			float last_pos = pre_positions[pre_positions.size()-1];
			for (UInt i=0; i<missing_points + from_max_to_right_ + dev_num_elements - num_elements_; ++i)
			{
				cuda_positions_[from_max_to_left_+data_length_+i] = last_pos + (i+1)*min_spacing_;
			};

			#ifdef DEBUG_FEATUREFINDER				
			std::stringstream name; name << "cuda_input_" << scan.getRT() << ".out\0"; 
				std::fstream outfile(name.str().c_str(), std::ios::out);
				for (size_t i=0; i<overall_size_; ++i)
					outfile << cuda_positions_[i] << " " << cuda_intensities_[i] << std::endl;
				outfile.close();
			#endif

			
			dimBlock_ = dim3 (Constants::CUDA_BLOCK_SIZE_MAX);
			to_compute_ = problem_size;

			//std::cout << "BlockSize: " << dimBlock_.x << "\t to load: " << to_load_  << "\t to compute: " << to_compute_ << std::endl; 

			//dimGrid_ = dim3 ((size_t)ceil(overall_size_/(float)problem_size)-1);
			//dimGrid_ = dim3 ((size_t)ceil((overall_size_-(dev_num_elements - num_elements_))/(float)problem_size));
			dimGrid_ = dim3 ((data_length_+missing_points)/problem_size);
			//std::cout << "griddim: " << dimGrid_.x << "\t blockdim: " << dimBlock_.x <<  std::endl;

			(cudaMalloc(&cuda_device_posindices_sorted_, overall_size_*sizeof(int)));

			(cudaMalloc(&cuda_device_pos_, overall_size_*sizeof(float)));
			(cudaMemcpy(cuda_device_pos_, &(cuda_positions_[0]), overall_size_*sizeof(float), cudaMemcpyHostToDevice));

			(cudaMalloc(&cuda_device_intens_, overall_size_*sizeof(float)));
			(cudaMemcpy(cuda_device_intens_, &(cuda_intensities_[0]), overall_size_*sizeof(float), cudaMemcpyHostToDevice));

			(cudaMalloc(&cuda_device_trans_intens_, overall_size_*sizeof(float)));
			
			(cudaMalloc(&cuda_device_fwd2_, overall_size_*sizeof(float)));

			(cudaMalloc(&cuda_device_trans_intens_sorted_, overall_size_*sizeof(float)));

			c_sorted_candidate_.resize (overall_size_);
			scores_.resize(data_length_);
			zeros_.resize(overall_size_);
			memset (&zeros_[0], 0., overall_size_*sizeof(float)); 

			(cudaMalloc(&cuda_device_scores_, overall_size_*sizeof(float)));
		
			return (Constants::CUDA_INIT_SUCCESS);
		}

		template <typename PeakType>
		void IsotopeWaveletTransform<PeakType>::getCudaTransforms (MSSpectrum<PeakType> &c_trans, const UInt c) 
		{
			std::vector<float> res (overall_size_, 0);
			(cudaMemcpy(cuda_device_trans_intens_, &zeros_[0], overall_size_*sizeof(float), cudaMemcpyHostToDevice));	
			(cudaMemcpy(cuda_device_fwd2_, &zeros_[0], overall_size_*sizeof(float), cudaMemcpyHostToDevice));	
			
			//std::cout << from_max_to_left_ << "\t" << from_max_to_right_ << "\t" << to_load_ << "\t" << to_compute_ << "\t" << data_length_ << std::endl;
			getExternalCudaTransforms (dimGrid_, dimBlock_, (float*)cuda_device_pos_, (float*)cuda_device_intens_, from_max_to_left_, from_max_to_right_, (float*)cuda_device_trans_intens_, 
				c+1, to_load_, to_compute_, peak_cutoff_intercept_, peak_cutoff_slope_, data_length_, (float*)cuda_device_fwd2_);
		
			//(cudaMemcpy(cuda_device_trans_intens_sorted_, cuda_device_trans_intens_, overall_size_*sizeof(float), cudaMemcpyDeviceToDevice));
			(cudaMemcpy(cuda_device_trans_intens_sorted_, cuda_device_fwd2_, overall_size_*sizeof(float), cudaMemcpyDeviceToDevice));
			(cudaThreadSynchronize());	

			(cudaMemcpy(&(res[0]), (float*)cuda_device_trans_intens_+from_max_to_left_, data_length_*sizeof(float), cudaMemcpyDeviceToHost));
			
			for (UInt i=0; i<data_length_; ++i)
			{
				c_trans[i].setIntensity (res[i]);
			};
		}


		template <typename PeakType>
		int IsotopeWaveletTransform<PeakType>::cudaSort (MSSpectrum<PeakType>& sorted)
		{

			(cudaMemcpy(cuda_device_posindices_sorted_, &indices_[0], overall_size_*sizeof(int), cudaMemcpyHostToDevice));
			Int gpu_index = sortOnDevice((float*)cuda_device_trans_intens_sorted_, (int*) cuda_device_posindices_sorted_, overall_size_, 0);
			(cudaMemcpy(h_data_, (float*)cuda_device_trans_intens_sorted_+gpu_index, sizeof(float) * (overall_size_-gpu_index), cudaMemcpyDeviceToHost));
			(cudaMemcpy(h_pos_, (int*)cuda_device_posindices_sorted_+gpu_index, sizeof(int) * (overall_size_-gpu_index), cudaMemcpyDeviceToHost));

			for (UInt i=0; i<(overall_size_-gpu_index); ++i)
			{
				sorted[i].setIntensity (h_data_[i]);
				sorted[i].setMZ (cuda_positions_[h_pos_[i]]);
			};
		
			return (gpu_index);
		}


		template <typename PeakType>
		void IsotopeWaveletTransform<PeakType>::identifyCudaCharges (const MSSpectrum<PeakType>& candidates,
			const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt c, const DoubleReal ampl_cutoff, const bool use_cmarr)
		{
			UInt index, MZ_start, MZ_end;
			typename MSSpectrum<PeakType>::iterator iter, bound_iter;
			typename MSSpectrum<PeakType>::const_iterator iter_start, iter_end, iter_p, iter2, seed_iter;
			DoubleReal mz_cutoff, seed_mz, c_av_intens=0, c_score=0, c_sd_intens=0, threshold=0, help_mz;
						
			Int gpu_index = cudaSort (c_sorted_candidate_), c_index;
			if (gpu_index < 0) //the transform produced non-exploitable data
			{
				return;
			};
	

			std::vector<UInt> processed (data_length_, 0);
			if (ampl_cutoff < 0)
			{
				threshold=0;
			}
			else
			{			
				c_av_intens = getAvIntens_ (candidates);
				c_sd_intens = getSdIntens_ (candidates, c_av_intens);
				threshold=ampl_cutoff*c_sd_intens + c_av_intens;
			};		
		
			Int num_of_scores = overall_size_-gpu_index;
		
			(cudaMemcpy(cuda_device_scores_, &zeros_[0], num_of_scores*sizeof(float), cudaMemcpyHostToDevice));
		
			//std::cout << "num_of_scores: " << num_of_scores << "\t" << c+1 << "\t" << gpu_index << std::endl;

			scoreOnDevice ((int*)cuda_device_posindices_sorted_, (float*)cuda_device_trans_intens_,  (float*)cuda_device_pos_, (float*)cuda_device_scores_, 
				c, num_of_scores, overall_size_, peak_cutoff_intercept_, peak_cutoff_slope_, max_num_peaks_per_pattern_);
			
			(cudaMemcpy(&scores_[0], cuda_device_scores_, num_of_scores*sizeof(float), cudaMemcpyDeviceToHost));

			std::vector<float>::iterator score_iter;
			for (c_index = overall_size_-gpu_index-1, score_iter = scores_.begin()+num_of_scores-1; c_index >= 0; --c_index, --score_iter)
			{				
				seed_mz = c_sorted_candidate_[c_index].getMZ();
				seed_iter = ref.MZBegin(seed_mz);
				index = distance(ref.begin(), seed_iter);

				if (seed_iter == ref.end() || processed[distance(ref.begin(), seed_iter)] || index <= 0)
				{
					continue;
				};
				
				mz_cutoff = getMzPeakCutOffAtMonoPos(seed_mz, (UInt) c+1);//;getPeakCutOff (seed_mz, c+1);
				//Mark the region as processed
				//Do not move this further down, since we have to mark this as processed in any case, 
				//even when score <=0; otherwise we would look around the maximum's position unless 
				//any significant point is found
				iter_start =ref.MZBegin(seed_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));
				//iter_end = ref.MZEnd(seed_mz+(peak_cutoff-1)-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));
				iter_end = ref.MZEnd(seed_mz+mz_cutoff/(c+1.));
			
				if (iter_end == ref.end())
				{
					--iter_end;
				};
			
				MZ_start = distance (ref.begin(), iter_start);
				MZ_end = distance (ref.begin(), iter_end);

				memset(&(processed[MZ_start]), 1, sizeof(UInt)*(MZ_end-MZ_start+1));

				c_score = *score_iter;
				if (c_score <=  c_sorted_candidate_[c_index].getIntensity()+ threshold)
				{
					continue;
				};


				//Push the seed into its corresponding box (or create a new one, if necessary)
				//Do ***NOT*** move this further down!
				push2TmpBox_ (seed_mz, scan_index, c, c_score, c_sorted_candidate_[c_index].getIntensity(), ref.getRT(), MZ_start, MZ_end);
					

				//Push neighboring peaks to compute finally a derivative over the isotope pattern envelope				
				help_mz = seed_mz - Constants::IW_NEUTRON_MASS/(c+1.);
				iter2 = candidates.MZBegin (help_mz);

				if (iter2 == candidates.end() || iter2 == candidates.begin())
				{
					continue;
				};

				if (fabs(iter2->getMZ()-seed_mz) > 0.5*Constants::IW_NEUTRON_MASS/(c+1.))
				//In the other case, we are too close to the peak, leading to incorrect derivatives.
				{
					if (iter2 != candidates.end())
					{
						//push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, iter2->getIntensity(), candidates.getRT(), MZ_start, MZ_end);
						push2TmpBox_ (help_mz, scan_index, c, 0, getLinearInterpolation(iter2-1, help_mz, iter2), candidates.getRT(), MZ_start, MZ_end);
					};
				};


				help_mz = seed_mz + Constants::IW_NEUTRON_MASS/(c+1.);
				iter2 = candidates.MZBegin (help_mz);

				if (iter2 == candidates.end() || iter2 == candidates.begin())
				{
					continue;
				};

				if (fabs(iter2->getMZ()-seed_mz) > 0.5*Constants::IW_NEUTRON_MASS/(c+1.))
				//In the other case, we are too close to the peak, leading to incorrect derivatives.
				{
					if (iter2 != candidates.end())
					{
						//push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, iter2->getIntensity(), candidates.getRT(), MZ_start, MZ_end);
						push2TmpBox_ (help_mz, scan_index, c, 0, getLinearInterpolation((iter2-1), help_mz, (iter2)), candidates.getRT(), MZ_start, MZ_end);
					};
				};
			};
			
			clusterCudaSeeds_(candidates, ref, scan_index, c, use_cmarr);
		}
	#endif



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getCMarrTransforms (const MSSpectrum<PeakType>& scan, 
		std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge) 
	{
		UInt scan_size = scan.size(), wavelet_length=0, new_length=0, old_length=0, peak_cutoff=0, c_mz_index=0;
		av_MZ_spacing_ = getAvMZSpacing_(scan);
		
		DoubleReal cum_spacing, //Helping variables
			max_w_monoi_intens=0, //The position of the monoisotopic peak within the coordinate sys. of the wavelet 
			sums=0, //Helping variables
			max_position_scan=0, //The position of the data point (within the scan) we want to align with
			align_offset; //Correction term; shifts the wavelet to get the desired alignment
		UInt c=0, k=0;
		DoubleReal c_charge; //DoubleReal, since we will oven divide by c_charge 
		typename MSSpectrum<PeakType>::const_iterator tmp_iter;

		//The upcoming variable is necessary to capture strange effects in special types of unequally spaced data sets.
		//Imagine some wholes in the m/z range (points the mass spectrometer did not sample). If they become larger than 
		//Constants::IW_QUARTER_NEUTRON_MASS (considering the case of charge 1), several data points will share the same max_position, 
		//causing the upcoming code to crash since suddenly some m/z positions will occur twice. The interval of multiple 
		//occurring points is stored by multiple_s and implicitly by i.
		std::vector<int> multiple_s (max_charge, -1);
		std::vector<DoubleReal> last_max_position_scan (max_charge, -1);
				
		//Determine the wavelet_length in advance since it is not depending on the m/z position
		//(as it is for the isotope wavelet)
		wavelet_length = (UInt) floor(Constants::MARR_WAVELET_CUTOFF*CoupledMarrWavelet::getSigma()/av_MZ_spacing_);
	
		if (wavelet_length > scan_size || wavelet_length <= 0)
		{		
			std::cerr << "Unable to sample an appropriate wavelet." << std::endl;
			std::cerr << "This is most probably no bug, since: cutoff=" << Constants::MARR_WAVELET_CUTOFF*CoupledMarrWavelet::getSigma()
				<< " and av m/z spacing=" << av_MZ_spacing_ << std::endl;
			std::cerr << "@Rene: replace this error message by some error handling stuff!" << std::endl;
			std::cerr << "Returning." << std::endl;
			return;
		};
					
		new_length = 2*wavelet_length+1;
		if (wavelet_length != old_length)
		{
			psi_.resize (new_length, 0);
			prod_.resize (new_length, 0);
			xs_.resize (new_length, 0);
			c_mzs_.resize (new_length+1, 0);
			c_spacings_.resize (new_length, 0);
			old_length = new_length;
		};

		
		//std::cout << "WARNING: THIS TRANSFORM ONLY WORKS FOR CHARGE 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;	
		for (UInt i=0; i<wavelet_length; ++i)
		{	
			transforms[0][i].setIntensity(0);
		};

		//Starting convolution
		for (UInt i=wavelet_length; i<scan_size; ++i)
		{
			//Now, let's sample the wavelets
			for (c=0; c<max_charge; ++c)
			{	
				c_charge=c+1;
				cum_spacing=0;				
				max_w_monoi_intens=i; //This is the position of the monoisotopic peak (centered)
				//Here, the position is equal to the i					
				align_offset = 0;
		

				//The upcoming variable holds the position of the spectrum that is aligned with the monoisotopic 
				//maximum of the wavelet. We do not add the overall correction term for the left shift at this point, 
				//since we will get trouble by the Constants::IW_NEUTRON_MASS and the resulting numerical instabilities. 
				//We will add this correcting term at the end of the whole processing.
				//std::cout << "Setting max postion to: " << scan[i].getMZ() << "\t" << j << std::endl;
				max_position_scan = scan[i].getMZ();
				cum_spacing = align_offset;
				
				sums=0;
				std::vector<std::pair<UInt, DoubleReal> > weights = IsotopeWavelet::getAveragine (scan[i].getMZ(), &peak_cutoff);
				std::vector<DoubleReal> helpv;
				for (UInt r=0; r<weights.size(); ++r)
				{		
					helpv.push_back(weights[r].second);
				};

				//std::stringstream name; name << "cmarr_wavelet_" << scan[i].getMZ() << "_" << c+1 << ".dat\0"; 
				//std::ofstream ofile (name.str().c_str());
				for (UInt r=0; r<peak_cutoff; ++r)
				{
					tmp_iter = scan.MZBegin(scan[i].getMZ()+r*Constants::IW_NEUTRON_MASS/c_charge);
					if (tmp_iter == scan.end())
					{
						break;
					};

					//c_mz_index = distance(scan.begin(), tmp_iter);
					c_mz_index = distance(scan.begin()+i, tmp_iter);

					//Sampling the wavelet 
					sampleTheCMarrWavelet_ (scan, wavelet_length, c_mz_index+i, (UInt) c_charge);

					k=0;
					//std::cout << wavelet_length << "\t" << i << "\t" << c_mz_index << "\t" << r << std::endl;
					//std::cout << "j: " << i+c_mz_index-wavelet_length << "\t" << scan_size << "\t" << new_length << std::endl;
					for (UInt j=i+c_mz_index-wavelet_length; j<scan_size && k<new_length; ++j, ++k)
					{
						//ofile << scan[j].getMZ() << "\t" <<  helpv[r]*psi_[k] << std::endl;
						prod_[k] = helpv[r]*scan[j].getIntensity()*psi_[k];
						xs_[k] = scan[j].getMZ();
					};

					if (k< new_length) // I.e. we have an overlapping wavelet
					{
						sums = 0;
						if (i==0)
						{
							std::cout << "ARGS: overlapping wavelet and i-1 is 0 !!!!!!!! r: " << r << std::endl;
							std::cout << c_mz_index << "\t" << wavelet_length << std::endl;
							std::cout << k << "\t" << new_length << std::endl;
						}
						max_position_scan = transforms[c][i-1].getMZ()+av_MZ_spacing_;
						break;
					}
					else
					{
						//sums += chordTrapezoidRule_ (xs_, prod_);
						for (UInt k=0; k<prod_.size(); ++k)
						{
							sums += prod_[k];
						};
					};
				};

				//ofile.close();

				//Store the current convolution result
				PeakType& c_peak1 (transforms[c][i]);
				c_peak1.setIntensity(sums);
				c_peak1.setMZ(max_position_scan);
				transforms[c][i].setIntensity(sums);
				//std::cout << "Positions: " << scan[i].getMZ() << " " << max_position_scan << std::endl;
				transforms[c][i].setMZ(max_position_scan);	
			}
		}

		#ifdef DEBUG_FEATUREFINDER
			for (c=0; c<max_charge; ++c)
			{
				std::stringstream name; name << "trans_cmarr_" << scan.getRT() << "_" << scan.begin()->getMZ() <<  "_" << c+1 << "_s" << CoupledMarrWavelet::getSigma() << ".dat\0"; 
				std::ofstream ofile (name.str().c_str());
				for (UInt i=0; i<transforms[c].size(); ++i)
					ofile << transforms[c][i].getMZ() << "\t" <<  transforms[c][i].getIntensity() << std::endl;
				ofile.close();
			};
		#endif

		return;
	}

	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates,
		const MSSpectrum<PeakType>& ref, const UInt scan_index, const DoubleReal ampl_cutoff, const bool use_cmarr)
	{
		UInt cands_size=candidates.size(), signal_size=candidates[0].size(); 
		typename ConstRefVector<MSSpectrum<PeakType> >::iterator iter;
		typename MSSpectrum<PeakType>::const_iterator iter_start, iter_end, iter_p, seed_iter, iter2;
		DoubleReal mz_cutoff, seed_mz, c_av_intens=0, c_score=0, c_sd_intens=0, threshold=0, help_mz, share, share_pos, bwd, fwd;
	 	UInt MZ_start, MZ_end;
		
		//For all charges do ...
		for (Size c=0; c<cands_size; ++c)
		{

			MSSpectrum<PeakType> diffed (candidates[c]);
			diffed[0].setIntensity(0); diffed[signal_size-1].setIntensity(0);
			for (unsigned int i=0; i<signal_size-2; ++i)
			{
				share = diffed[i+1].getIntensity(), share_pos = diffed[i+1].getMZ();
				bwd = (share-diffed[i].getIntensity())/(share_pos-diffed[i].getMZ());
				fwd = (diffed[i+2].getIntensity()-share)/(diffed[i+2].getMZ()-share_pos);
				if (bwd<0 || fwd>0)
				{
					diffed[i+1].setIntensity(0);
				};
			};

			//DPeakConstReferenceArray<MSSpectrum<PeakType> > c_sorted_candidate_ (candidates[c].begin(), candidates[c].end());
			ConstRefVector<MSSpectrum<PeakType> > c_sorted_candidate_ (diffed.begin(), diffed.end());

			//Sort the transform in descending order according to the intensities present in the transform
			c_sorted_candidate_.sortByIntensity();
			//sort (c_sorted_candidate_.begin(), c_sorted_candidate_.end(), intensityComparator<PeakType>);

			std::vector<UInt> processed (signal_size, 0);

			if (ampl_cutoff < 0)
			{
				threshold=0;
			}
			else
			{
				c_av_intens = getAvIntens_ (candidates[c]);
				c_sd_intens = getSdIntens_ (candidates[c], c_av_intens);
				threshold=ampl_cutoff*c_sd_intens + c_av_intens;
			};

			for (iter=c_sorted_candidate_.end()-1; iter != c_sorted_candidate_.begin(); --iter)
			{					
				if (iter->getIntensity() <= 0)
				{
					break;
				};
		
				seed_mz = iter->getMZ();
				seed_iter = ref.MZBegin(seed_mz);

				if (seed_iter == ref.end() || processed[distance(ref.begin(), seed_iter)])
				{
					continue;
				};
				 
				mz_cutoff = getMzPeakCutOffAtMonoPos (seed_mz, c+1.);
				//Mark the region as processed
				//Do not move this further down, since we have to mark this as processed in any case,
				//even when score <=0; otherwise we would look around the maximum's position unless
				//any significant point is found
				iter_start = ref.MZBegin(seed_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));
				iter_end = ref.MZEnd(seed_mz+mz_cutoff/(c+1.));
				if (iter_end == ref.end())
				{
					--iter_end;
				};
								
				MZ_start = distance (ref.begin(), iter_start);
				MZ_end = distance (ref.begin(), iter_end);
		
				memset (&(processed[MZ_start]), 1, sizeof(UInt)*(MZ_end-MZ_start+1));

				c_score = scoreThis_ (candidates[c], getNumPeakCutOff(seed_mz, c+1.), seed_mz, c, iter->getIntensity(), threshold);

				if (c_score <= 0)
				{
					continue;
				};


				//Push the seed into its corresponding box (or create a new one, if necessary)
				//Do ***NOT*** move this further down!
				
				push2TmpBox_ (seed_mz, scan_index, c, c_score, iter->getIntensity(), ref.getRT(), MZ_start, MZ_end);
	
				help_mz = seed_mz - Constants::IW_NEUTRON_MASS/(c+1.);
				iter2 = candidates[c].MZBegin (help_mz);
				if (iter2 == candidates[c].end() || iter2 == candidates[c].begin())
				{
					continue;
				};

				if (fabs(iter2->getMZ()-seed_mz) > 0.5*Constants::IW_NEUTRON_MASS/(c+1.))
				//In the other case, we are too close to the peak, leading to incorrect derivatives.
				{
					if (iter2 != candidates[c].end())
					{
						push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, getLinearInterpolation(iter2-1, help_mz, iter2), ref.getRT(), MZ_start, MZ_end);
					};
				};
			
				help_mz = seed_mz + Constants::IW_NEUTRON_MASS/(c+1.);
					}
				iter2 = candidates[c].MZBegin (help_mz);
				if (iter2 == candidates[c].end()|| iter2 == candidates[c].begin())
				{
					continue;
				};

				if (fabs(iter2->getMZ()-seed_mz) > 0.5*Constants::IW_NEUTRON_MASS/(c+1.))
				//In the other case, we are too close to the peak, leading to incorrect derivatives.
				{
					if (iter2 != candidates[c].end())
					{
						push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, getLinearInterpolation(iter2-1, help_mz, iter2), ref.getRT(), MZ_start, MZ_end);

		clusterSeeds_(candidates, ref, scan_index, candidates.size(), use_cmarr);
	}

	
	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::estimatePeakCutOffs_ (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge)
	{
		std::vector<DoubleReal> x, y;
		UInt peak_cutoff=0;
		DoubleReal step = 100.0;
		while (x.size() < 3) //otherwise: no fit possible.
			IsotopeWavelet::getAveragine (i, &peak_cutoff);
			x.push_back (i);
			y.push_back (peak_cutoff);
		};

		if (x.size() <3)
		{
			x.clear(); y.clear();
			DoubleReal step = (max_mz*max_charge-min_mz)/3.;
			for (DoubleReal i=min_mz; i<max_mz*max_charge; i+=step)
			{
				IsotopeWavelet::getAveragine (i, &peak_cutoff);
				x.push_back (i);
				y.push_back (peak_cutoff);
			};
			step /= 2.0;
		};

		Math::LinearRegression regress;
		regress.computeRegression (0.95, x.begin(), x.end(), y.begin());
		peak_cutoff_intercept_ = regress.getIntercept();
		peak_cutoff_slope_ = regress.getSlope();
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::sampleTheIsotopeWavelet (const MSSpectrum<PeakType>& scan, const UInt wavelet_length, 
		const UInt mz_index, const DoubleReal offset, const UInt charge)
	{
		UInt scan_size = scan.size();
		//Usually we should compute something like: c_pos*charge-mode*charge*PROTON_MASS, 
		//where mode denotes the recording mode +1, -1
		//normally this has only a tiny effect on lambda
		DoubleReal c_pos=scan[mz_index].getMZ(), lambda=IsotopeWavelet::getLambdaQ(c_pos*charge);

		if (mz_index+wavelet_length >= scan_size)
		{
			psi_ = std::vector<double> (wavelet_length, 0);
			return;
		}

		DoubleReal cum_spacing=offset, max_spacing=offset;
		c_mzs_[0] = scan[mz_index].getMZ();
		for (Size j=1; j<wavelet_length+1; ++j)
		{
			c_mzs_[j] = scan[mz_index+j].getMZ();
			c_spacings_[j-1] = c_mzs_[j]-c_mzs_[j-1];
			c_spacings_[j-1] = (c_spacings_[j-1] > 0) ? c_spacings_[j-1] : av_MZ_spacing_;
			max_spacing += c_spacings_[j-1];
		}

		//Building up (sampling) the wavelet
		DoubleReal tz1;
		DoubleReal inv_table_steps = IsotopeWavelet::getInvTableSteps();
		DoubleReal max_tz1 = max_spacing*charge+1;

		if ( ceil(max_tz1*inv_table_steps) < IsotopeWavelet::getGammaTableMaxIndex() &&  ceil(lambda*inv_table_steps) < IsotopeWavelet::getExpTableMaxIndex())
		{
			for (UInt j=0; j<wavelet_length; ++j)
			{
				tz1=cum_spacing*charge+1;
				//psi_[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambda (lambda, tz1) : 0;
				//psi_[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambdaExact (lambda, tz1) : 0;
				cum_spacing += c_spacings_[j];
			}
		}
		else
		{
			for (UInt j=0; j<wavelet_length; ++j)
			{
				tz1=cum_spacing*charge+1;
				//psi_[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambdaExtrapol (lambda, tz1) : 0;
				//psi_[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambdaExact (lambda, tz1) : 0;
				cum_spacing += c_spacings_[j];
			};
		};

		#ifdef DEBUG_FEATUREFINDER
			if (trunc(c_mzs_[0]) == 573)// || trunc(c_mzs_[0]) == 1000 || trunc(c_mzs_[0]) == 1808 || trunc(c_mzs_[0]) == 2000 || trunc(c_mzs_[0]) == 3000)
			{
				std::stringstream stream; stream << "wavelet_" << c_mzs_[0] << "_" << charge << ".dat\0";
				std::ofstream ofile (stream.str().c_str());
				for (UInt i=0; i<wavelet_length; ++i)
				{
					ofile << scan[mz_index+i].getMZ() << "\t" << psi_[i] << std::endl;
				}
				ofile.close();
			};
		#endif

	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::sampleTheCMarrWavelet_ (const MSSpectrum<PeakType>& scan, const Int wavelet_length, 
		const Int mz_index, const UInt charge)
	{
		Int scan_size = scan.size();

		if (mz_index+wavelet_length >= scan_size || mz_index-wavelet_length<=0)
		{
			psi_ = std::vector<double> (2*wavelet_length+1, 0);
			return;
		}

		c_mzs_[wavelet_length] = scan[mz_index].getMZ();
		for (Int j=0; j<wavelet_length; ++j) //right side
		{		
			c_mzs_[wavelet_length+j+1] = scan[mz_index+j+1].getMZ();

			//std::cout << "Spacings R: " <<  c_mzs_[wavelet_length+j] << "\t" << c_mzs_[wavelet_length+j+1] << 
			//	"\t" << c_mzs_[wavelet_length+j+1]-c_mzs_[wavelet_length+j] << std::endl;

			c_spacings_[wavelet_length+j] = c_mzs_[wavelet_length+j+1]-c_mzs_[wavelet_length+j];
			c_spacings_[wavelet_length+j] = (c_spacings_[wavelet_length+j] > 0) ? c_spacings_[wavelet_length+j] : av_MZ_spacing_;
		}

		c_mzs_[0] = scan[mz_index-wavelet_length].getMZ();
		for (Int j=0; j<wavelet_length; ++j) //left side
		{
			c_mzs_[j+1] = scan[mz_index-wavelet_length+j+1].getMZ();

			//std::cout << "Spacings L: " <<  c_mzs_[j] << "\t" << c_mzs_[j+1] <<  "\t" << c_mzs_[j+1]-c_mzs_[j] <<std::endl;


			c_spacings_[j] = c_mzs_[j+1]-c_mzs_[j];
			c_spacings_[j] = (c_spacings_[j] > 0) ? c_spacings_[j] : av_MZ_spacing_;
		}


		//std::cout << "wavelet_length: " << wavelet_length << "\t" << scan[mz_index].getMZ() << "\t" << scan[mz_index-wavelet_length].getMZ() << std::endl;
		DoubleReal cum_spacing=scan[mz_index-wavelet_length].getMZ()-scan[mz_index].getMZ();
		//Building up (sampling) the wavelet
		DoubleReal t2z2, t;
		Int j=0;
		for (; j<2*wavelet_length+1; ++j)
		{
			t=cum_spacing;
			t2z2 = pow(cum_spacing*charge,2);
			//std::cout << ::std::setprecision(8) << std::fixed << cum_spacing << "\t" << cum_spacing+scan[mz_index].getMZ() << "\t" << c_spacings_[j] << std::endl;
			psi_[j] = CoupledMarrWavelet::getValueByMass (t2z2, charge);
			//std::cout << "\t" << psi_[j] << std::endl;
			cum_spacing += c_spacings_[j];
		}


		#ifdef DEBUG_FEATUREFINDER
			std::stringstream stream; stream << "cmarr_wavelet_" << scan[mz_index].getMZ() << "_" << charge << ".dat\0"; 
			std::ofstream ofile (stream.str().c_str());
			for (Int i=0; i<2*wavelet_length+1; ++i)
			{
				ofile << scan[mz_index+i-wavelet_length].getMZ() << "\t" << psi_[i] << std::endl;
			}
			ofile.close();
		#endif
	}


	
	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThis_ (const MSSpectrum<PeakType>& candidate,
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff)
	{
		DoubleReal c_score=0, c_val;
		typename MSSpectrum<PeakType>::const_iterator c_left_iter2, c_right_iter2;
		Int signal_size = candidate.size();

		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff-1) -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> positions (end);
		for (Int i=0; i<end; ++i)
		{
			positions[i] =  seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-(i+1)*Constants::IW_HALF_NEUTRON_MASS)/((DoubleReal)c+1);
		};

		Int start_index = distance(candidate.begin(), candidate.MZBegin (positions[0]))-1;
		for (Int v=1; v<=end; ++v, ++p_h_ind)
		{		
			do 
			{
				if (start_index < signal_size-1) ++start_index;
				else return (0);
			} while (candidate[start_index].getMZ() < positions[v-1]);
			
			c_left_iter2 = candidate.begin()+start_index;
			--c_left_iter2;	
			c_right_iter2 = c_left_iter2+1;

			c_val = c_left_iter2->getIntensity() + (c_right_iter2->getIntensity() - c_left_iter2->getIntensity())/(c_right_iter2->getMZ() - c_left_iter2->getMZ()) * (positions[v-1]-c_left_iter2->getMZ()); 	
			
			if (p_h_ind%2 == 1) //I.e. a whole
			{
				c_score -= c_val;
			}
			else
			{
				c_score +=c_val;
			};	

			start_index = distance(candidate.begin(), c_left_iter2);
		};
	
		#ifdef DEBUG_FEATUREFINDER	
			std::ofstream ofile_score ("scores.dat", ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", ios::app);
			ofile_score << c_check_point << "\t" << c_score << std::endl;
			ofile_score.close();
			ofile_check_score.close();
		#endif

		if (c_score <= ampl_cutoff+intens)
		{
			return(0);
		};

		return (c_score);
	}


	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThisAlternative_ (const MSSpectrum<PeakType>& candidate, 
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) 
	{
		const UInt c_check_point = 1606, c_check_charge = 1; 

		IsotopeDistribution::ContainerType heights = IsotopeWavelet::getAveragine (seed_mz*(c+1.), NULL);
		DoubleReal max=0;
		for (unsigned int i=0; i<heights.size(); ++i)
		{
			max = std::max (heights[i].second, max);
		};
		for (unsigned int i=0; i<heights.size(); ++i)
		{
			heights[i].second /= max/intens;
		};

		std::vector<DoubleReal> theo_conv (4*peak_cutoff-2), pos(4*peak_cutoff-2), real_conv; DoubleReal c_val, c_pos, c_tz, lambda = IsotopeWavelet::getLambdaQ(seed_mz*(c+1));
	
		for (unsigned int i=0; i<peak_cutoff; ++i) //wavelet
		{
			c_val = 0;
			c_pos = seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-(i)*Constants::IW_NEUTRON_MASS)/((DoubleReal)c+1);
			pos[i*2+1] = c_pos;
			for (unsigned int j=0; j<=i; ++j) //pattern
			{
				c_tz = Constants::IW_QUARTER_NEUTRON_MASS+(peak_cutoff-1-(i-j))*Constants::IW_NEUTRON_MASS;
				if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
					std::cout << ::std::setprecision(8) << std::fixed << j << "\t" << c_tz << "\t\t" << heights[j].second << "\t" 
						<< (c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0) << std::endl;
				c_val += heights[j].second*(c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0);
			};
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
			{
				std::cout << "***\t" << c_pos << "\t" << c_val << std::endl;
			};	

			theo_conv[i*2+1] = c_val;
		};
		for (unsigned int i=0; i<peak_cutoff; ++i) //wavelet
		{
			c_val = 0;
			c_pos = seed_mz-((peak_cutoff-0.5)*Constants::IW_NEUTRON_MASS-(i)*Constants::IW_NEUTRON_MASS)/((DoubleReal)c+1);
			pos[i*2] = c_pos;
			for (unsigned int j=0; j<=i; ++j) //pattern
			{
				c_tz = Constants::IW_HALF_NEUTRON_MASS+Constants::IW_QUARTER_NEUTRON_MASS+(peak_cutoff-1-(i-j))*Constants::IW_NEUTRON_MASS;
				if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
					std::cout << ::std::setprecision(8) << std::fixed << j << "\t" << c_tz << "\t\t" << heights[j].second << "\t" 
						<< (c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0) << std::endl;
				c_val += heights[j].second*(c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0);
			};
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
			{
				std::cout << "***\t" << c_pos << "\t" << c_val << std::endl;
			};	

			theo_conv[i*2] = c_val;
		};

		for (unsigned int j=1; j<peak_cutoff; ++j) //pattern
		{
			c_val = 0;
			c_pos = seed_mz+(j*Constants::IW_NEUTRON_MASS)/((DoubleReal)c+1);
			pos[2*peak_cutoff+(j-1)*2+1] = c_pos;
			for (unsigned int i=0; i<peak_cutoff-j; ++i) //wavelet
			{
				c_tz = Constants::IW_QUARTER_NEUTRON_MASS+(i)*Constants::IW_NEUTRON_MASS;
				if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
					std::cout << ::std::setprecision(8) << std::fixed << j+i << "\t" << c_tz << "\t\t" << heights[j+i].second << "\t"
						<<  (c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0) << std::endl;
				c_val += heights[j+i].second*(c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0);
			};
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
			{
				std::cout <<  "***\t" << c_pos << "\t" << c_val << std::endl;
			};	
			theo_conv[2*peak_cutoff+(j-1)*2+1] = c_val;
		};

		for (unsigned int j=1; j<peak_cutoff; ++j) //pattern
		{
			c_val = 0;
			c_pos = seed_mz+(-Constants::IW_HALF_NEUTRON_MASS+j*Constants::IW_NEUTRON_MASS)/((DoubleReal)c+1);
			pos[2*peak_cutoff+(j-1)*2] = c_pos;
			for (unsigned int i=0; i<peak_cutoff-j; ++i) //wavelet
			{
				c_tz = Constants::IW_HALF_NEUTRON_MASS+Constants::IW_QUARTER_NEUTRON_MASS+(i)*Constants::IW_NEUTRON_MASS;
				if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
					std::cout << ::std::setprecision(8) << std::fixed << j+i << "\t" << c_tz << "\t\t" << heights[j+i].second << "\t"
						<<  (c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0) << std::endl;
				c_val += heights[j+i].second*(c_tz > 0 ? IsotopeWavelet::getValueByLambdaExact (lambda, c_tz+1) : 0);
			};
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
			{
				std::cout <<  "***\t" << c_pos << "\t" << c_val << std::endl;
			};	
			theo_conv[2*peak_cutoff+(j-1)*2] = c_val;
		};

		std::stringstream stream; stream << "score_" << seed_mz << "_" << c+1 << "\0";
		std::ofstream ofile (stream.str().c_str());
		for (unsigned int i=0; i<theo_conv.size(); ++i)
		{
			ofile << pos[i] << "\t" << theo_conv[i] << std::endl;
		};
		ofile.close();

		DoubleReal c_theo_score=0, c_score=0;

		for (unsigned int i=0; i<theo_conv.size(); ++i)
		{
			c_theo_score += fabs(theo_conv[i]);
		};



		typename MSSpectrum<PeakType>::const_iterator c_left_iter2, c_right_iter2;
		Int signal_size = candidate.size();

		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff)-2; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> positions (end);
		for (Int i=0; i<end; ++i)
		{
			positions[i] =  seed_mz-((peak_cutoff)*Constants::IW_NEUTRON_MASS-(i+1)*Constants::IW_HALF_NEUTRON_MASS)/((DoubleReal)c+1);
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
			{
				std::cout << positions[i] << std::endl;
			};	
		};

		Int start_index = distance(candidate.begin(), candidate.MZBegin (positions[0]))-1;
		for (Int v=1; v<=end; ++v, ++p_h_ind)
		{		
			do 
			{
				if (start_index < signal_size-1) ++start_index;
				else return (0);
			} while (candidate[start_index].getMZ() < positions[v-1]);
			
			c_left_iter2 = candidate.begin()+start_index;
			--c_left_iter2;	
			c_right_iter2 = c_left_iter2+1;

			c_val = c_left_iter2->getIntensity() + (c_right_iter2->getIntensity() - c_left_iter2->getIntensity())/(c_right_iter2->getMZ() - c_left_iter2->getMZ()) * (positions[v-1]-c_left_iter2->getMZ()); 	
			
			if (p_h_ind%2 == 1) //I.e. a whole
			{
				c_score -= c_val;
				real_conv.push_back(c_val);
			}
			else
			{
				c_score +=c_val;
				real_conv.push_back(c_val);
			};	

			start_index = distance(candidate.begin(), c_left_iter2);
		};

		if (c_score <= 0 || c_theo_score <= 0)
		{
			return (0);
		};

		DoubleReal pearson=0, sum_x=0, sum_y=0, sum_x2=0, sum_y2=0, sum_xy=0;
		for (unsigned int i=0; i< theo_conv.size(); ++i)
		{	
			if (trunc(seed_mz) == c_check_point && c+1 ==c_check_charge)
				std::cout << "corr: " <<  pos[i] << "\t" << theo_conv[i] << "\t" << real_conv[i] << std::endl;
			sum_xy += theo_conv[i]*real_conv[i];
			sum_x += theo_conv[i];
			sum_y += real_conv[i];
			sum_x2 += theo_conv[i]*theo_conv[i];
			sum_y2 += real_conv[i]*real_conv[i];
		};

		pearson = ( theo_conv.size() * sum_xy - sum_x*sum_y ) / (sqrt(theo_conv.size()*sum_x2 - sum_x*sum_x) * sqrt(theo_conv.size()*sum_y2 - sum_y*sum_y));

		DoubleReal final_percent = fabs(c_score-c_theo_score)/(2.*(c_score+c_theo_score))*100;
	
		//#ifdef DEBUG_FEATUREFINDER	
			std::ofstream ofile_score ("scores.dat", std::ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", std::ios::app);
			ofile_score << seed_mz << "\t" << c_theo_score << "\t" << c_score << "\t" << final_percent << "\t" << pearson << std::endl;
			ofile_score.close();
			ofile_check_score.close();
		//#endif

		if (pearson < ampl_cutoff)
		{
			return (0);
		};

		/*if (c_score <= ampl_cutoff+intens)
		{
			return(0);
		};*/

		return (c_score);
	}


	/*template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::cudaScoreThis_ (const MSSpectrum<PeakType>& candidate, 
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) 
	{
		DoubleReal c_score=0, c_check_point=-1, c_val;
		std::pair<int, int> c_between_left, c_between_right;
		typename MSSpectrum<PeakType>::const_iterator c_left_iter, c_right_iter, c_iter;

		DoubleReal leftBound, rightBound;
		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff-1) -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		UInt i=0; bool broken=false; DoubleReal dist, lambda, y1, y2, help;
		for (Int v=1; v<end; ++v, ++p_h_ind)
		{
			c_check_point = seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-v*Constants::IW_HALF_NEUTRON_MASS)/((DoubleReal)c+1);

			leftBound = c_check_point;
			rightBound = c_check_point;

			c_left_iter = candidate.MZBegin (c_check_point); 			
			//ugly, but the only way to check it I guess
			if (c_left_iter == candidate.begin()) continue;
			if ( c_left_iter == candidate.end()) return(0);		
			while (c_left_iter->getMZ() > c_check_point)
			{
				if (--c_left_iter == candidate.begin() && c_left_iter->getMZ() > c_check_point)
				{
					broken=true;
					break;
				};
			};
			if (broken)
			{
				broken=false;	
				continue;
			};

			c_right_iter = c_left_iter+1;
			if (c_right_iter == candidate.end()) continue;

			//linear interpolation
			c_val = c_left_iter->getIntensity() + (c_right_iter->getIntensity() - c_left_iter->getIntensity())/(c_right_iter->getMZ() - c_left_iter->getMZ()) * (c_check_point-c_left_iter->getMZ()); 		
			
			if (p_h_ind%2 == 1) //I.e. a whole
			{
				c_score -= c_val;
			}
			else
			{
				c_score +=c_val;
			};
		};
	
		#ifdef DEBUG_FEATUREFINDER
			std::ofstream ofile_score ("scores.dat", ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", ios::app);
			ofile_score << c_check_point << "\t" << c_score << std::endl;
			ofile_score.close();
			ofile_check_score.close();
		#endif

		if (c_score <= ampl_cutoff+intens)
		{
			return(0);
		};

		return (c_score);
	}*/


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvMZSpacing_ (const MSSpectrum<PeakType>& scan)//, Int start_index, Int end_index) 
	{
		std::vector<DoubleReal> diffs (scan.size()-1, 0);
		for (UInt i=0; i<scan.size()-1; ++i)
		{
			 diffs[i]= scan[i+1].getMZ() - scan[i].getMZ();
		};

		sort(diffs.begin(), diffs.end());
		DoubleReal av_MZ_spacing=0;
		for (UInt i=0; i<diffs.size()/2; ++i)
		{
			av_MZ_spacing += diffs[i];
		};

		return (av_MZ_spacing / (diffs.size()/2));
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvIntens_ (const MSSpectrum<PeakType>& scan)
	{
		DoubleReal av_intens=0;
		for (Size i=0; i<scan.size(); ++i)
		{
			if (scan[i].getIntensity() >= 0)
			{
				av_intens += scan[i].getIntensity();
			}
		};
		return (av_intens / (double)scan.size());
	}

	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getPointerAvIntens_ (const MSSpectrum<PeakType*>& scan) 
	{ 
		DoubleReal av_intens=0;
		for (UInt i=0; i<scan.size(); ++i)
		{
			if (scan[i]->getIntensity() >= 0)
			{
				av_intens += scan[i]->getIntensity();
			}
		};
		return (av_intens / (double)scan.size());
	}
	
	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean)
	{
		DoubleReal res=0, intens;
		for (Size i=0; i<scan.size(); ++i)
		{
			if (scan[i].getIntensity() >= 0)
			{
				intens = scan[i].getIntensity();
				res += (intens-mean)*(intens-mean);
			};
		};
		return (sqrt(res/(double)(scan.size()-1)));
	}
	
	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getPointerSdIntens_ (const MSSpectrum<PeakType*>& scan, const DoubleReal mean) 
	{
		DoubleReal res=0, intens;
		for (UInt i=0; i<scan.size(); ++i)
		{		
			if (scan[i]->getIntensity() >= 0)
			{
				intens = scan[i]->getIntensity();
				res += (intens-mean)*(intens-mean);
			};
		};
		return (sqrt(res/(double)(scan.size()-1)));
	}

	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getCubicInterpolatedValue_ (const std::vector<DoubleReal>& x,
		const DoubleReal xi, const std::vector<DoubleReal>& y)
	{
		gsl_spline_init (spline_, &x[0], &y[0], x.size());
		DoubleReal yi = gsl_spline_eval (spline_, xi, acc_);
		return (yi);
	}


	template <typename PeakType>
	std::pair<int, int> IsotopeWaveletTransform<PeakType>::getNearBys_ (const MSSpectrum<PeakType>& signal, const DoubleReal mz,
		const UInt start) const
	{
		for (Size i=start; i<signal.size(); ++i)
		{
			if (signal[i].getMZ() > mz)
			{
				if (i>start) //everything's fine
				{
					return (std::pair<int, int> (i-1, i));
				}
				else //wrong boundaries!
				{
					break;
				};
			};
		};

		//not found
		return (std::pair<int, int> (-1, -1));
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::push2Box_ (const DoubleReal mz, const UInt scan, UInt c, 
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end)
	{
		//std::cout << "in final push2Box: " << mz << std::endl;
	
		typename std::multimap<DoubleReal, Box>::iterator upper_iter = open_boxes_.upper_bound(mz);
		typename std::multimap<DoubleReal, Box>::iterator lower_iter;

		lower_iter = open_boxes_.lower_bound(mz);
		if (lower_iter != open_boxes_.end())
		{
			//Ugly, but necessary due to the implementation of STL lower_bound
			if (mz != lower_iter->first && lower_iter != open_boxes_.begin())
			{
				--lower_iter;
			};
		};

		typename std::multimap<DoubleReal, Box>::iterator insert_iter;
		bool create_new_box=true;
		if (lower_iter == open_boxes_.end()) //I.e. there is no open Box for that mz position
		{
			//There is another special case to be considered here:
			//Assume that the current box contains only a single element that is (slightly) smaller than the new mz value,
			//then the lower bound for the new mz value is box.end and this would usually force a new entry
			if (!open_boxes_.empty())
			{
				if (fabs((--lower_iter)->first - mz) < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) //matching box
				{
					create_new_box=false;
					insert_iter = lower_iter;
				};
			}
			else
			{
				create_new_box=true;
			}
		}
		else
		{
			if (upper_iter == open_boxes_.end() && fabs(lower_iter->first - mz) < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) //Found matching Box
			{
				insert_iter = lower_iter;
				create_new_box=false;
			}
			else
			{
				create_new_box=true;
			};
		};


		if (upper_iter != open_boxes_.end() && lower_iter != open_boxes_.end())
		{
			//Here is the question if you should figure out the smallest charge .... and then

			//Figure out which entry is closer to m/z
			DoubleReal dist_lower = fabs(lower_iter->first - mz);
			DoubleReal dist_upper = fabs(upper_iter->first - mz);
			dist_lower = (dist_lower < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) ? dist_upper : INT_MAX;

			if (dist_lower>=Constants::IW_HALF_NEUTRON_MASS/(c+1.0) && dist_upper>=Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) // they are both too far away
			{
				create_new_box=true;
			}
			else
			{
				insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;
				create_new_box=false;
			};
		};

		BoxElement element;
		element.c = c; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end;


		if (create_new_box == false)
		{
			double max=intens;
			typename Box::iterator box_iter;
			for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
			{
				if (box_iter->second.max_intens > max)
				{
					max=box_iter->second.max_intens;
				};
			};
			if (max==intens)
			{
				for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
				{
					box_iter->second.max_intens=intens;
				};
			};
			element.max_intens=max;

			std::pair<UInt, BoxElement> help2 (scan, element);
			insert_iter->second.insert (help2);

			//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
			Box replacement (insert_iter->second);

			//We cannot divide both m/z by 2, since we already inserted some m/zs whose weight would be lowered.
			//Also note that we already inserted the new entry, leading to size-1.
			DoubleReal c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;
			c_mz /= ((DoubleReal) insert_iter->second.size());

			//Now let's remove the old and insert the new one
			open_boxes_.erase (insert_iter);
			std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help3 (c_mz, replacement);
			open_boxes_.insert (help3);
		}
		else
		{
			element.max_intens=intens;
			std::pair<UInt, BoxElement> help2 (scan, element);
			std::multimap<UInt, BoxElement> help3;
			help3.insert (help2);
			std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help4 (mz, help3);
			open_boxes_.insert (help4);
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::push2TmpBox_ (const DoubleReal mz, const UInt scan, UInt c, 
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end)
	{
		//std::cout << "Pushing to tmp box: " << mz << "\t scan " << scan+1 << "\t charge " << c+1 << "\t score " << score << std::endl;

		std::multimap<DoubleReal, Box>& tmp_box (tmp_boxes_->at(c));

		lower_iter = tmp_box.lower_bound(mz);
		if (lower_iter != tmp_box.end())
		{
			//Ugly, but necessary due to the implementation of STL lower_bound
			if (mz != lower_iter->first && lower_iter != tmp_box.begin())
			{
				--lower_iter;
			};			
		};

		typename std::multimap<DoubleReal, Box>::iterator insert_iter;
		bool create_new_box=true;
		if (lower_iter == tmp_box.end()) //I.e. there is no tmp Box for that mz position
		{
			//There is another special case to be considered here:
			//Assume that the current box contains only a single element that is (slightly) smaller than the new mz value,
			//then the lower bound for the new mz value is box.end and this would usually force a new entry
			if (!tmp_box.empty())
			{
				if (fabs((--lower_iter)->first - mz) < Constants::IW_HALF_NEUTRON_MASS/(c+1.)) //matching box
				{
					create_new_box=false;
					insert_iter = lower_iter;
				};
			}
			else
			{
				create_new_box=true;
			}
		}
		else
		{
			if (upper_iter == tmp_box.end() && fabs(lower_iter->first - mz) < Constants::IW_HALF_NEUTRON_MASS/(c+1.)) //Found matching Box
			{
				insert_iter = lower_iter;
				create_new_box=false;
			}
			else
			{
				create_new_box=true;
			};
		};


		if (upper_iter != tmp_box.end() && lower_iter != tmp_box.end())
		{
			//Figure out which entry is closer to m/z
			DoubleReal dist_lower = fabs(lower_iter->first - mz);
			DoubleReal dist_upper = fabs(upper_iter->first - mz);
			dist_lower = (dist_lower < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) ? dist_upper : INT_MAX;

			if (dist_lower>=Constants::IW_HALF_NEUTRON_MASS/(c+1.0) && dist_upper>=Constants::IW_HALF_NEUTRON_MASS/(c+1.0)) // they are both too far away
			{
				create_new_box=true;
			}
			else
			{
				insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;
				create_new_box=false;
			};
		};

		BoxElement element;
		element.c = c; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end;

		if (create_new_box == false)
		{
			double max=intens;
			typename Box::iterator box_iter;
			for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
				if (box_iter->second.max_intens > max)
					max=box_iter->second.max_intens;
			if (max==intens)
			for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
				box_iter->second.max_intens=intens;
			element.max_intens=max;

			std::pair<UInt, BoxElement> help2 (scan, element);
			insert_iter->second.insert (help2);

			//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
			Box replacement (insert_iter->second);

			//We cannot divide both m/z by 2, since we already inserted some m/zs whose weight would be lowered.
			//Also note that we already inserted the new entry, leading to size-1.
			DoubleReal c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;
			c_mz /= ((DoubleReal) insert_iter->second.size());

			//Now let's remove the old and insert the new one
			tmp_box.erase (insert_iter);
			std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help3 (c_mz, replacement);
			tmp_box.insert (help3);
		}
		else
		{
			element.max_intens=intens;
			std::pair<UInt, BoxElement> help2 (scan, element);
			std::multimap<UInt, BoxElement> help3;
			help3.insert (help2);

			std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help4 (mz, help3);
			tmp_box.insert (help4);
		};
	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::updateBoxStates (const MSExperiment<PeakType>& map, const Size scan_index, const UInt RT_interleave,
		const UInt RT_votes_cutoff)
	{
		typename std::multimap<DoubleReal, Box>::iterator iter, iter2;

		for (iter=open_boxes_.begin(); iter!=open_boxes_.end(); )
		{
			//For each Box we need to figure out, if and when the last RT value has been inserted
			//If the Box his unchanged since RT_interleave_ scans, we will close the Box.
			UInt lastScan = (--(iter->second.end()))->first;
			if (scan_index - lastScan > RT_interleave) //I.e. close the box!
			{
				iter2 = iter;
				++iter2;
				//Please do **NOT** simplify the upcoming lines.
				//The 'obvious' overhead is necessary since the object represented by iter might be erased
				//by push2Box which might be called by extendBox_.
				if (iter->second.size() >= RT_votes_cutoff)
				{
					extendBox_ (map, iter->second);
					iter = iter2;
					closed_boxes_.insert (*(--iter));
				}
				open_boxes_.erase (iter);
				iter=iter2;
			}
			else
			{
				++iter;
			};
		};

	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::extendBox_ (const MSExperiment<PeakType>& map, const Box box)
	{
		#ifdef DEBUG_FEATUREFINDER
			std::cout << "**** CHECKING FOR BOX EXTENSIONS ****" << std::endl;
		#endif

		//Determining the elution profile
		typename Box::const_iterator iter;
		std::vector<DoubleReal> elution_profile (box.size());
		UInt index=0;
		for (iter=box.begin(); iter != box.end(); ++iter, ++index)
		{
			for (Size i=iter->second.MZ_begin; i!= iter->second.MZ_end; ++i)
			{
				elution_profile[index] += map[iter->second.RT_index][i].getIntensity();
			};
			elution_profile[index] /= iter->second.MZ_end-iter->second.MZ_begin+1.;
		};

		DoubleReal max=0;
		Int max_index=INT_MIN;
		for (Size i=0; i<elution_profile.size(); ++i)
		{
			if (elution_profile[i] > max)
			{
				max_index = i;
				max = elution_profile[i];
			};
		};

		Int max_extension = (Int)(elution_profile.size()) - 2*max_index;

		DoubleReal av_elution=0;
		for (Size i=0; i<elution_profile.size(); ++i)
		{
			av_elution += elution_profile[i];
		};
		av_elution /= (DoubleReal)elution_profile.size();

		DoubleReal sd_elution=0;
		for (Size i=0; i<elution_profile.size(); ++i)
		{
			sd_elution += (av_elution-elution_profile[i])*(av_elution-elution_profile[i]);
		};
		sd_elution /= (DoubleReal)(elution_profile.size()-1);
		sd_elution = sqrt(sd_elution);

		//Determine average m/z monoisotopic pos
	 	DoubleReal av_mz=0;
		for (iter=box.begin(); iter != box.end(); ++iter, ++index)
		{
			av_mz += iter->second.mz;
			#ifdef DEBUG_FEATUREFINDER
				std::cout << iter->second.RT << "\t" << iter->second.mz << "\t" << iter->second.c+1 << std::endl;
			#endif
		};
		av_mz /= (DoubleReal)box.size();


		//Boundary check
		if ((Int)(box.begin()->second.RT_index)-1 < 0)
		{
			return;
		};

		UInt pre_index =  box.begin()->second.RT_index-1;
		typename MSSpectrum<PeakType>::const_iterator c_iter =	map[pre_index].MZBegin(av_mz);
		DoubleReal pre_elution=0;
		
		DoubleReal mz_start = map[pre_index+1][box.begin()->second.MZ_begin].getMZ();
		DoubleReal mz_end = map[pre_index+1][box.begin()->second.MZ_end].getMZ();
		
		typename MSSpectrum<PeakType>::const_iterator mz_start_iter = map[pre_index].MZBegin(mz_start), mz_end_iter = map[pre_index].MZBegin(mz_end);
		for (typename MSSpectrum<PeakType>::const_iterator mz_iter=mz_start_iter; mz_iter != mz_end_iter; ++mz_iter)
		{
			pre_elution += mz_iter->getIntensity();
		};


		//Do we need to extend at all?
		if (pre_elution <= av_elution-2*sd_elution)
		{
			return;
		};

		Int c_index = max_extension;
		Int first_index = box.begin()->second.RT_index;
		for (Int i=1; i<max_extension; ++i)
		{
			c_index = first_index-i;
			if (c_index < 0)
			{
				break;
			};

			//CHECK Majority vote for charge???????????????
			#ifdef DEBUG_FEATUREFINDER
				std::cout << box.begin()->second.RT << "\t" << av_mz << "\t" << box.begin()->second.c+1 << "\t" << " extending the box " << std::endl;
			#endif

			push2Box_ (av_mz, c_index, box.begin()->second.c, box.begin()->second.score, c_iter->getIntensity(),
				map[c_index].getRT(), box.begin()->second.MZ_begin, box.begin()->second.MZ_end);
		};
	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::clusterSeeds_ (const std::vector<MSSpectrum<PeakType> >& candidates,
		const MSSpectrum<PeakType>& ref,  const UInt scan_index, const UInt max_charge, const bool use_cmarr) 
	{
		typename std::multimap<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
	 	DoubleReal c_mz, c_RT, av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, count=0, av_abs_intens;
		UInt num_o_feature, l_mz, r_mz;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (UInt c=0; c<max_charge; ++c)
		{
			std::vector<BoxElement> final_box;
			for (iter=tmp_boxes_->at(c).begin(); iter!=tmp_boxes_->at(c).end(); ++iter)
			{
				Box& c_box = iter->second;
				av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, count=0, l_mz=INT_MAX, r_mz=0, av_abs_intens=0;
				//Now, let's get the RT boundaries for the box
				for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
				{
					c_mz = box_iter->second.mz;
					c_RT = box_iter->second.RT;
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_abs_intens += fabs(box_iter->second.intens);
					av_mz += c_mz*box_iter->second.intens;

					if (l_mz > box_iter->second.MZ_begin) l_mz=box_iter->second.MZ_begin;
					if (r_mz < box_iter->second.MZ_end) r_mz=box_iter->second.MZ_end;

					++count;
				};

				if (count == 0)
				{
					continue;
				};

				av_max_intens = c_box.begin()->second.max_intens;
				av_intens /= count;
				//in contrast to the key entry of tmp_box_, this mz average is weighted by intensity
				av_mz /= av_abs_intens;
				av_score /= count;
				av_RT = c_box.begin()->second.RT;

				BoxElement c_box_element;
				c_box_element.mz = av_mz;
				c_box_element.c = c;
				c_box_element.score = av_score;
				c_box_element.intens = av_intens;
				av_abs_intens += fabs(box_iter->second.intens);
				c_box_element.max_intens = av_max_intens;
				c_box_element.RT=av_RT;
				final_box.push_back(c_box_element);
			};

			num_o_feature = final_box.size();
			if (num_o_feature == 0)
			{
				tmp_boxes_->at(c).clear();
				return;
			};

			//Computing the derivative
			std::vector<DoubleReal> bwd_diffs(num_o_feature, 0);

			bwd_diffs[0]=0;
			for (Size i=1; i<num_o_feature; ++i)
			{
				bwd_diffs[i] = (final_box[i].max_intens-final_box[i-1].max_intens)/(final_box[i].mz-final_box[i-1].mz);
			};

			#ifdef DEBUG_FEATUREFINDER
				std::ofstream ofile_bwd ("bwd.dat");
				for (UInt i=0; i<num_o_feature; ++i)
				{
					ofile_bwd << final_box[i].mz << "\t" << bwd_diffs[i] << std::endl;
				};
				ofile_bwd.close();
			#endif

			for (UInt i=0; i<num_o_feature-1; ++i)
			{
				while (i<num_o_feature-1)
				{
					if(final_box[i].score>0) //this has been an helping point
						break;
					++i;
				};

				if (bwd_diffs[i]>0 && bwd_diffs[i+1]<0)
				{
					checkPosition_ (candidates[c], ref, final_box[i].mz, final_box[i].c, scan_index, use_cmarr);	
					continue;
				};
			};
			tmp_boxes_->at(c).clear();
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::clusterCudaSeeds_ (const MSSpectrum<PeakType>& candidates, 
		const MSSpectrum<PeakType>& ref,  const UInt scan_index, const UInt c, const bool use_cmarr) 
	{
		//std::cout << "entering clusterCudaSeeds" << std::endl;
		typename std::map<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
		std::vector<BoxElement> final_box;
	 	DoubleReal c_mz, c_RT, av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, av_abs_intens=0, count=0;
		UInt num_o_feature, l_mz, r_mz;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=tmp_boxes_->at(c).begin(); iter!=tmp_boxes_->at(c).end(); ++iter)
		{	
			//std::cout << "*************" << std::endl;
			
			Box& c_box = iter->second;
			av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, av_abs_intens=0, count=0, l_mz=INT_MAX, r_mz=0;
			//Now, let's get the RT boundaries for the box
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				//std::cout << "in clustering: " << box_iter->second.mz	<< "\t" << box_iter->second.RT << "\t" <<  box_iter->second.intens << std::endl;	
	
				//Do not use this code!
				//In the case of a very low-resolved spectrum the use of additional points improves 
				//the estimated position of the isotope pattern and hence eases to track its
				//retention profile.
				/*if (box_iter->second.score == 0)
				{
					continue;
				};*/

				c_mz = box_iter->second.mz;
				c_RT = box_iter->second.RT;
				av_score += box_iter->second.score;
				av_intens += box_iter->second.intens;
				av_abs_intens += fabs(box_iter->second.intens);
				av_mz += c_mz*fabs(box_iter->second.intens);

				if (l_mz > box_iter->second.MZ_begin) l_mz=box_iter->second.MZ_begin;
				if (r_mz < box_iter->second.MZ_end) r_mz=box_iter->second.MZ_end;

				++count;
			};

			if (count == 0)
			{	
				continue;
			};

			av_max_intens = c_box.begin()->second.max_intens;
			av_intens /= count;		
			av_score /= count; 
			//in contrast to the key entry of tmp_box_, this mz average is weighted by intensity
			av_mz /= av_abs_intens;
			//std::cout << "av_mz is now: " << av_mz << std::endl;

			av_RT = c_box.begin()->second.RT;

			BoxElement c_box_element;
			c_box_element.mz = av_mz;
			c_box_element.c = c;
			c_box_element.score = av_score;
			c_box_element.intens = av_intens;
			c_box_element.max_intens = av_max_intens;
			c_box_element.RT=av_RT;
			final_box.push_back(c_box_element);
		};	

		num_o_feature = final_box.size();
		if (num_o_feature == 0)
		{
			tmp_boxes_->at(c).clear();
			return;
		};

		//Computing the derivatives
		std::vector<DoubleReal> bwd_diffs(num_o_feature, 0);

		bwd_diffs[0]=0;
		for (UInt i=1; i<num_o_feature; ++i)
		{
			bwd_diffs[i] = (final_box[i].max_intens-final_box[i-1].max_intens)/(final_box[i].mz-final_box[i-1].mz);
		};		

		#ifdef DEBUG_FEATUREFINDER
			std::ofstream ofile_bwd ("bwd.dat");
			for (UInt i=0; i<num_o_feature; ++i)
			{
				ofile_bwd << final_box[i].mz << "\t" << bwd_diffs[i] << std::endl;
			};
			ofile_bwd.close();
		#endif
			

		//std::cout << "av_mz before loop: " << av_mz << std::endl;
		for (UInt i=0; i<num_o_feature-1; ++i)
		{	
			while (i<num_o_feature-1)
			{
				if(final_box[i].score>0) //this has been an helping point
					break;
				++i;
			};
			
			//std::cout << "Before dev check: " << final_box[i].mz << std::endl;

			if (bwd_diffs[i]>0 && bwd_diffs[i+1]<0)
			{					
				//std::cout << "Passing to check box: " << final_box[i].mz << std::endl;

				checkPosition_ (candidates, ref, final_box[i].mz, final_box[i].c, scan_index, use_cmarr);	
				continue;
			};
		};

		tmp_boxes_->at(c).clear();
	}



	template <typename PeakType>
	FeatureMap<Feature> IsotopeWaveletTransform<PeakType>::mapSeeds2Features
		(const MSExperiment<PeakType>& map, const UInt max_charge, const UInt RT_votes_cutoff)
	{
		FeatureMap<Feature> feature_map;
		typename std::multimap<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
		UInt best_charge_index; DoubleReal best_charge_score, c_mz, c_RT; UInt c_charge;
		DoubleReal av_intens=0, av_score=0, av_mz=0, av_RT=0, av_max_intens=0; DoubleReal mz_cutoff;
		ConvexHull2D c_conv_hull;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=closed_boxes_.begin(); iter!=closed_boxes_.end(); ++iter)
		{
			Box& c_box = iter->second;
			std::vector<DoubleReal> charge_votes (max_charge, 0), charge_binary_votes (max_charge, 0);

			//Let's first determine the charge
			//Therefor, we can use two types of votes: qualitative ones (charge_binary_votes) or quantitative ones (charge_votes)
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				charge_votes[box_iter->second.c] += box_iter->second.score;
				++charge_binary_votes[box_iter->second.c];
			};

			//... determining the best fitting charge
			best_charge_index=0; best_charge_score=0;
			for (Size i=0; i<max_charge; ++i)
			{
				if (charge_votes[i] > best_charge_score)
				{
					best_charge_index = i;
					best_charge_score = charge_votes[i];
				};
			};

			//Pattern found in too few RT scan
			if (charge_binary_votes[best_charge_index] < RT_votes_cutoff && RT_votes_cutoff <= map.size())
			{
				continue;
			}

			c_charge = best_charge_index + 1; //that's the finally predicted charge state for the pattern

			av_intens=0, av_score=0, av_mz=0, av_RT=0, av_max_intens=0;
			//Now, let's get the RT boundaries for the box
			std::vector<DPosition<2> > point_set;
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				c_mz = box_iter->second.mz;
				c_RT = box_iter->second.RT;

				mz_cutoff = getMzPeakCutOffAtMonoPos (c_mz, c_charge);

				point_set.push_back (DPosition<2> (c_RT, c_mz - Constants::IW_QUARTER_NEUTRON_MASS/(DoubleReal)c_charge)); 
				//-1 since we are already at the first peak and +0.75, since this includes the last peak of the wavelet as a whole
				point_set.push_back (DPosition<2> (c_RT, c_mz + mz_cutoff/(DoubleReal)c_charge)); 
				if (best_charge_index == box_iter->second.c)
				{
					av_max_intens += box_iter->second.max_intens;
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_mz += c_mz*box_iter->second.intens;
				};
				av_RT += c_RT;
			};
			av_intens /= (DoubleReal)charge_binary_votes[best_charge_index];
			av_max_intens /= (DoubleReal)charge_binary_votes[best_charge_index];

			av_mz /= av_intens*(DoubleReal)charge_binary_votes[best_charge_index];
			av_score /= (DoubleReal)charge_binary_votes[best_charge_index];
			av_RT /= (DoubleReal)c_box.size();

			Feature c_feature;
			c_conv_hull = point_set;
			c_feature.setCharge (c_charge);
			c_feature.setConvexHulls (std::vector<ConvexHull2D> (1, c_conv_hull));
			c_feature.setMZ (av_mz);
			c_feature.setIntensity (av_max_intens);
			c_feature.setRT (av_RT);
			c_feature.setQuality (1, av_score);
			feature_map.push_back (c_feature);
		};

		return (feature_map);
	}



	template <typename PeakType>
	bool IsotopeWaveletTransform<PeakType>::checkPosition_ (const MSSpectrum<PeakType>& candidate,
		const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, const UInt c, const UInt scan_index, const bool use_cmarr)
	{
		typename MSSpectrum<PeakType>::const_iterator iter; 
		UInt peak_cutoff;
		peak_cutoff = getNumPeakCutOff (seed_mz, c+1);

		iter = candidate.MZBegin(seed_mz);
		//we can ignore those cases
		if (iter==candidate.begin() || iter==candidate.end()) 
		{
			return (false);
		};

		std::pair<DoubleReal, DoubleReal> reals;
		//Correct the position
		if (use_cmarr)
		{
			reals = correctMZ_ (ref, iter->getMZ(), c);
		}
		else
		{
			reals = std::pair<DoubleReal, DoubleReal> (seed_mz, ref.MZBegin(seed_mz)->getIntensity()); 
		};
		DoubleReal real_mz = reals.first, real_intens = reals.second;		

		if (real_mz <= 0 || real_intens <= 0)
		{
			return (false);
		};
		
		DoubleReal c_score = scoreThis_ (candidate, peak_cutoff, real_mz, c, iter->getIntensity(), 0);
		//DoubleReal c_score = scoreThisAlternative_ (candidate, peak_cutoff, real_mz, c, iter->getIntensity(), 0);

		if (c_score <= 0)
		{
			return (false);
		};

		DoubleReal mz_cutoff = getMzPeakCutOffAtMonoPos (real_mz, c+1.);
		typename MSSpectrum<PeakType>::const_iterator real_l_MZ_iter = ref.MZBegin(real_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));		
		typename MSSpectrum<PeakType>::const_iterator real_r_MZ_iter = ref.MZBegin(real_mz+mz_cutoff/(c+1.));
		if (real_r_MZ_iter == ref.end())
		{
			--real_r_MZ_iter;
		};
	

		UInt real_mz_begin = distance (ref.begin(), real_l_MZ_iter);
		UInt real_mz_end = distance (ref.begin(), real_r_MZ_iter);

		push2Box_ (real_mz, scan_index, c, c_score, real_intens, ref.getRT(), real_mz_begin, real_mz_end);
		return (true);
	}


	template <typename PeakType>
	std::pair<DoubleReal, DoubleReal> IsotopeWaveletTransform<PeakType>::correctMZ_ (const MSSpectrum<PeakType>& ref, DoubleReal c_mz, const DoubleReal c)
	{
		UInt peak_cutoff = getNumPeakCutOff (c_mz, 1);

		typename MSSpectrum<PeakType>::const_iterator liter = ref.MZBegin(c_mz-(peak_cutoff+1)*Constants::IW_NEUTRON_MASS);
		typename MSSpectrum<PeakType>::const_iterator riter = ref.MZBegin(c_mz+(2*peak_cutoff+1)*Constants::IW_NEUTRON_MASS);
		
		//Coarse structured coupled Marr wavelet
		MSSpectrum<PeakType> spec;
		spec.assign (liter, riter);
		CoupledMarrWavelet::setSigma(Constants::IW_QUARTER_NEUTRON_MASS);
		std::vector<MSSpectrum<PeakType> > pwts (1, spec);	
		this->getCMarrTransforms (spec, pwts, 1);
		MSSpectrum<PeakType> c_spec = pwts[0];

		DoubleReal max=INT_MIN;
		typename MSSpectrum<PeakType>::const_iterator s_iter, max_iter;
		for (s_iter=c_spec.begin(); s_iter!=c_spec.end(); ++s_iter)
		{
			if (s_iter->getMZ() < c_mz - 1.5*Constants::IW_NEUTRON_MASS)
			{
				continue;
			};

			if (s_iter->getMZ() > c_mz + 1.5*Constants::IW_NEUTRON_MASS)
			{
				break;
			};

			if(s_iter->getIntensity() > max)
			{	
				max = s_iter->getIntensity();
				max_iter = s_iter;
			};
		};

		if (max == INT_MIN)
		{	
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << Constants::IW_QUARTER_NEUTRON_MASS << ")> rejected: no peak found" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};

		DoubleReal ppms = getPPMs(peptideMassRule(c_mz), max_iter->getMZ());
		if (ppms >= Constants::PEPTIDE_MASS_RULE_THEO_PPM_BOUND)
		{
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << Constants::IW_QUARTER_NEUTRON_MASS << ")> rejected: ppm too large \t" << ppms << " (rule: " << peptideMassRule(c_mz) << " got: " << max_iter->getMZ() << ")" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};
		
		/*DoubleReal ppms = getPPMs(peptideMassRule(c_mz), c_mz);
		if (ppms >= THEO_PPM_BOUND)
		{
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << Constants::IW_QUARTER_NEUTRON_MASS << ")> rejected: ppm too large \t" << ppms << " (rule: " << peptideMassRule(c_mz) << " got: " << max_iter->getMZ() << ")" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};*/

		//Check for non-interleaving pattern overlap to the left (this is the by far easiest case to detect)
		/*typename MSSpectrum<PeakType>::const_iterator off1_r_iter = c_spec.MZBegin(max_iter->getMZ()+Constants::IW_NEUTRON_MASS);
		typename MSSpectrum<PeakType>::const_iterator off1_l_iter = c_spec.MZBegin(max_iter->getMZ()-Constants::IW_NEUTRON_MASS);

		if (off1_r_iter != c_spec.end() && off1_l_iter != c_spec.end())
		{
			DoubleReal perc_exp = fabs(off1_l_iter->getIntensity()-off1_r_iter->getIntensity())/(0.5*(off1_l_iter->getIntensity()+off1_r_iter->getIntensity()));
			if (perc_exp > 0.1)
			{
				std::cout << "potential overlap \t" << perc_exp << "\t(" << c_mz << ") ... ";
				typename MSSpectrum<PeakType>::const_iterator off2_l_iter = c_spec.MZBegin(max_iter->getMZ()-2*Constants::IW_NEUTRON_MASS);
				std::cout << "mzs: " << *off2_l_iter << "\t" << *off1_r_iter << std::endl;
				DoubleReal perc = fabs(off2_l_iter->getIntensity()-off1_r_iter->getIntensity())/(0.5*(off2_l_iter->getIntensity()+off1_r_iter->getIntensity()));
				if (perc_exp > perc)
				{
					std::cout << "yes " << perc << std::endl; 	
					
					typename MSSpectrum<PeakType>::const_iterator real_l_MZ_iter = ref.MZBegin(off1_l_iter->getMZ()-QUARTER_Constants::IW_NEUTRON_MASS/(c+1.));		
					typename MSSpectrum<PeakType>::const_iterator real_r_MZ_iter = ref.MZBegin(off1_l_iter->getMZ()+getPeakCutOff(off1_l_iter->getMZ(), c+1)*Constants::IW_NEUTRON_MASS/(c+1.));

					UInt real_mz_begin = distance (ref.begin(), real_l_MZ_iter);
					UInt real_mz_end = distance (ref.begin(), real_r_MZ_iter);

					push2Box_ (off1_l_iter->getMZ(), 0, c, 2000, off1_l_iter->getIntensity(), ref.getRT(), real_mz_begin, real_mz_end);
				}
				else
				{
					std::cout << "no " << perc << std::endl;
				};
			};
		};*/


		//Fine structured coupled Marr wavelet
		spec.assign (liter, riter);
		CoupledMarrWavelet::setSigma(sigma_);
		pwts = std::vector<MSSpectrum<PeakType> > (1, spec);	
		this->getCMarrTransforms (spec, pwts, 1);
		c_spec = pwts[0];		
	
		std::vector<DoubleReal> bwd_diffs (pwts[0].size(), 0);
		for (UInt i=1; i<pwts[0].size(); ++i)
		{
			bwd_diffs[i] = (pwts[0][i].getIntensity()-pwts[0][i-1].getIntensity())/(pwts[0][i].getMZ()-pwts[0][i-1].getMZ());
		};

		std::vector<UInt> max_indices; // if entry i, then the max should be apparent by i+0.5 index
		for (UInt i=1; i<pwts[0].size()-1; ++i)
		{		
			if (pwts[0][i].getMZ() < c_mz - Constants::IW_QUARTER_NEUTRON_MASS)
			{
				continue;
			};

			if (pwts[0][i].getMZ() > c_mz)
			{
				break;
			};

			if (bwd_diffs[i] > 0 && bwd_diffs[i+1] < 0 && pwts[0][i].getIntensity() > 0)
			{
				max_indices.push_back (i);
			};
		};

		if (max_indices.empty())
		{
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << sigma_ << ")> rejected: no peak found" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};

		Int min_index = INT_MAX; DoubleReal min_ppm = INT_MAX;
		
		std::vector<DoubleReal> dist_to_theo_mz (max_indices.size());
		for (UInt i=0; i<max_indices.size(); ++i)
		{
			dist_to_theo_mz[i] = getPPMs(peptideMassRule(c_mz),  pwts[0][max_indices[i]].getMZ());
			if (dist_to_theo_mz[i] < min_ppm)
			{
				min_ppm = dist_to_theo_mz[i];
				min_index = i;
			};
		};

		if (pwts[0][max_indices[min_index]].getIntensity() <= 0)
		{
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << sigma_ << "))> rejected: picked negative peak" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};

		if (min_ppm >= Constants::PEPTIDE_MASS_RULE_THEO_PPM_BOUND)
		{			
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << sigma_ << ")> rejected: ppm too large \t" << min_ppm << " (rule: " 
				<< peptideMassRule(c_mz) << " got: " << pwts[0][max_indices[min_index]].getMZ() << ")" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};

		std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << sigma_ << ")> ACCEPT \t" << min_ppm << " (rule: " 
			<< peptideMassRule(c_mz) << " got: " << pwts[0][max_indices[min_index]].getMZ() << ")" << std::endl;

		c_mz=pwts[0][max_indices[min_index]].getMZ();

		return (std::pair<DoubleReal, DoubleReal> (c_mz, ref.MZBegin(c_mz)->getIntensity()));
	
	}

} //namespace

#endif
