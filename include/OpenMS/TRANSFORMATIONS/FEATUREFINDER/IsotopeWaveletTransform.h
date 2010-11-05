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
// $Maintainer: Rene Hussong$
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETTRANSFORM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETTRANSFORM_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletConstants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/TransPeak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <cmath>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>

#ifdef OPENMS_HAS_CUDA
	#include <cuda.h>
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
 		* If you just want to find features using the isotope wavelet, take a look at the FeatureFinderAlgorithmIsotopeWavelet class. Usually, you only
 		* have to consider the class at hand if you plan to change the basic implementation of the transform.  */ 
	template <typename PeakType>
	class IsotopeWaveletTransform
	{
		public:


			/** @brief Internally used data structure. */
			struct BoxElement
			{
				DoubleReal mz; //<The monoisotopic position
				UInt c; //<Note, this is not the charge (it is charge-1!!!)
				DoubleReal score; //<The associated score
				DoubleReal intens; //<The transformed intensity at the monoisotopic mass
				DoubleReal ref_intens;
				DoubleReal RT; //<The elution time (not the scan index)
				UInt RT_index; //<The elution time (map) index
				UInt MZ_begin; //<Index
				UInt MZ_end; //<Index
			};

			typedef std::multimap<UInt, BoxElement> Box; ///<Key: RT index, value: BoxElement

			
			/** @brief Internally (only by GPUs) used data structure . 
 				*	It allows efficient data exchange between CPU and GPU and avoids unnecessary memory moves. 
				*	The class is tailored on the isotope wavelet transform and is in general not applicable on similar - but different - situations. */ 
			class TransSpectrum
			{
				friend class IsotopeWaveletTransform;

				public:

					/** Default constructor */
					TransSpectrum()
						: reference_(NULL), trans_intens_(NULL)
					{
					};

					/** Copy constructor */
					TransSpectrum(const MSSpectrum<PeakType>* reference)
						: reference_(reference)
					{
						trans_intens_ = new std::vector<float> (reference_->size(), 0.0);
					};

					/** Destructor */
					virtual ~TransSpectrum()
					{
						delete (trans_intens_);
					};

					virtual void destroy ()
					{
						delete (trans_intens_);
						trans_intens_=NULL;
						delete (reference_);
						trans_intens_=NULL;
					};

					/** Returns the RT value (not the index) of the associated scan. */
					inline DoubleReal getRT () const 
					{ 
						return (reference_->getRT());
					};

					/** Returns the mass-over-charge ratio at index @p i. */
					inline DoubleReal getMZ (const UInt i) const
					{ 
						return ((*reference_)[i].getMZ());
					};
					
					/** Returns the reference (non-transformed) intensity at index @p i. */
					inline DoubleReal getRefIntensity (const UInt i) const
					{ 
						return ((*reference_)[i].getIntensity());
					};
				
					/** Returns the transformed intensity at index @p i. */
					inline DoubleReal getTransIntensity (const UInt i) const
					{ 
						return ((*trans_intens_)[i]);
					};
						
					/** Stores the intensity value @p i of the transform at position @p i. */
					inline void setTransIntensity (const UInt i, const DoubleReal intens)
					{ 
						(*trans_intens_)[i] = intens;
					};

					/** Returns the size of spectra. */
					inline Size size () const
					{
						return (trans_intens_->size());
					};

					/** Returns a pointer to the reference spectrum. */
					inline const MSSpectrum<PeakType>* getRefSpectrum ()
					{
						return (reference_);
					};

					/** Returns a pointer to the reference spectrum. */
					inline const MSSpectrum<PeakType>* getRefSpectrum () const
					{
						return (reference_);
					};

					/** Attention: iterations will only performed over the reference spectrum. 
 						* You will have to use the "distance"-function in order to get the corresponding entry of the transform. */
					inline typename MSSpectrum<PeakType>::const_iterator MZBegin (const DoubleReal mz) const
					{
						return (reference_->MZBegin(mz));
					};
										
					/** Attention: iterations will only performed over the reference spectrum. 
 						* You will have to use the "distance"-function in order to get the corresponding entry of the transform. */
					inline typename MSSpectrum<PeakType>::const_iterator MZEnd (const DoubleReal mz) const
					{
						return (reference_->MZEnd(mz));
					};
										
					/** Attention: iterations will only performed over the reference spectrum. 
 						* You will have to use the "distance"-function in order to get the corresponding entry of the transform. */
					inline typename MSSpectrum<PeakType>::const_iterator end () const
					{
						return (reference_->end());
					};
					
					/** Attention: iterations will only performed over the reference spectrum. 
 						* You will have to use the "distance"-function in order to get the corresponding entry of the transform. */
					inline typename MSSpectrum<PeakType>::const_iterator begin () const
					{
						return (reference_->begin());
					};


				protected:

					const MSSpectrum<PeakType>* reference_; //<The reference spectrum
					std::vector<float>* trans_intens_; //<The intensities of the transform
					
			}; 


			
			/** @brief Constructor.
 				*
 				* @param min_mz The smallest m/z value occurring in your map.
 				* @param max_mz The largest m/z value occurring in your map.
 				* @param max_charge The highest charge state you would like to consider. */
			//IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const UInt max_scan_size=0);
			IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const UInt max_scan_size=0, const bool use_cuda=false, const bool hr_data=false);

			/** @brief Destructor. */
			virtual ~IsotopeWaveletTransform () ;
	

			/** @brief Computes the isotope wavelet transform of charge state @p c. 
 				* @param c_trans The transform.
 				* @param c_ref The reference spectrum. 
 				* @oaram c The charge state minus 1 (e.g. c=2 means charge state 3) at which you want to compute the transform. */	
			virtual void getTransform (MSSpectrum<PeakType>& c_trans, const MSSpectrum<PeakType>& c_ref, const UInt c);

			/** @brief Computes the isotope wavelet transform of charge state @p c. 
 				* @param c_trans The transform.
 				* @param c_ref The reference spectrum. 
 				* @oaram c The charge state minus 1 (e.g. c=2 means charge state 3) at which you want to compute the transform. */	
			virtual void getTransformHighRes (MSSpectrum<PeakType>& c_trans, const MSSpectrum<PeakType>& c_ref, const UInt c);

			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
 				* process is stored internally. Important: Before calling this function, apply updateRanges() to the original map.
 				*
 				* @param candidates A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
 				* charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function @see getTransforms.)
 				* @param ref The reference scan (the untransformed raw data) corresponding to @p candidates.
				* @param c The corrsponding charge state minus 1 (e.g. c=2 means charge state 3)
				* @param scan_index The index of the scan (w.r.t. to some map) currently under consideration.
				* @param ampl_cutoff The thresholding parameter. This parameter is the only (and hence a really important)
 				* parameter of the isotope wavelet transform. On the basis of @p ampl_cutoff the program tries to distinguish between
 				* noise and signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
 				* line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
 				* get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
 				* makes sense to start @p ampl_cutoff=0 or even @p ampl_cutoff=-1,
 				* indicating no thresholding at all. Note that also ampl_cutoff=0 triggers (a moderate) thresholding based on the
 				* average intensity in the wavelet transform. 
 				* @param check_PPMs If enabled, the algorithm will check each monoisotopic mass candidate for its plausibility
 				* by computing the ppm difference between this mass and the averagine model. */ 
			virtual void identifyCharge (const MSSpectrum<PeakType>& candidates, const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt c, 
				const DoubleReal ampl_cutoff, const bool check_PPMs) ;
							
			virtual void initializeScan (const MSSpectrum<PeakType>& c_ref); 
			
			#ifdef OPENMS_HAS_CUDA
				/** @brief Sets up all necessary arrays with correct boundaries and 'worst-case' sizes. 
 					* @param scan The scan under consideration. */	
				virtual int initializeScanCuda (const MSSpectrum<PeakType>& scan, const UInt c=0); 
			
				/** @brief Clean up. */	
				virtual void finalizeScanCuda ();
						
				/** @brief Computes The isotope wavelet transform of charge state (@p c+1) on a CUDA compatible GPU. 
 					* @param c_trans Contains the reference spectrum (already by call) as well as the transformed intensities. 
 					* @param c The charge state minus 1 (e.g. c=2 means charge state 3)*/
				virtual void getTransformCuda (TransSpectrum &c_trans, const UInt c);
				
				/** @brief Essentially the same as its namesake CPU-version, but on a CUDA compatible GPU device. 
  			* @param candidates A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
 				* charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function @see getTransforms.)    
 				* @param c The corrsponding charge state minus 1 (e.g. c=2 means charge state 3)
				* @param scan_index The index of the scan (w.r.t. to some map) currently under consideration. 
 				* @param ampl_cutoff The thresholding parameter. This parameter is the only (and hence a really important)
 				* parameter of the isotope wavelet transform. On the basis of @p ampl_cutoff the program tries to distinguish between 
 				* noise and signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
 				* line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
 				* get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
 				* makes sense to start @p ampl_cutoff=0 or even @p ampl_cutoff=-1, 
 				* indicating no thresholding at all. Note that also ampl_cutoff=0 triggers (a moderate) thresholding based on the 
 				* average intensity in the wavelet transform.  				
 				* * @param check_PPMs If enabled, the algorithm will check each monoisotopic mass candidate for its plausibility
 				* by computing the ppm difference between this mass and the averagine model. */ 
				virtual void identifyChargeCuda (const TransSpectrum& candidates, const UInt scan_index, const UInt c, 
					const DoubleReal ampl_cutoff, const bool check_PPMs);

				/** Sorts the associated spectrum @p by increasing intensities. 
 					* @param sorted The spectrum to be sorted. */	
				virtual int sortCuda (MSSpectrum<PeakType>& sorted);
			#endif	
	

			/** @brief A function keeping track of currently open and closed sweep line boxes.
 				* This function is used by the isotope wavelet feature finder and must be called for each processed scan.
 				* @param map The original map containing the data set to be analyzed.
				* @param scan_index The index of the scan currently under consideration w.r.t. its MS map.
 				* This information is necessary to sweep across the map after each scan has been evaluated.
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class. */
			void updateBoxStates (const MSExperiment<PeakType>& map, const Size scan_index, const UInt RT_interleave,
				const UInt RT_votes_cutoff, const Int front_bound=-1, const Int end_bound=-1) ;

		
			void mergeFeatures (IsotopeWaveletTransform<PeakType>* later_iwt, const UInt RT_interleave, const UInt RT_votes_cutoff); 
	

			/** @brief Filters the candidates further more and maps the internally used data structures to the OpenMS framework.
 				* @param map The original map containing the data set to be analyzed.
 				* @param max_charge The maximal charge state under consideration.
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class.*/
			FeatureMap<Feature> mapSeeds2Features (const MSExperiment<PeakType>& map, const UInt RT_votes_cutoff) ;

			/** @brief Returns the closed boxes. */
			virtual std::multimap<DoubleReal, Box> getClosedBoxes ()
				{ return (closed_boxes_); };
		

			/** @brief Computes a linear (intensity) interpolation. 
 				* @param left_iter The point left to the query. 
 				* @param mz_pos The query point.
 				* @param right_iter The point right to the query. */	
			inline DoubleReal getLinearInterpolation (const typename MSSpectrum<PeakType>::const_iterator& left_iter, const DoubleReal mz_pos, const typename MSSpectrum<PeakType>::const_iterator& right_iter)
			{
				return (left_iter->getIntensity() + (right_iter->getIntensity() - left_iter->getIntensity())/(right_iter->getMZ() - left_iter->getMZ()) * (mz_pos-left_iter->getMZ())); 
			};
			
			/** @brief Computes a linear (intensity) interpolation. 
 				* @param mz_a The m/z value of the point left to the query.  
 				* @param mz_a The intensity value of the point left to the query. 
 				* @param mz_pos The query point.				
 				* @param mz_b The m/z value of the point right to the query.  
 				* @param intens_b The intensity value of the point left to the query. */ 
			inline DoubleReal getLinearInterpolation (const DoubleReal mz_a, const DoubleReal intens_a, const DoubleReal mz_pos, const DoubleReal mz_b, const DoubleReal intens_b)
			{
				return (intens_a + (intens_b - intens_a)/(mz_b - mz_a) * (mz_pos-mz_a)); 
			};

			inline DoubleReal getSigma () const
			{
				return (sigma_);
			};

			inline void setSigma (const DoubleReal sigma)
			{
				sigma_ = sigma;
			};

			virtual void computeMinSpacing (const MSSpectrum<PeakType>& c_ref);
							
			inline DoubleReal getMinSpacing () const
			{
				return (min_spacing_);
			};			
			
			inline Int getMaxScanSize () const
			{
				return (max_scan_size_);
			};


		protected:


			/** @brief Default Constructor.
 				* @note Provided just for inheritance reasons. You should always use the other constructor. */
			IsotopeWaveletTransform () ;


			inline void sampleTheCMarrWavelet_ (const MSSpectrum<PeakType>& scan, const Int wavelet_length, const Int mz_index, const UInt charge);	
			

			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern.
				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined.
 				* @param ampl_cutoff The threshold. */
			virtual DoubleReal scoreThis_ (const TransSpectrum& candidate, const UInt peak_cutoff, 
				const DoubleReal seed_mz, const UInt c, const DoubleReal ampl_cutoff) ;
						
			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score 
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern.
				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined. 
 				* @param ampl_cutoff The threshold. */
			virtual DoubleReal scoreThis_ (const MSSpectrum<PeakType>& candidate, const UInt peak_cutoff, 
				const DoubleReal seed_mz, const UInt c, const DoubleReal ampl_cutoff) ;	


			/** @brief A ugly but necessary function to handle "off-by-1-Dalton predictions" due to idiosyncrasies of the data set
 				* (in comparison to the averagine model)
 				* @param candidate The wavelet transformed spectrum containing the candidate.
 				* @param ref The original spectrum containing the candidate.
 				* @param seed_mz The m/z position of the candidate pattern.
 				* @param c The predicted charge state minus 1 (e.g. c=2 means charge state 3) of the candidate.
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). */
			virtual bool checkPositionForPlausibility_ (const TransSpectrum& candidate, const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, 
				const UInt c, const UInt scan_index, const bool check_PPMs, const DoubleReal transintens) ;
			
			/** @brief A ugly but necessary function to handle "off-by-1-Dalton predictions" due to idiosyncrasies of the data set
 				* (in comparison to the averagine model)
 				* @param candidate The wavelet transformed spectrum containing the candidate. 
 				* @param ref The original spectrum containing the candidate.
 				* @param seed_mz The m/z position of the candidate pattern.
 				* @param c The predicted charge state minus 1 (e.g. c=2 means charge state 3) of the candidate.
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). */
			virtual bool checkPositionForPlausibility_ (const MSSpectrum<PeakType>& candidate, const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, 
				const UInt c, const UInt scan_index, const bool check_PPMs, const DoubleReal transintens) ;
			
			virtual std::pair<DoubleReal, DoubleReal> checkPPMTheoModel_ (const MSSpectrum<PeakType>& ref, const DoubleReal c_mz, const UInt c) ;


			/** @brief Computes the average (transformed) intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getAvIntens_ (const TransSpectrum& scan);
			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getAvIntens_ (const MSSpectrum<PeakType>& scan); 	

			/** @brief Computes the standard deviation (neglecting negative values) of the (transformed) intensities of @p scan. */
			inline DoubleReal getSdIntens_ (const TransSpectrum& scan, const DoubleReal mean) ;
			/** @brief Computes the standard deviation (neglecting negative values) of the intensities of @p scan. */
			inline DoubleReal getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) ;

			/** @brief Inserts a potential isotopic pattern into an open box or - if no such box exists - creates a new one.
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the post-processing (sweep lining).
 				* @param charge The estimated charge state minus 1 (e.g. c=2 means charge state 3) of the pattern. 
 				* @param score The pattern's score.
 				* @param intens The intensity at the monoisotopic peak.
 				* @param rt The retention time of the scan (similar to @p scan, but here: no index, but the real value).
 				* @param MZ_begin The starting index of the pattern (m/z) w.r.t. the current scan.
 				* @param MZ_end The end index (w.r.t. the monoisotopic position!) of the pattern (m/z) w.r.t. the current scan. */
			virtual void push2Box_ (const DoubleReal mz, const UInt scan, UInt c, const DoubleReal score, 
				const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end, const DoubleReal ref_intens) ;

			/** @brief Essentially the same function as @see push2Box_.
 				* In contrast to @see push2Box this function stores its candidates only temporarily. In particular, this
 				* function is only used within a single scan transform. After the wavelet transform is computed on
 				* that scan, all candidates are pushed by this function and finally clustered together by @see clusterSeeds_.
 				* Afterwards, a final push by @see push2Box_ is performed storing the clustered candidates.
 				*
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the post-processing (sweep lining).
 				* @param charge The estimated charge state minus 1 (e.g. c=2 means charge state 3) of the pattern. 
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


			/** @brief Clusters the seeds stored by push2TmpBox_.
 				* @param candidates A isotope wavelet transformed spectrum.
 				* @param ref The corresponding original spectrum (w.r.t. @p candidates).
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). */ 
			void clusterSeeds_ (const TransSpectrum& candidates, const MSSpectrum<PeakType>& ref, 
				const UInt scan_index, const UInt c, const bool check_PPMs) ;
		
			/** @brief Clusters the seeds stored by push2TmpBox_.
 				* @param candidates A isotope wavelet transformed spectrum. 
 				* @param ref The corresponding original spectrum (w.r.t. @p candidates). 
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). */ 
			virtual void clusterSeeds_ (const MSSpectrum<PeakType>& candidates, const MSSpectrum<PeakType>& ref, 
				const UInt scan_index, const UInt c, const bool check_PPMs) ;


			/** @brief A currently still necessary function that extends the box @p box in order to capture also
 				* signals whose isotopic pattern is nearly diminishing 
 				* @param map The experimental map.
 				* @param box The box to be extended. */
			void extendBox_ (const MSExperiment<PeakType>& map, const Box box);

			/** @brief Returns the monoisotopic mass (with corresponding decimal values) we would expect at @p c_mass. 
 				* @param c_mass The mass for which we would like to know the averagine decimal places. */
			inline DoubleReal peptideMassRule (const DoubleReal c_mass) const
			{
				DoubleReal correction_fac = c_mass / Constants::PEPTIDE_MASS_RULE_BOUND;
				DoubleReal old_frac_mass = c_mass - (Int)(c_mass);
				DoubleReal new_mass = ((Int)(c_mass))* (1.+Constants::PEPTIDE_MASS_RULE_FACTOR)-(Int)(correction_fac);
				DoubleReal new_frac_mass = new_mass - (Int)(new_mass);
				
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

			/** @brief Returns the parts-per-million deviation of the masses.
 				* @param mass_a The first mass. 
 				* @param mass_b The second mass. */
			inline DoubleReal getPPMs (const DoubleReal mass_a, const DoubleReal mass_b) const
			{
				return (fabs(mass_a-mass_b)/(0.5*(mass_a+mass_b))*1e6);
			};


			//internally used data structures for the sweep line algorithm
			std::multimap<DoubleReal, Box> open_boxes_, closed_boxes_, end_boxes_, front_boxes_;	//DoubleReal = average m/z position
			std::vector<std::multimap<DoubleReal, Box> >* tmp_boxes_; //for each charge we need a separate container

			DoubleReal av_MZ_spacing_, sigma_;  
			std::vector<DoubleReal> c_mzs_, c_spacings_, psi_, prod_, xs_;
			std::vector<DoubleReal> interpol_xs_, interpol_ys_;

			UInt max_scan_size_, max_num_peaks_per_pattern_, max_charge_;
			bool hr_data_;
			Int from_max_to_left_, from_max_to_right_;
			std::vector<int> indices_;
			
			MSSpectrum<PeakType> c_sorted_candidate_;
			DoubleReal min_spacing_, max_mz_cutoff_;
			std::vector<float> scores_, zeros_;	//hier muessen bestimmt noch einige sachen ins ifdef

			#ifdef OPENMS_HAS_CUDA
				float *h_data_;		
				int *h_pos_;
				UInt largest_array_size_, overall_size_, block_size_, data_length_, to_load_, to_compute_;
				Int num_elements_;
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
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box> > (1);
		av_MZ_spacing_=1;
		max_scan_size_ = 0;
		max_mz_cutoff_ = 3;
		max_num_peaks_per_pattern_ = 3;
		#ifdef OPENMS_HAS_CUDA 
			largest_array_size_ = 0;
			num_elements_ = 0;
			h_data_ = NULL;
			h_pos_ = NULL;
		#endif
	}

	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const UInt max_scan_size, const bool use_cuda, const bool hr_data) 
	{
		max_charge_ = max_charge;
		max_scan_size_ = max_scan_size;
		hr_data_ = hr_data;
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box> > (max_charge);
		if (max_scan_size <= 0) //only important for the CPU
		{
			IsotopeWavelet::init (max_mz, max_charge);			
		};

		av_MZ_spacing_=1;
		max_mz_cutoff_ =  IsotopeWavelet::getMzPeakCutOffAtMonoPos(max_mz, max_charge);
		max_num_peaks_per_pattern_=  IsotopeWavelet::getNumPeakCutOff(max_mz, max_charge);
		
		#ifdef OPENMS_HAS_CUDA 
			if (use_cuda > 0) //only important for the GPU
			{	
				if (hr_data_)
				{
					max_scan_size_ = max_scan_size*(4*(max_mz_cutoff_-1)-1);
				};
				largest_array_size_ =  pow(2, ceil(log(max_scan_size_ +  Constants::CUDA_EXTENDED_BLOCK_SIZE_MAX)/log(2.0)));
				
				cuda_positions_.reserve(largest_array_size_);
		  	cuda_intensities_.reserve(largest_array_size_);	
				indices_.resize (largest_array_size_);
				for (UInt q=0; q<largest_array_size_; ++q)
				{
					indices_[q] = q;
				}; 

				h_data_ = (float*) malloc (largest_array_size_*sizeof(float));		
				h_pos_ = (int*) malloc (largest_array_size_*sizeof(int)); 	
			}
			else
			{			
				h_data_ = NULL;
				h_pos_ = NULL;
				Int size_estimate ((Int)ceil(max_scan_size_/(max_mz-min_mz)));
				Int to_reserve  ((Int)ceil(size_estimate*max_num_peaks_per_pattern_*Constants::IW_NEUTRON_MASS));
				psi_.reserve (to_reserve); //The wavelet
				prod_.reserve (to_reserve); 
				xs_.reserve (to_reserve);
				interpol_xs_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
				interpol_ys_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
			}
		#else
			Int size_estimate ((Int)ceil(max_scan_size_/(max_mz-min_mz)));
			Int to_reserve  ((Int)ceil(size_estimate*max_num_peaks_per_pattern_*Constants::IW_NEUTRON_MASS));
			psi_.reserve (to_reserve); //The wavelet
			prod_.reserve (to_reserve); 
			xs_.reserve (to_reserve);
			interpol_xs_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
			interpol_ys_.resize(Constants::DEFAULT_NUM_OF_INTERPOLATION_POINTS);
		#endif
	}


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::~IsotopeWaveletTransform ()
	{
		#ifdef OPENMS_HAS_CUDA
			if (h_data_ != NULL) free (h_data_);
			if (h_pos_ != NULL) free (h_pos_);			
			h_data_ = NULL;
			h_pos_ = NULL;
		#endif	
	
		delete (tmp_boxes_);
	}
	


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getTransform (MSSpectrum<PeakType>& c_trans, const MSSpectrum<PeakType>& c_ref, const UInt c)
	{
		Int spec_size ((Int)c_ref.size());
		//in the very unlikely case that size_t will not fit to int anymore this will be a problem of course
		//for the sake of simplicity (we need here a signed int) we do not cast at every following comparison individually
		UInt charge = c+1;
		DoubleReal value, T_boundary_left, T_boundary_right, old, c_diff, current, old_pos, my_local_MZ, my_local_lambda, origin, c_mz;

		for(Int my_local_pos=0; my_local_pos<spec_size; ++my_local_pos)
		{
			value=0; T_boundary_left = 0, T_boundary_right = IsotopeWavelet::getMzPeakCutOffAtMonoPos(c_ref[my_local_pos].getMZ(), charge)/(DoubleReal)charge;
			old=0; old_pos=(my_local_pos-from_max_to_left_-1>=0) ? c_ref[my_local_pos-from_max_to_left_-1].getMZ() : c_ref[0].getMZ()-min_spacing_;
			my_local_MZ = c_ref[my_local_pos].getMZ(); my_local_lambda = IsotopeWavelet::getLambdaL(my_local_MZ*charge);
			c_diff=0;
			origin = -my_local_MZ+Constants::IW_QUARTER_NEUTRON_MASS/(DoubleReal)charge;

			for (Int current_conv_pos =  std::max(0, my_local_pos-from_max_to_left_); c_diff < T_boundary_right; ++current_conv_pos)
			{
				if (current_conv_pos >= spec_size)
				{
					value += 0.5*old*min_spacing_;
					break;
				};

				c_mz = c_ref[current_conv_pos].getMZ();
				c_diff = c_mz + origin;

				//Attention! The +1. has nothing to do with the charge, it is caused by the wavelet's formula (tz1).
				current = c_diff > T_boundary_left && c_diff <= T_boundary_right ? IsotopeWavelet::getValueByLambda (my_local_lambda, c_diff*charge+1.)*c_ref[current_conv_pos].getIntensity() : 0;

				value += 0.5*(current + old)*(c_mz-old_pos);

				old = current;
				old_pos = c_mz;
			};		



			c_trans[my_local_pos].setIntensity (value);
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getTransformHighRes (MSSpectrum<PeakType>& c_trans, const MSSpectrum<PeakType>& c_ref, const UInt c)
	{
		Int spec_size ((Int)c_ref.size());
		//in the very unlikely case that size_t will not fit to int anymore this will be a problem of course
		//for the sake of simplicity (we need here a signed int) we do not cast at every following comparison individually
		UInt charge = c+1;
		DoubleReal value, T_boundary_left, T_boundary_right, c_diff, current, my_local_MZ, my_local_lambda, origin, c_mz;

		for(Int my_local_pos=0; my_local_pos<spec_size; ++my_local_pos)
		{
			value=0; T_boundary_left = 0, T_boundary_right = IsotopeWavelet::getMzPeakCutOffAtMonoPos(c_ref[my_local_pos].getMZ(), charge)/(DoubleReal)charge;


			my_local_MZ = c_ref[my_local_pos].getMZ(); my_local_lambda = IsotopeWavelet::getLambdaL(my_local_MZ*charge);
			c_diff=0;
			origin = -my_local_MZ+Constants::IW_QUARTER_NEUTRON_MASS/(DoubleReal)charge;

			for (Int current_conv_pos =  std::max(0, my_local_pos-from_max_to_left_); c_diff < T_boundary_right; ++current_conv_pos)
			{
				if (current_conv_pos >= spec_size)
				{
					break;
				};

				c_mz = c_ref[current_conv_pos].getMZ();
				c_diff = c_mz + origin;

				//Attention! The +1. has nothing to do with the charge, it is caused by the wavelet's formula (tz1).
				current = c_diff > T_boundary_left && c_diff <= T_boundary_right ? IsotopeWavelet::getValueByLambda (my_local_lambda, c_diff*charge+1.)*c_ref[current_conv_pos].getIntensity() : 0;
			


				value += current;
			};		

			c_trans[my_local_pos].setIntensity (value);
		};
	}

	
	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::initializeScan (const MSSpectrum<PeakType>& c_ref) 
	{
		computeMinSpacing(c_ref);

		Int wavelet_length = (UInt) ceil(max_mz_cutoff_/min_spacing_);
		Int max_index = (UInt) (Constants::IW_QUARTER_NEUTRON_MASS/min_spacing_);
		from_max_to_left_ = max_index;
		from_max_to_right_ = wavelet_length-1-from_max_to_left_;
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::computeMinSpacing (const MSSpectrum<PeakType>& c_ref) 
	{
		min_spacing_ = INT_MAX;
		for (UInt c_conv_pos=1; c_conv_pos<c_ref.size(); ++c_conv_pos) 
		{
			min_spacing_ = std::min (min_spacing_, c_ref[c_conv_pos].getMZ()-c_ref[c_conv_pos-1].getMZ()); 
		};
	}


	#ifdef OPENMS_HAS_CUDA
		template <typename PeakType>
		void IsotopeWaveletTransform<PeakType>::finalizeScanCuda ()
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
		int IsotopeWaveletTransform<PeakType>::initializeScanCuda (const MSSpectrum<PeakType>& scan, const UInt c) 
		{
			data_length_ = scan.size();		

			std::vector<float> pre_positions (data_length_), pre_intensities (data_length_);
			float c_spacing;
			min_spacing_=INT_MAX;
			pre_positions[0] = scan[0].getMZ();
			pre_intensities[0] = scan[0].getIntensity();
			for (UInt i=1; i<data_length_; ++i)
			{
				pre_positions[i] = scan[i].getMZ();
				c_spacing = pre_positions[i]-pre_positions[i-1];
				if (c_spacing < min_spacing_)
				{
					min_spacing_ = c_spacing;
				};	
				pre_intensities[i] = scan[i].getIntensity();
			};
			if (min_spacing_ == INT_MAX) //spectrum consists of a single data point
			{
				std::cout << "Scan consits of a single point. Unable to compute transform." << std::endl;
				return (Constants::CUDA_INIT_FAIL);
			};
			
			//Estimating the wave_length ...
			UInt wavelet_length=0, quarter_length=0;
			if(hr_data_) //We have to check this separately, because the simply estimation for LowRes data is destroyed by large gaps
			{
				UInt c_mz_cutoff;
				typename MSSpectrum<PeakType>::const_iterator start_iter, end_iter;
				for (UInt i=0; i<data_length_; ++i)
				{
					c_mz_cutoff =  IsotopeWavelet::getMzPeakCutOffAtMonoPos(scan[i].getMZ(), c+1);
					start_iter = scan.MZEnd(scan[i].getMZ());
					end_iter = scan.MZBegin(scan[i].getMZ()+c_mz_cutoff);
					wavelet_length = std::max ((long int) wavelet_length, distance(start_iter, end_iter)+1);
					end_iter = scan.MZEnd(scan[i].getMZ()-Constants::IW_QUARTER_NEUTRON_MASS/DoubleReal(c+1.));
					quarter_length = std::max ((long int) quarter_length, distance(end_iter, start_iter)+1);
				}
			}
			else
			{
				wavelet_length = (UInt) ceil(max_mz_cutoff_/min_spacing_);
			};
			//... done

			if (wavelet_length > data_length_ || wavelet_length == 1) //==1, because of 'ceil'
			{
				std::cout << "Warning: wavelet is larger than the number of data points. This might (!) severely affect the transform."<< std::endl;
				std::cout << "Warning/Error generated at scan with RT " << scan.getRT() << "." << std::endl;
			}; 

			UInt max_index;
			if (hr_data_)
			{
				max_index = quarter_length;
			}
			else
			{
				max_index = (UInt) (Constants::IW_QUARTER_NEUTRON_MASS/min_spacing_);
			}
			from_max_to_left_ = max_index;
			from_max_to_right_ = wavelet_length-1-from_max_to_left_;

			Int problem_size = Constants::CUDA_BLOCK_SIZE_MAX;
			to_load_ = problem_size + from_max_to_left_ + from_max_to_right_;	
			
			UInt missing_points = problem_size - (data_length_ % problem_size);
			overall_size_ = wavelet_length-1+data_length_+missing_points;

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

			cuda_intensities_.resize (overall_size_, 0); cuda_positions_.resize (overall_size_, 0);
			//Pad the values to the left; the positions should not matter if the values are zero
			float first_pos = pre_positions[0];
			for (Int i=0; i<from_max_to_left_; ++i)
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

			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET				
				std::stringstream name; name << "cuda_input_" << scan.getRT() << ".out\0"; 
				std::fstream outfile(name.str().c_str(), std::ios::out);
				for (size_t i=0; i<overall_size_; ++i)
					outfile << cuda_positions_[i] << " " << cuda_intensities_[i] << std::endl;
				outfile.close();
			#endif

			
			dimBlock_ = dim3 (Constants::CUDA_BLOCK_SIZE_MAX);
			to_compute_ = problem_size;
			dimGrid_ = dim3 ((data_length_+missing_points)/problem_size);

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
		void IsotopeWaveletTransform<PeakType>::getTransformCuda (TransSpectrum &c_trans, const UInt c) 
		{
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET				
				std::cout << "res in vector" << std::endl;
				std::vector<float> res (overall_size_, 0);
			#endif
			(cudaMemcpy(cuda_device_trans_intens_, &zeros_[0], overall_size_*sizeof(float), cudaMemcpyHostToDevice));	
			
			(cudaMemcpy(cuda_device_fwd2_, &zeros_[0], overall_size_*sizeof(float), cudaMemcpyHostToDevice));				
			getExternalCudaTransforms (dimGrid_, dimBlock_, (float*)cuda_device_pos_, (float*)cuda_device_intens_, from_max_to_left_, from_max_to_right_, (float*)cuda_device_trans_intens_, 
				c+1, to_load_, to_compute_, data_length_, (float*)cuda_device_fwd2_, hr_data_);			
			
			(cudaMemcpy(cuda_device_trans_intens_sorted_, cuda_device_fwd2_, overall_size_*sizeof(float), cudaMemcpyDeviceToDevice));
			
			(cudaMemcpy(&((*c_trans.trans_intens_)[0]), (float*)cuda_device_trans_intens_+from_max_to_left_, data_length_*sizeof(float), cudaMemcpyDeviceToHost));
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				(cudaMemcpy(&(res[0]), (float*)cuda_device_trans_intens_+from_max_to_left_, data_length_*sizeof(float), cudaMemcpyDeviceToHost));
				for (UInt i=0; i<data_length_; ++i)
				{
					c_trans.setTransIntensity (i, res[i]);
				};
			#endif
		}


		template <typename PeakType>
		int IsotopeWaveletTransform<PeakType>::sortCuda (MSSpectrum<PeakType>& sorted)
		{
			(cudaMemcpy(cuda_device_posindices_sorted_, &indices_[0], overall_size_*sizeof(int), cudaMemcpyHostToDevice));
			Int gpu_index = sortOnDevice((float*)cuda_device_trans_intens_sorted_, (int*) cuda_device_posindices_sorted_, overall_size_, 0);
			
			if (gpu_index < 0) //i.e., there is no positive intensity value at all
			{
				return (gpu_index);
			};
			
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
		void IsotopeWaveletTransform<PeakType>::identifyChargeCuda (const TransSpectrum& candidates,
			const UInt scan_index, const UInt c, const DoubleReal ampl_cutoff, const bool check_PPMs)
		{
			const MSSpectrum<PeakType>& ref (*candidates.getRefSpectrum()); 
			UInt index, MZ_start, MZ_end;
			typename MSSpectrum<PeakType>::iterator iter, bound_iter;
			typename MSSpectrum<PeakType>::const_iterator iter_start, iter_end, iter_p, iter2, seed_iter;
			DoubleReal mz_cutoff, seed_mz, c_av_intens=0, c_score=0, c_sd_intens=0, threshold=0, help_mz;
						
			Int gpu_index = sortCuda (c_sorted_candidate_), c_index;
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
		
			scoreOnDevice ((int*)cuda_device_posindices_sorted_, (float*)cuda_device_trans_intens_,  (float*)cuda_device_pos_, (float*)cuda_device_scores_, 
				c, num_of_scores, overall_size_, max_num_peaks_per_pattern_, threshold);
			
			(cudaMemcpy(&scores_[0], cuda_device_scores_, num_of_scores*sizeof(float), cudaMemcpyDeviceToHost));
	
			std::vector<float>::iterator score_iter;
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				std::stringstream stream;
				stream << "sorted_gpu_" << candidates.getRT() << "_" << c+1 << ".trans\0"; 
				std::ofstream ofile (stream.str().c_str());
				for (c_index = overall_size_-gpu_index-1, score_iter = scores_.begin()+num_of_scores-1; c_index >= 0; --c_index, --score_iter)
				{
					ofile << c_sorted_candidate_[c_index].getMZ() << "\t" << c_sorted_candidate_[c_index].getIntensity() << "\t" << *score_iter << std::endl;
				};
				ofile.close();
			#endif

			for (c_index = overall_size_-gpu_index-1, score_iter = scores_.begin()+num_of_scores-1; c_index >= 0; --c_index, --score_iter)
			{				
				seed_mz = c_sorted_candidate_[c_index].getMZ();
			
				//We can replace the following two lines ...
				//seed_iter = ref.MZBegin(seed_mz);
				//index = distance(ref.begin(), seed_iter);
				//... with:		
				index = h_pos_[c_index]-from_max_to_left_;
				seed_iter = ref.begin()+index;

				if (seed_iter == ref.end() || processed[distance(ref.begin(), seed_iter)] || index <= 0)
				{				
					continue;
				};
				
				mz_cutoff = IsotopeWavelet::getMzPeakCutOffAtMonoPos(seed_mz, c+1);
				//Mark the region as processed
				//Do not move this further down, since we have to mark this as processed in any case, 
				//even when score <=0; otherwise we would look around the maximum's position unless 
				//any significant point is found
				iter_start =ref.MZBegin(ref.begin(), seed_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.), seed_iter);
				iter_end = ref.MZEnd(seed_iter, seed_mz+mz_cutoff/(c+1.), ref.end());
			
				if (iter_end == ref.end())
				{
					--iter_end;
				};
			
				MZ_start = distance (ref.begin(), iter_start);
				MZ_end = distance (ref.begin(), iter_end);

				memset(&(processed[MZ_start]), 1, sizeof(UInt)*(MZ_end-MZ_start+1));

				c_score = *score_iter;
					
				if (c_score <= 0)
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
						UInt dist = distance(candidates.begin(), iter2);
						push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, 
							getLinearInterpolation ((iter2-1)->getMZ(), candidates.getTransIntensity(dist-1), help_mz, iter2->getMZ(), candidates.getTransIntensity(dist)),
							candidates.getRT(), MZ_start, MZ_end);
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
						UInt dist = distance(candidates.begin(), iter2);
						push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, 
							getLinearInterpolation ((iter2-1)->getMZ(), candidates.getTransIntensity(dist-1), help_mz, iter2->getMZ(), candidates.getTransIntensity(dist)),
							candidates.getRT(), MZ_start, MZ_end);
					};
				};
			};
			
			clusterSeeds_(candidates, ref, scan_index, c, check_PPMs);
		}
	#endif


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::identifyCharge (const MSSpectrum<PeakType>& candidates,
		const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt c, const DoubleReal ampl_cutoff, const bool check_PPMs)
	{
		Size scan_size (candidates.size()); 
		typename ConstRefVector<MSSpectrum<PeakType> >::iterator iter;
		typename MSSpectrum<PeakType>::const_iterator iter_start, iter_end, iter_p, seed_iter, iter2;
		DoubleReal mz_cutoff, seed_mz, c_av_intens=0, c_score=0, c_sd_intens=0, threshold=0, help_mz, share, share_pos, bwd, fwd;
		UInt MZ_start, MZ_end;
		
		MSSpectrum<PeakType> diffed (candidates);
		diffed[0].setIntensity(0); diffed[scan_size-1].setIntensity(0);

		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			std::stringstream stream;
			stream << "diffed_" << ref.getRT() << "_" << c+1 << ".trans\0"; 
			std::ofstream ofile (stream.str().c_str());
		#endif

		if (!hr_data_)//LowRes data
		{
			for (UInt i=0; i<scan_size-2; ++i)
			{
				share = candidates[i+1].getIntensity(), share_pos = candidates[i+1].getMZ();
				bwd = (share-candidates[i].getIntensity())/(share_pos-candidates[i].getMZ());
				fwd = (candidates[i+2].getIntensity()-share)/(candidates[i+2].getMZ()-share_pos);
				
				if (!(bwd>=0 && fwd<=0) || share > ref[i+1].getIntensity())
				{
					diffed[i+1].setIntensity(0);
				};
				
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
					ofile << diffed[i+1].getMZ() << "\t" <<  diffed[i+1].getIntensity() << std::endl;
				#endif
			};
		}
		else//HighRes data
		{
			for (UInt i=0; i<scan_size-2; ++i)
			{
				share = candidates[i+1].getIntensity(), share_pos = candidates[i+1].getMZ();
				bwd = (share-candidates[i].getIntensity())/(share_pos-candidates[i].getMZ());
				fwd = (candidates[i+2].getIntensity()-share)/(candidates[i+2].getMZ()-share_pos);
				
				if (!(bwd>=0 && fwd<=0))
				{
					diffed[i+1].setIntensity(0);
				};
				
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
					ofile << diffed[i+1].getMZ() << "\t" <<  diffed[i+1].getIntensity() << std::endl;
				#endif
			};
		}	
		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			ofile.close();
		#endif

		ConstRefVector<MSSpectrum<PeakType> > c_sorted_candidate_ (diffed.begin(), diffed.end());

		//Sort the transform in descending order according to the intensities present in the transform 	
		c_sorted_candidate_.sortByIntensity();
	
		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			std::stringstream stream2;
			stream2 << "sorted_cpu_" << candidates.getRT() << "_" << c+1 << ".trans\0"; 
			std::ofstream ofile2 (stream2.str().c_str());
			for (iter=c_sorted_candidate_.end()-1; iter != c_sorted_candidate_.begin(); --iter)
			{
				ofile2 << iter->getMZ() << "\t" << iter->getIntensity() << std::endl;
			};
			ofile2.close();
		#endif

		std::vector<UInt> processed (scan_size, 0);

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
			 
			mz_cutoff = IsotopeWavelet::getMzPeakCutOffAtMonoPos (seed_mz, c+1);
			//Mark the region as processed
			//Do not move this further down, since we have to mark this as processed in any case, 
			//even when score <=0; otherwise we would look around the maximum's position unless 
			//any significant point is found
			iter_start = ref.MZBegin(ref.begin(), seed_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.), seed_iter);
			iter_end = ref.MZEnd(seed_iter, seed_mz+mz_cutoff/(c+1.), ref.end());
			if (iter_end == ref.end())
			{
				--iter_end;
			};
							
			MZ_start = distance (ref.begin(), iter_start);
			MZ_end = distance (ref.begin(), iter_end);
	
			memset (&(processed[MZ_start]), 1, sizeof(UInt)*(MZ_end-MZ_start+1));

			c_score = scoreThis_ (candidates, IsotopeWavelet::getNumPeakCutOff(seed_mz*(c+1.)), seed_mz, c, threshold);
		
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
				if (trunc(seed_mz*10) == 5373)
					std::cout << seed_mz << "\t" << c_score << std::endl;
			#endif
	
			if (c_score <= 0)
			{
				continue;
			};

			//Push the seed into its corresponding box (or create a new one, if necessary)
			//Do ***NOT*** move this further down!
			push2TmpBox_ (seed_mz, scan_index, c, c_score, iter->getIntensity(), ref.getRT(), MZ_start, MZ_end);

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
					push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, getLinearInterpolation(iter2-1, help_mz, iter2), ref.getRT(), MZ_start, MZ_end);
				};
			};
		
			help_mz = seed_mz + Constants::IW_NEUTRON_MASS/(c+1.);
			iter2 = candidates.MZBegin (help_mz);
			if (iter2 == candidates.end()|| iter2 == candidates.begin())
			{
				continue;
			};

			if (fabs(iter2->getMZ()-seed_mz) > 0.5*Constants::IW_NEUTRON_MASS/(c+1.))
			//In the other case, we are too close to the peak, leading to incorrect derivatives.
			{
				if (iter2 != candidates.end())
				{
					push2TmpBox_ (iter2->getMZ(), scan_index, c, 0, getLinearInterpolation(iter2-1, help_mz, iter2), ref.getRT(), MZ_start, MZ_end);
				};
			};
		};	

		clusterSeeds_(candidates, ref, scan_index, c, check_PPMs);
	}

	#if defined(OPENMS_HAS_TBB) && defined (OPENMS_HAS_CUDA)
		template <typename PeakType>
		void IsotopeWaveletTransform<PeakType>::mergeFeatures (IsotopeWaveletTransform<PeakType>* later_iwt, const UInt RT_interleave, const UInt RT_votes_cutoff)
		{
			typename std::multimap<DoubleReal, Box>::iterator front_iter, end_iter, best_match, help_iter;

			//First of all do the trivial part of the merge
			for (end_iter=later_iwt->closed_boxes_.begin(); end_iter!=later_iwt->closed_boxes_.end(); ++end_iter)
			{
				closed_boxes_.insert (*end_iter);
			};

			typename std::multimap<DoubleReal, Box>& end_container (this->end_boxes_);
			typename std::multimap<DoubleReal, Box>& front_container (later_iwt->front_boxes_);
	
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
				std::cout << "FontBox: " << front_container.size() << std::endl;
				for (front_iter=front_container.begin(); front_iter != front_container.end(); ++front_iter)
					std::cout << front_iter->first << "\t" << front_iter->second.size() << std::endl;

				std::cout << "EndBox: " << end_container.size() << std::endl;
				for (front_iter=end_container.begin(); front_iter != end_container.end(); ++front_iter)
					std::cout << front_iter->first << "\t" << front_iter->second.size() << std::endl;
			#endif	

			typename std::multimap<UInt, BoxElement>::iterator biter;

			DoubleReal best_dist, c_dist; UInt c;
			//Now, try to find matching boxes for the rest
			for (front_iter=front_container.begin(); front_iter != front_container.end(); )
			{
				best_match = end_container.end(); best_dist = INT_MAX;
				//This is everything else than efficient, but both containers should be very small in size
				for (end_iter=end_container.begin(); end_iter != end_container.end(); ++end_iter)
				{
					c=0;
					for (biter=front_iter->second.begin(); biter != front_iter->second.end(); ++biter)
					{
						c=std::max (c, biter->second.c);
					}; 
					
					#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
						std::cout << "Trying to match: " << end_iter->first << " to front " <<  front_iter->first << std::endl;
					#endif

					c_dist = fabs(end_iter->first - front_iter->first); 
					if (c_dist < Constants::IW_HALF_NEUTRON_MASS/(c+1.) && c_dist < best_dist)
					{
						#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
							std::cout << "best match " << front_iter->second.begin()->first << "\t" << (--(end_iter->second.end()))->first << std::endl;
						#endif
						if ((front_iter->second.begin()->first - (--(end_iter->second.end()))->first) <= RT_interleave+1)
						//otherwise, there are too many blank scans in between
						{	
							best_match = end_iter;
							best_dist = c_dist;
						};
					};
				};
				if (best_match == end_container.end()) //No matching pair found
				{
					if (front_iter->second.size() >= RT_votes_cutoff)
					{
						closed_boxes_.insert (*front_iter);
						//extendBox_ (map, front_iter->second);
					};
					++front_iter;
				}		
				else //That's the funny part
				{
					#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
						std::cout << "Merging the boxes: " << front_iter->first  << "\t" << best_match->first << std::endl;
					#endif

					front_iter->second.insert (best_match->second.begin(), best_match->second.end());
					Box replacement (front_iter->second);	

					//We cannot divide both m/z by 2, since we already inserted some m/zs whose weight would be lowered.
					DoubleReal c_mz = front_iter->first * (front_iter->second.size()-best_match->second.size()) + best_match->first*best_match->second.size();	
					c_mz /= ((DoubleReal) (front_iter->second.size()));		

					help_iter = front_iter;
					++help_iter;
					std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help3 (c_mz, replacement);
					closed_boxes_.insert (help3);
					//extendBox_ (map, help3.second);				
					front_container.erase (front_iter);
					end_container.erase (best_match);
					front_iter = help_iter;
				};
			};

			//Merge the rest in end_container
			for (end_iter=end_container.begin(); end_iter != end_container.end(); ++end_iter)
			{
				if (end_iter->second.size() >= RT_votes_cutoff)
				{
					closed_boxes_.insert (*end_iter);
					//extendBox_ (map, end_iter->second);
				};
			};
		}
	#endif


	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThis_ (const MSSpectrum<PeakType>& candidate,
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal ampl_cutoff)
	{
		DoubleReal c_score=0, c_val;
		typename MSSpectrum<PeakType>::const_iterator c_left_iter2, c_right_iter2;
		Int signal_size ((Int)candidate.size());
		//in the very unlikely case that size_t will not fit to int anymore this will be a problem of course
		//for the sake of simplicity (we need here a signed int) we do not cast at every following comparison individually

		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff-1) -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> positions (end);
		for (Int i=0; i<end; ++i)
		{
			positions[i] =  seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-(i+1)*Constants::IW_HALF_NEUTRON_MASS)/((DoubleReal)c+1);
		};

		DoubleReal l_score=0, mid_val=0;
		Int start_index = distance(candidate.begin(), candidate.MZBegin (positions[0]))-1;
		for (Int v=1; v<=end; ++v, ++p_h_ind)
		{		
			do 
			{
				if (start_index < signal_size-1) ++start_index;
				else break;
			} while (candidate[start_index].getMZ() < positions[v-1]);
			
			if (start_index <= 0 || start_index>=signal_size-1) //unable to interpolate
			{
				continue;
			};

			c_left_iter2 = candidate.begin()+start_index-1;
			c_right_iter2 = c_left_iter2+1;

			c_val = c_left_iter2->getIntensity() + (c_right_iter2->getIntensity() - c_left_iter2->getIntensity())/(c_right_iter2->getMZ() - c_left_iter2->getMZ()) * (positions[v-1]-c_left_iter2->getMZ()); 	
			
			if (v == (int)(ceil(end/2.)))
			{
				l_score = c_score;
				mid_val = c_val;
			};

			if (p_h_ind%2 == 1) //I.e. a whole
			{
				c_score -= c_val;
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
					if (trunc(seed_mz*10) == 4897) std::cout << -c_val << std::endl;	
				#endif
			}
			else
			{
				c_score +=c_val;
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
					if (trunc(seed_mz*10) == 4897) std::cout << c_val << std::endl;	
				#endif
			};


			start_index = distance(candidate.begin(), c_left_iter2);
		};
	
		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
			std::ofstream ofile_score ("scores.dat", ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", ios::app);
			ofile_score.close();
			ofile_check_score.close();
		#endif

		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
			if (trunc(seed_mz*10) == 4897)
				std::cout << "final_score: " <<  seed_mz << "\t" << c_score << "\t l_score: " << l_score << "\t" << c_score-l_score-mid_val << "\t" <<  c_score-mid_val << "\t" << ampl_cutoff << std::endl;
		#endif

		if (l_score <=0 || c_score-l_score-mid_val <= 0 || c_score-mid_val <= ampl_cutoff)
		{
			return(0);
		};

		return (c_score);
	}


	
	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThis_ (const TransSpectrum& candidate, 
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal ampl_cutoff) 
	{
		DoubleReal c_score=0, c_val;
		typename MSSpectrum<PeakType>::const_iterator c_left_iter2, c_right_iter2;
		Int signal_size((Int)candidate.size());

		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff-1) -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> positions (end);
		for (Int i=0; i<end; ++i)
		{
			positions[i] =  seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-(i+1)*Constants::IW_HALF_NEUTRON_MASS)/((DoubleReal)c+1);
		};

		DoubleReal l_score=0, mid_val=0;
		Int start_index = distance(candidate.begin(), candidate.MZBegin (positions[0]))-1;
		for (Int v=1; v<=end; ++v, ++p_h_ind)
		{		
			do 
			{
				if (start_index < signal_size-1) ++start_index;
				else break;
			} while (candidate.getMZ(start_index) < positions[v-1]);
			
			if (start_index <= 0 || start_index>=signal_size-1) //unable to interpolate
			{
				continue;
			};

			c_left_iter2 = candidate.begin()+start_index-1;
			c_right_iter2 = c_left_iter2+1;

			c_val = candidate.getTransIntensity(start_index-1) + (candidate.getTransIntensity(start_index)-candidate.getTransIntensity(start_index-1))/(c_right_iter2->getMZ()-c_left_iter2->getMZ()) * (positions[v-1]-c_left_iter2->getMZ()); 
			if (v == (int)(ceil(end/2.)))
			{
				l_score = c_score;
				mid_val = c_val;
			};

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
	
		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
			std::ofstream ofile_score ("scores.dat", ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", ios::app);
			ofile_score << c_check_point << "\t" << c_score << std::endl;
			ofile_score.close();
			ofile_check_score.close();
		#endif

		if (l_score <=0 || c_score-l_score-mid_val <= 0 || c_score-mid_val <= ampl_cutoff)
		{
			return(0);
		};

		return (c_score);
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::clusterSeeds_ (const MSSpectrum<PeakType>& candidate, 
		const MSSpectrum<PeakType>& ref,  const UInt scan_index, const UInt c, const bool check_PPMs) 
	{
		typename std::multimap<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
		std::vector<BoxElement> final_box;
	 	DoubleReal c_mz, av_score=0, av_mz=0, av_intens=0, av_abs_intens=0, count=0;
		DoubleReal virtual_av_score=0, virtual_av_mz=0, virtual_av_intens=0, virtual_av_abs_intens=0, virtual_count=0;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=tmp_boxes_->at(c).begin(); iter!=tmp_boxes_->at(c).end(); ++iter)
		{	

			Box& c_box = iter->second;
			av_score=0, av_mz=0, av_intens=0, av_abs_intens=0, count=0;
			virtual_av_score=0, virtual_av_mz=0, virtual_av_intens=0, virtual_av_abs_intens=0, virtual_count=0;

			//Now, let's get the RT boundaries for the box
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				if (box_iter->second.score == 0) //virtual helping point
				{
					if (count != 0) continue; //it is in any way not pure virtual

					c_mz = box_iter->second.mz;
					virtual_av_intens += box_iter->second.intens;
					virtual_av_abs_intens += fabs(box_iter->second.intens);
					virtual_av_mz += c_mz*fabs(box_iter->second.intens);
					++virtual_count;
				}
				else
				{
					c_mz = box_iter->second.mz;
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_abs_intens += fabs(box_iter->second.intens);
					av_mz += c_mz*fabs(box_iter->second.intens);
					++count;
				};
			};

			if (count == 0) //pure virtual helping box
			{
				av_intens = virtual_av_intens/virtual_count;		
				av_score = 0; 
				av_mz = virtual_av_mz/virtual_av_abs_intens;
			}
			else
			{
				av_intens /= count;		
				av_score /= count; 
				av_mz /= av_abs_intens;
			};

			BoxElement c_box_element;
			c_box_element.mz = av_mz;
			c_box_element.c = c;
			c_box_element.score = av_score;
			c_box_element.intens = av_intens;
		
			c_box_element.RT=c_box.begin()->second.RT;
			final_box.push_back(c_box_element);
		};	

		Size num_o_feature = final_box.size();				
		if (num_o_feature == 0)
		{
			tmp_boxes_->at(c).clear();
			return;
		};

		//Computing the derivatives
		std::vector<DoubleReal> bwd_diffs(num_o_feature, 0);

		bwd_diffs[0]=0;
		for (Size i=1; i<num_o_feature; ++i)
		{
			bwd_diffs[i] = (final_box[i].intens-final_box[i-1].intens)/(final_box[i].mz-final_box[i-1].mz);
		};		

		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			std::ofstream ofile_bwd ("bwd_cpu.dat");
			for (Size i=0; i<num_o_feature; ++i)
			{
				ofile_bwd << final_box[i].mz << "\t" << bwd_diffs[i] << std::endl;
			};
			ofile_bwd.close();
		#endif
			

		for (Size i=0; i<num_o_feature-1; ++i)
		{	
			while (i<num_o_feature-2)
			{
				if(final_box[i].score>0) //this has been an helping point
					break;
				++i;
			};
	
			if (bwd_diffs[i]>0 && bwd_diffs[i+1]<0)
			{	
				checkPositionForPlausibility_ (candidate, ref, final_box[i].mz, final_box[i].c, scan_index, check_PPMs, final_box[i].intens);	
				continue;
			};
		};

		tmp_boxes_->at(c).clear();
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvIntens_ (const MSSpectrum<PeakType>& scan) 
	{ 
		DoubleReal av_intens=0;
		for (UInt i=0; i<scan.size(); ++i)
		{
			if (scan[i].getIntensity() >= 0)
			{
				av_intens += scan[i].getIntensity();
			}
		};
		return (av_intens / (double)scan.size());
	}	


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) 
	{
		DoubleReal res=0, intens;
		for (UInt i=0; i<scan.size(); ++i)
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
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvIntens_ (const TransSpectrum& scan) 
	{ 
		DoubleReal av_intens=0;
		for (UInt i=0; i<scan.size(); ++i)
		{
			if (scan.getTransIntensity(i) >= 0)
			{
				av_intens += scan.getTransIntensity(i);
			}
		};
		return (av_intens / (double)scan.size());
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getSdIntens_ (const TransSpectrum& scan, const DoubleReal mean) 
	{
		DoubleReal res=0, intens;
		for (UInt i=0; i<scan.size(); ++i)
		{		
			if (scan.getTransIntensity(i) >= 0)
			{
				intens = scan.getTransIntensity(i);
				res += (intens-mean)*(intens-mean);
			};
		};
		return (sqrt(res/(double)(scan.size()-1)));
	}

	

	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::push2Box_ (const DoubleReal mz, const UInt scan, UInt c, 
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end, DoubleReal ref_intens)
	{
		const DoubleReal dist_constraint (Constants::IW_HALF_NEUTRON_MASS/(DoubleReal)max_charge_);

		typename std::multimap<DoubleReal, Box>::iterator upper_iter (open_boxes_.upper_bound(mz));
		typename std::multimap<DoubleReal, Box>::iterator lower_iter ( open_boxes_.lower_bound(mz));

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
				if (fabs((--lower_iter)->first - mz) < dist_constraint) //matching box
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
			if (upper_iter == open_boxes_.end() && fabs(lower_iter->first - mz) < dist_constraint) //Found matching Box
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
			dist_lower = (dist_lower < dist_constraint) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < dist_constraint) ? dist_upper : INT_MAX;

			if (dist_lower>=dist_constraint && dist_upper>=dist_constraint) // they are both too far away
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
		element.c = c; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens; element.ref_intens=ref_intens;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end;


		if (create_new_box == false)
		{
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
		const DoubleReal dist_constraint (Constants::IW_HALF_NEUTRON_MASS/(DoubleReal)max_charge_);

		std::multimap<DoubleReal, Box>& tmp_box (tmp_boxes_->at(c));
		typename std::multimap<DoubleReal, Box>::iterator upper_iter (tmp_box.upper_bound(mz));
		typename std::multimap<DoubleReal, Box>::iterator lower_iter (tmp_box.lower_bound(mz)); 
	
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
				if (fabs((--lower_iter)->first - mz) < dist_constraint) //matching box
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
			if (upper_iter == tmp_box.end() && fabs(lower_iter->first - mz) < dist_constraint) //Found matching Box
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
			dist_lower = (dist_lower < dist_constraint) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < dist_constraint) ? dist_upper : INT_MAX;

			if (dist_lower>=dist_constraint && dist_upper>=dist_constraint) // they are both too far away
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
		element.c = c; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens; element.ref_intens=-1000;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end;

		if (create_new_box == false)
		{
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
			std::pair<UInt, BoxElement> help2 (scan, element);
			std::multimap<UInt, BoxElement> help3;
			help3.insert (help2);

			std::pair<DoubleReal, std::multimap<UInt, BoxElement> > help4 (mz, help3);
			tmp_box.insert (help4);
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::updateBoxStates (const MSExperiment<PeakType>& map, const Size scan_index, const UInt RT_interleave,
		const UInt RT_votes_cutoff, const Int front_bound, const Int end_bound)
	{
		typename std::multimap<DoubleReal, Box>::iterator iter, iter2;
	
		if ((Int)scan_index == end_bound && end_bound != (Int)map.size()-1)
		{
			for (iter=open_boxes_.begin(); iter!=open_boxes_.end(); ++iter)
			{	
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET	
					std::cout << "LOW THREAD insert in end_box " << iter->first << std::endl;
					typename Box::iterator dings;
					for (dings=iter->second.begin(); dings != iter->second.end(); ++dings)
						std::cout << map[dings->first].getRT() << "\t" << dings->second.c+1 <<  std::endl;  
				#endif
				end_boxes_.insert (*iter);
			};
			open_boxes_.clear();
			return;
		};

		for (iter=open_boxes_.begin(); iter!=open_boxes_.end(); )
		{
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				if (front_bound > 0)
				{
					std::cout << "HIGH THREAD open box. " << iter->first << "\t current scan index " << scan_index << "\t" << ((iter->second.begin()))->first << "\t of last scan " << map.size()-1 << "\t" << front_bound << std::endl; 
				};
			#endif

			//For each Box we need to figure out, if and when the last RT value has been inserted
			UInt lastScan = (--(iter->second.end()))->first;
			if (scan_index - lastScan > RT_interleave+1 || scan_index == map.size()-1) //I.e. close the box!
			{
				if (iter->second.begin()->first -front_bound <= RT_interleave+1 && front_bound > 0)
				{
					iter2=iter;
					++iter2;	
					#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
						std::cout << "HIGH THREAD insert in front_box " << iter->first << std::endl;
					#endif
					front_boxes_.insert (*iter);
					open_boxes_.erase (iter);
					iter=iter2;
					continue;
				};

				iter2 = iter;
				++iter2;
				//Please do **NOT** simplify the upcoming lines.
				//The 'obvious' overhead is necessary since the object represented by iter might be erased
				//by push2Box which might be called by extendBox_.
				if (iter->second.size() >= RT_votes_cutoff)
				{
					//extendBox_ (map, iter->second);
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
		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
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
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
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
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				std::cout << box.begin()->second.RT << "\t" << av_mz << "\t" << box.begin()->second.c+1 << "\t" << " extending the box " << std::endl;
			#endif

			push2Box_ (av_mz, c_index, box.begin()->second.c, box.begin()->second.score, c_iter->getIntensity(),
				map[c_index].getRT(), box.begin()->second.MZ_begin, box.begin()->second.MZ_end);
		};
	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::clusterSeeds_ (const TransSpectrum& candidates, 
		const MSSpectrum<PeakType>& ref,  const UInt scan_index, const UInt c, const bool check_PPMs) 
	{
		typename std::multimap<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
		std::vector<BoxElement> final_box;
	 	DoubleReal c_mz, av_score=0, av_mz=0, av_intens=0, av_abs_intens=0, count=0;
		DoubleReal virtual_av_score=0, virtual_av_mz=0, virtual_av_intens=0, virtual_av_abs_intens=0, virtual_count=0;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=tmp_boxes_->at(c).begin(); iter!=tmp_boxes_->at(c).end(); ++iter)
		{	
			Box& c_box = iter->second;
			av_score=0, av_mz=0, av_intens=0, av_abs_intens=0, count=0;
			virtual_av_score=0, virtual_av_mz=0, virtual_av_intens=0, virtual_av_abs_intens=0, virtual_count=0;

			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				if (box_iter->second.score == 0) //virtual helping point
				{
					if (count != 0) continue; //it is in any way not pure virtual

					c_mz = box_iter->second.mz;
					virtual_av_intens += box_iter->second.intens;
					virtual_av_abs_intens += fabs(box_iter->second.intens);
					virtual_av_mz += c_mz*fabs(box_iter->second.intens);
					++virtual_count;
				}
				else
				{
					c_mz = box_iter->second.mz;
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_abs_intens += fabs(box_iter->second.intens);
					av_mz += c_mz*fabs(box_iter->second.intens);
					++count;
				};
			};

			if (count == 0) //pure virtual helping box
			{
				av_intens = virtual_av_intens/virtual_count;		
				av_score = 0; 
				av_mz = virtual_av_mz/virtual_av_abs_intens;
			}
			else
			{
				av_intens /= count;		
				av_score /= count; 
				av_mz /= av_abs_intens;
			};
				
			BoxElement c_box_element;
			c_box_element.mz = av_mz;
			c_box_element.c = c;
			c_box_element.score = av_score;
			c_box_element.intens = av_intens;
		
			c_box_element.RT=c_box.begin()->second.RT;

			final_box.push_back(c_box_element);	
		};	

		UInt num_o_feature = final_box.size();			
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
			bwd_diffs[i] = (final_box[i].intens-final_box[i-1].intens)/(final_box[i].mz-final_box[i-1].mz);
		};		

		#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
			std::ofstream ofile_bwd ("bwd_gpu.dat");
			for (UInt i=0; i<num_o_feature; ++i)
			{
				ofile_bwd << final_box[i].mz << "\t" << bwd_diffs[i] << std::endl;
			};
			ofile_bwd.close();
		#endif
			
		for (UInt i=0; i<num_o_feature-1; ++i)
		{	
			while (i<num_o_feature-2)
			{
				if(final_box[i].score>0) //this has been an helping point
					break;
				++i;
			};
			
			if (bwd_diffs[i]>0 && bwd_diffs[i+1]<0)
			{					
				checkPositionForPlausibility_ (candidates, ref, final_box[i].mz, final_box[i].c, scan_index, check_PPMs, final_box[i].intens);	
				continue;
			};
		};

		tmp_boxes_->at(c).clear();
	}



	template <typename PeakType>
	FeatureMap<Feature> IsotopeWaveletTransform<PeakType>::mapSeeds2Features (const MSExperiment<PeakType>& map, const UInt RT_votes_cutoff) 
	{
		FeatureMap<Feature> feature_map;
		typename std::multimap<DoubleReal, Box>::iterator iter;
		typename Box::iterator box_iter;
		UInt best_charge_index; DoubleReal best_charge_score, c_mz, c_RT; UInt c_charge;
		DoubleReal av_intens=0, av_ref_intens=0, av_score=0, av_mz=0, av_RT=0, mz_cutoff;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=closed_boxes_.begin(); iter!=closed_boxes_.end(); ++iter)
		{
			Box& c_box = iter->second;
			std::vector<DoubleReal> charge_votes (max_charge_, 0), charge_binary_votes (max_charge_, 0);

			//Let's first determine the charge
			//Therefor, we can use two types of votes: qualitative ones (charge_binary_votes) or quantitative ones (charge_votes)
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
					if (trunc(box_iter->second.mz) == 369)
						std::cout << box_iter->second.c << "\t" <<  box_iter->second.intens << std::endl;
				#endif

				charge_votes[box_iter->second.c] += box_iter->second.intens;//score; Do not use score, can get problematic for charge state 2 vs 4
				++charge_binary_votes[box_iter->second.c];
			};

			//... determining the best fitting charge
			best_charge_index=0; best_charge_score=0;
			for (UInt i=0; i<max_charge_; ++i)
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
			};

			c_charge = best_charge_index + 1; //that's the finally predicted charge state for the pattern

			av_intens=0, av_ref_intens=0, av_score=0, av_mz=0, av_RT=0;
			//Now, let's get the RT boundaries for the box
			std::vector<DPosition<2> > point_set;
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				c_mz = box_iter->second.mz;
				c_RT = box_iter->second.RT;
   
				mz_cutoff = IsotopeWavelet::getMzPeakCutOffAtMonoPos (c_mz, c_charge);

				point_set.push_back (DPosition<2> (c_RT, c_mz - Constants::IW_QUARTER_NEUTRON_MASS/(DoubleReal)c_charge)); 
				//-1 since we are already at the first peak and +0.75, since this includes the last peak of the wavelet as a whole
				point_set.push_back (DPosition<2> (c_RT, c_mz + mz_cutoff/(DoubleReal)c_charge)); 
				if (best_charge_index == box_iter->second.c)
				{
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_ref_intens += box_iter->second.ref_intens;
					av_mz += c_mz*box_iter->second.intens;
				};
				av_RT += c_RT;
			};

			av_mz /= av_intens;
			av_ref_intens /= (DoubleReal)charge_binary_votes[best_charge_index];
			av_score /= (DoubleReal)charge_binary_votes[best_charge_index];
			av_RT /= (DoubleReal)c_box.size();
			
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				if (trunc(av_mz) == 369)
					std::cout << av_mz << "\t" << best_charge_index << "\t" << best_charge_score << std::endl;
			#endif

			Feature c_feature;
			ConvexHull2D c_conv_hull;
			c_conv_hull.addPoints (point_set);
			c_feature.setCharge (c_charge);
			c_feature.setConvexHulls (std::vector<ConvexHull2D> (1, c_conv_hull));

			//This makes the intensity value independent of the m/z (the lambda) value (Skellam distribution)
			DoubleReal lambda = IsotopeWavelet::getLambdaL(av_mz*c_charge);
			av_intens /= exp(-2*lambda) * boost::math::cyl_bessel_i(0, 2*lambda);
			//DoubleReal breen (0.000594*av_mz*c_charge-0.03091);
			//av_ref_intens /= exp(-breen);
	
			c_feature.setMZ (av_mz);
			c_feature.setIntensity (av_intens);
			c_feature.setRT (av_RT);
			c_feature.setOverallQuality (av_score);
			feature_map.push_back (c_feature);
		};

		return (feature_map);
	}


	template <typename PeakType>
	bool IsotopeWaveletTransform<PeakType>::checkPositionForPlausibility_ (const MSSpectrum<PeakType>& candidate,
		const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, const UInt c, const UInt scan_index, const bool check_PPMs, const DoubleReal transintens)
	{		
		typename MSSpectrum<PeakType>::const_iterator iter, ref_iter; 
		UInt peak_cutoff;
		peak_cutoff = IsotopeWavelet::getNumPeakCutOff (seed_mz, c+1);

		iter = candidate.MZBegin(seed_mz);
		//we can ignore those cases
		if (iter==candidate.begin() || iter==candidate.end()) 
		{
			return (false);
		};
			
		std::pair<DoubleReal, DoubleReal> reals;
		ref_iter =  ref.MZBegin(seed_mz);
		//Correct the position
		DoubleReal real_mz, real_intens;
		if (check_PPMs)
		{
			reals = checkPPMTheoModel_ (ref, iter->getMZ(), c);
			real_mz = reals.first, real_intens = reals.second;
			if (real_mz <= 0 || real_intens <= 0)
			{	
				typename MSSpectrum<PeakType>::const_iterator h_iter=ref_iter, hc_iter=ref_iter;
				while (h_iter != ref.begin())
				{	 
					--h_iter;
					if (h_iter->getIntensity() > hc_iter->getIntensity() || (h_iter->getIntensity() == hc_iter->getIntensity() && hc_iter->getIntensity() == 0))
					{
						hc_iter=h_iter;
					}
					else
					{
						break;
					};
			
					if (seed_mz - h_iter->getMZ() > Constants::IW_QUARTER_NEUTRON_MASS/(c+1.))
					{
						return (false);
					}
				}
				reals = checkPPMTheoModel_ (ref, h_iter->getMZ(), c);
				DoubleReal real_mz = reals.first, real_intens = reals.second;
				if (real_mz <= 0 || real_intens <= 0)
				{
					return (false);
				}
				real_mz = h_iter->getMZ();
				real_intens = h_iter->getIntensity();
			}
		}
		else
		{
			reals = std::pair<DoubleReal, DoubleReal> (seed_mz, ref_iter->getIntensity());
			real_mz = reals.first, real_intens = reals.second;

			if (real_mz <= 0 || real_intens <= 0)
			{
				typename MSSpectrum<PeakType>::const_iterator h_iter=ref_iter, hc_iter=ref_iter;
				while (h_iter != ref.begin())
				{	 
					--h_iter;
					if (h_iter->getIntensity() > hc_iter->getIntensity() || (h_iter->getIntensity() == hc_iter->getIntensity() && hc_iter->getIntensity() == 0))
					{
						hc_iter=h_iter;
					}
					else
					{
						break;
					};
			
					if (seed_mz - h_iter->getMZ() > Constants::IW_QUARTER_NEUTRON_MASS/(c+1.))
					{
						return (false);
					}
				}
				real_mz = h_iter->getMZ(), real_intens = h_iter->getIntensity();
				if (real_mz <= 0 || real_intens <= 0)
				{
					return (false);
				}
				real_mz = h_iter->getMZ();
				real_intens = h_iter->getIntensity();
			} 
		};

		DoubleReal c_score = scoreThis_ (candidate, peak_cutoff, real_mz, c, 0);

		if (c_score <= 0)
		{
			return (false);
		};

		DoubleReal mz_cutoff = IsotopeWavelet::getMzPeakCutOffAtMonoPos (real_mz, c+1);
		typename MSSpectrum<PeakType>::const_iterator real_l_MZ_iter = ref.MZBegin(real_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));		
		typename MSSpectrum<PeakType>::const_iterator real_r_MZ_iter = ref.MZBegin(real_l_MZ_iter, real_mz+mz_cutoff/(c+1.), ref.end());
		if (real_r_MZ_iter == ref.end())
		{
			--real_r_MZ_iter;
		};
	

		UInt real_mz_begin = distance (ref.begin(), real_l_MZ_iter);
		UInt real_mz_end = distance (ref.begin(), real_r_MZ_iter);

		push2Box_ (real_mz, scan_index, c, c_score, transintens, ref.getRT(), real_mz_begin, real_mz_end, real_intens);
		return (true);
	}

	template <typename PeakType>
	bool IsotopeWaveletTransform<PeakType>::checkPositionForPlausibility_ (const TransSpectrum& candidate,
		const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, const UInt c, const UInt scan_index, const bool check_PPMs, const DoubleReal transintens)
	{
		typename MSSpectrum<PeakType>::const_iterator iter, ref_iter; 
		UInt peak_cutoff;
		peak_cutoff = IsotopeWavelet::getNumPeakCutOff (seed_mz, c+1);

		iter = candidate.MZBegin(seed_mz);
		//we can ignore those cases
		if (iter==candidate.begin() || iter==candidate.end()) 
		{		
			return (false);
		};
			
		std::pair<DoubleReal, DoubleReal> reals;
		ref_iter =  ref.MZBegin(seed_mz);
		//Correct the position
		DoubleReal real_mz, real_intens;
		if (check_PPMs)
		{
			reals = checkPPMTheoModel_ (ref, iter->getMZ(), c);
			real_mz = reals.first, real_intens = reals.second;
			if (real_mz <= 0 || real_intens <= 0)
			{	
				typename MSSpectrum<PeakType>::const_iterator h_iter=ref_iter, hc_iter=ref_iter;
				while (h_iter != ref.begin())
				{	 
					--h_iter;
					if (h_iter->getIntensity() > hc_iter->getIntensity() || (h_iter->getIntensity() == hc_iter->getIntensity() && hc_iter->getIntensity() == 0))
					{
						hc_iter=h_iter;
					}
					else
					{
						break;
					};
			
					if (seed_mz - h_iter->getMZ() > Constants::IW_QUARTER_NEUTRON_MASS/(c+1.))
					{
						return (false);
					}
				}
				reals = checkPPMTheoModel_ (ref, h_iter->getMZ(), c);
				real_mz = reals.first, real_intens = reals.second;
				if (real_mz <= 0 || real_intens <= 0)
				{
					return (false);
				}
				real_mz = h_iter->getMZ();
				real_intens = h_iter->getIntensity();
			}
		}
		else
		{
			reals = std::pair<DoubleReal, DoubleReal> (seed_mz, ref_iter->getIntensity());
			real_mz = reals.first, real_intens = reals.second;
					
			if (real_mz <= 0 || real_intens <= 0)
			{
				typename MSSpectrum<PeakType>::const_iterator h_iter=ref_iter, hc_iter=ref_iter;
				while (h_iter != ref.begin())
				{	 
					--h_iter;
					if (h_iter->getIntensity() > hc_iter->getIntensity() || (h_iter->getIntensity() == hc_iter->getIntensity() && hc_iter->getIntensity() == 0))
					{
						hc_iter=h_iter;
					}
					else
					{
						break;
					};
			
					if (seed_mz - h_iter->getMZ() > Constants::IW_QUARTER_NEUTRON_MASS/(c+1.))
					{
						return (false);
					}
				}
				real_mz = h_iter->getMZ(), real_intens = h_iter->getIntensity();
				if (real_mz <= 0 || real_intens <= 0)
				{
					return (false);
				}
				real_mz = h_iter->getMZ();
				real_intens = h_iter->getIntensity();
			} 
		};
		
		DoubleReal c_score = scoreThis_ (candidate, peak_cutoff, real_mz, c, 0);

		if (c_score <= 0)
		{
			return (false);
		};

		DoubleReal mz_cutoff = IsotopeWavelet::getMzPeakCutOffAtMonoPos (real_mz, c+1);
		typename MSSpectrum<PeakType>::const_iterator real_l_MZ_iter = ref.MZBegin(real_mz-Constants::IW_QUARTER_NEUTRON_MASS/(c+1.));		
		typename MSSpectrum<PeakType>::const_iterator real_r_MZ_iter = ref.MZBegin(real_l_MZ_iter, real_mz+mz_cutoff/(c+1.), ref.end());
		if (real_r_MZ_iter == ref.end())
		{
			--real_r_MZ_iter;
		};
	

		UInt real_mz_begin = distance (ref.begin(), real_l_MZ_iter);
		UInt real_mz_end = distance (ref.begin(), real_r_MZ_iter);

		push2Box_ (real_mz, scan_index, c, c_score, transintens, ref.getRT(), real_mz_begin, real_mz_end, real_intens);
		return (true);
	}


	template <typename PeakType>
	std::pair<DoubleReal, DoubleReal> IsotopeWaveletTransform<PeakType>::checkPPMTheoModel_ (const MSSpectrum<PeakType>& ref, const DoubleReal c_mz, const UInt c)
	{
		DoubleReal mass = c_mz*(c+1)-Constants::IW_PROTON_MASS*(c);
		DoubleReal ppms = getPPMs(peptideMassRule(mass), mass);
		if (ppms >= Constants::PEPTIDE_MASS_RULE_THEO_PPM_BOUND)
		{
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << "ISO_WAVE" << ")> " << "REJECT \t" << ppms 
				<< " (rule: " << peptideMassRule(mass) << " got: " << mass << ")" << std::endl;
			return (std::pair<DoubleReal, DoubleReal> (-1,-1));
		};
			std::cout << ::std::setprecision(8) << std::fixed << c_mz << "\t =(" << "ISO_WAVE" << ")> " << "ACCEPT \t" << ppms << " (rule: " 
				<< peptideMassRule(mass) << " got: " << mass << ")" << std::endl;
		return (std::pair<DoubleReal, DoubleReal> (c_mz, ref.MZBegin(c_mz)->getIntensity()));	
	}

} //namespace

#endif
