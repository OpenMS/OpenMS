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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_TRANSFORM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_TRANSFORM_H


#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/Averagine.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <map>
#include <backward/hash_map.h>
#include <sstream>

#ifndef DEFAULT_HASH_PRECISION //setting the default value
#define DEFAULT_HASH_PRECISION 1000 //i.e. float precision (w.r.t. the m/z coordindate) is 3
#endif

#ifndef DEFAULT_NUM_OF_INTERPOLATION_POINTS
#define DEFAULT_NUM_OF_INTERPOLATION_POINTS 3 
#endif

#ifndef DEBUG_MALDI
#define DEBUG_MALDI 1
#endif


namespace OpenMS
{
	/** @brief A class implementing the isotope wavelet transform. 
 		* If you just want to find features using the isotope wavelet, take a look at the IsotopeWaveletFF class. Usually, you only
 		* have to consider this class if you plan to change the basic implementation of the transform. 
 		*
 		*	@ingroup FeatureFinder */ 
	template <typename PeakType>
	class IsotopeWaveletTransform
	{
		public:
			
			/** @brief Default Constructor. */
			IsotopeWaveletTransform () throw();

			/** @brief Destructor. */
			virtual ~IsotopeWaveletTransform () throw();
		
			/** @brief Computes the discrete-time continuous wavelet transform simultaneously for several charges.
 				* 
 				* The function computes the isotope wavelet transformed versions of @p scan. 
 				* The transform is determined for several charge states (up to charge @p max_charge) at the same time.
 				* Hence, the user has to ensure that the size of @p transforms equals to @p max_charge and that each spectrum in 
 				* @p transforms has the same length as @p scan.
 				*
 				* @param scan The MS scan you wish to transform.
 				* @param transforms A vector (with indices running from 0 to @p max_charge-1) of MS spectra (each of the size of @p scan).
 				* The code will NOT check the allocated memory (the sizes) for transforms and its entries.  
 				* @param max_charge The maximal charge state that is considered.
 				* @return The transformed spectra. Entry i in this vector corresponds to the "i+1"-charge-state-transform 
 				* of @p scan. */ 	
			static void getTransforms (const MSSpectrum<PeakType>& scan, 
				std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge) throw ();


			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
 				* process is stored internally. Important: Before calling this function, apply updateRanges() to the original map.
 				*
 				* @param candidates A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
 				* charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function getTransforms.)    
 				* @param scan The index of the scan (w.r.t. to some map) currently under consideration. 
 				* @param ampl_cutoff The thresholding parameter. This parameter is the only (and hence a really important)
 				* parameter of the isotope wavelet transform. On the basis of @p ampl_cutoff the program tries to distinguish between 
 				* and noise signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
 				* line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
 				* get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
 				* makes sense to start with @p ampl_cutoff=0. */ 
			virtual void identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates,
				const UInt scan_index, const DoubleReal ampl_cutoff=0) throw ();
			
			
			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score 
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param start_index One of the indices enclosing the interesting region.
 				* @param end_index One of the indices enclosing the interesting region.
 				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined. 
 				* @param intens The intensity of the transform at @p seed_mz.
 				* @param The threshold. */
			virtual DoubleReal scoreThis (const MSSpectrum<PeakType>& candidate, 
				const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff=0) throw ();
			

			/** @brief Returns a floating point precision, which is internally used for a hash map. 
				*	The default value is 1000, which corresponds to a floating point precision of 3. Usually, this value can be left unchanged. 		
				*	If you are considering some very high resulted spectra (e.g. FT-ICR) you might increase this value to 4 or 5. The parameter
				*	is used to mark regions of interest for which an isotopic pattern has already been found. Hence, these regions will not 
				*	be considered during the rest of the feature finding algorithm. */
			inline DoubleReal getHashPrecision () const throw ()
			{ 
				return (hash_precision_); 
			}

			/** @brief See getHashPrecision. */
			void setHashPrecision (const DoubleReal hash_precision) throw ()
				{ hash_precision_ = hash_precision; }
			
			/** @brief A function keeping track of currently open and closed sweep line boxes. 
 				* This function is used by the isotope wavelet feature finder and must be called for each processed scan. 
 				* @param c_scan_number The index of the scan currently under consideration w.r.t. its MS map. 
 				* This information is necessary to sweep across the map affter each scan has been evaluated. 
 				* @param See the IsotopeWaveletFF class. 
 				* @param The sweep line cutoff parameter. See also IsotopeWaveletFF class.*/	
			void updateBoxStates (const UInt c_scan_number, const UInt RT_interleave, 
					const UInt RT_votes_cutoff) throw ();

			/** @brief Maps the internally used data structures to the OpenMS framework. 
 				* @param max_charge The maximal charge state under consideration. 
 				* @param The sweep line cutoff parameter. See also IsotopeWaveletFF class.*/
			FeatureMap<Feature> mapSeeds2Features (const UInt max_charge, const UInt RT_votes_cutoff) throw ();
	

		protected:			

			struct BoxElement_
			{			
				DoubleReal mz;
				UInt c; //Note, this is not the charge (it is charge-1!!!)
				DoubleReal score;
				DoubleReal intens;
				DoubleReal RT; //The elution time (not the scan index)
			};

			typedef std::map<UInt, BoxElement_> Box_; ///<Key: RT index, value: BoxElement_

			/** @brief Samples the wavelet at discrete time points, s.t. they match automatically to the m/z positions provided
 				* in @p scan and returns the discrete values of psi in @psi.  
 				*	@param mz_index The start index of @p scan for which the the wavelet should be adapted.  
 				* @param offset Is the offset the wavelet function needs to be aligned with a signal point.
 				* @param z The charge the wavelet function should adapt (corresponds to z in the paper).
 				* @param av_MZ_spacing The average spacing in the m/z domain.
 				* @param psi The sampled values.
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */ 
			static void sampleTheWavelet_ (const MSSpectrum<PeakType>& scan, const UInt mz_index, 
				const DoubleReal offset, const UInt z, const DoubleReal av_MZ_spacing, std::vector<DoubleReal>& psi, 
				const Int mode=+1) throw ();		

			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getAvintens_ (const MSSpectrum<PeakType>& scan) throw (); 		
			
			/** @brief Computes the standard deviation (neglecting negative values) of the intensity of @p scan. */
			inline DoubleReal getSdintens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) throw ();

			/** @brief A wrapper function to the GSL interpolation routine. */
			DoubleReal getCubicInterpolatedValue_ (const std::vector<DoubleReal>& x, const DoubleReal xi, const std::vector<DoubleReal>& y) throw ();

			/** @brief A function to map m/z values to m/z indices. In particular useful if you know already the
 				* approximate position of the corresponding entry which can be indicated by @p start. */ 		
			inline std::pair<int, int> getNearBys_ (const MSSpectrum<PeakType>& signal, const DoubleReal mz, 
				const UInt start=0) const throw ();

			/** @brief Inserts a potential isotopic pattern into an open box or - if no such box exists - creates a new one.
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the postprocessing (sweep lining). 
 				* @param The estimated charge state of the pattern. 
 				* @param The pattern's score. 
 				* @param The intensity at the monoisotopic peak. 
 				* @param The retention time of the scan (similar to @p scan, but here: no index, but the real value). */
			virtual void push2Box_ (const DoubleReal mz, const UInt scan, UInt charge, const DoubleReal score, 
				const DoubleReal intens, const DoubleReal rt) throw ();

			/** @brief Computes the average MZ spacing of @p scan in the range @p start_index to @p end_index. */
			static DoubleReal getAvMZSpacing_ (const MSSpectrum<PeakType>& scan, Int start_index=0, Int end_index=-1) throw ();


			//A floating point precision for an internally used hash map
			DoubleReal hash_precision_;	

			//internally used data structures for the sweep line algorithm	
			std::map<DoubleReal, Box_> open_boxes_, closed_boxes_;	//DoubleReal = average m/z position
		
			gsl_interp_accel* acc_;
			gsl_spline* spline_; 
	};







	template <typename PeakType>	
	bool comparator (const PeakType& a, const PeakType& b)
	{
		return (a.getIntensity() > b.getIntensity());
	}		


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform () throw() : hash_precision_ (DEFAULT_HASH_PRECISION)
	{	
		acc_ = gsl_interp_accel_alloc ();
		spline_ = gsl_spline_alloc (gsl_interp_cspline, DEFAULT_NUM_OF_INTERPOLATION_POINTS); 
	}


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::~IsotopeWaveletTransform () throw()
	{
		gsl_interp_accel_free (acc_);
		gsl_spline_free (spline_);
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getTransforms (const MSSpectrum<PeakType>& scan, 
		std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge) throw ()
	{	
		UInt scan_size = scan.size();
		DoubleReal av_MZ_spacing = getAvMZSpacing_(scan);
		UInt peak_cutoff = IsotopeWavelet::getPeakCutoff();
		UInt wavelet_length = (UInt) trunc(peak_cutoff/av_MZ_spacing);	
		std::vector<DoubleReal> psi (wavelet_length, 0); //The wavelet
		
		if (scan_size < wavelet_length)
		{
			return;
		};
		
		DoubleReal cum_spacing, c_spacing, //Helping variables
			max_w_monoi_intens=QUARTER_NEUTRON_MASS, //The position of the monoisotopic peak within the coordinate sys. of the wavelet 
			sums=0, //Helping variables
			max_position_scan=0, //The position of the data point (within the scan) we want to align with
			align_offset, //Correction term; shifts the wavelet to get the desired alignment
			last;
		UInt c=0, k=0, j=0;
		DoubleReal c_charge; //DoubleReal, since we will oven devide by c_charge 
		
		//The upcoming variable is necessary to capture strange effects in special types of unequally spaced data sets.
		//Imagine some wholes in the m/z range (points the mass spectrometer did not sample). If they become larger than 
		//0.25*NEUTRON_MASS (considering the case of charge 1), several data points wll share the same max_position, 
		//causing the upcoming code to crash since suddenly some m/z pom/z positions willl occure twice. The interval of multiple 
		//occuring points is stored by multiple_s and implicitly by i.
		std::vector<int> multiple_s (max_charge,-1);
		std::vector<DoubleReal> last_max_position_scan (max_charge, -1);
		bool repair=false;

		//Starting convolution
		for (UInt i=0; i<scan_size; ++i)
		{
			//Now, let's sample the wavelets
			for (c=0; c<max_charge; ++c)
			{	
				c_charge=c+1;
				cum_spacing=0;				
				max_w_monoi_intens=QUARTER_NEUTRON_MASS/c_charge; //This is the position of the monoistopic peak (centered)
				
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
						cum_spacing += av_MZ_spacing;
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
				//since we will get trouble by the NEUTRON_MASS and the resulting numerical instabilities. 
				//We will add this correcting term at the end of the whole processing.
				max_position_scan = scan[(i+j)%scan_size].getMZ();
				if (max_position_scan == last_max_position_scan[c]) //Uuups, multiple times the same m/z coordinate
				{
					if (multiple_s[c] < 0) //This is the first entry where this artifact occured
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
					
				//Sampling the wavelet 
				sampleTheWavelet_ (scan, i, cum_spacing, (UInt) c_charge, av_MZ_spacing, psi);
					
				//The convolution
				k=0; sums=0;
				for (UInt j=i; j<scan_size && k<wavelet_length; ++j, ++k)
				{
					sums += scan[j].getIntensity()*psi[k];
				};

				if (k< wavelet_length) // I.e. we have an overlapping wavelet
				{
					sums=0; // => We can absolutely neglect this feature since it is too near at the boundary.
					max_position_scan = transforms[c][i-1].getMZ()+av_MZ_spacing;
				};

				//Store the current convolution result
				PeakType& c_peak1 (transforms[c][i]);
				c_peak1.setIntensity(sums);
				c_peak1.setMZ(max_position_scan);
				//transforms[c][i].setIntensity(sums);
				//transforms[c][i].setMZ(max_position_scan);	

				if (repair)
				{		
					//std::cout << "repair" << std::endl;
					UInt noi2interpol = i - multiple_s[c]; //NOT +1

					//The special case if we are the boundary (exactly the last point in the spectrum)
					if (i == scan_size-1)
					{
						//We do not care about the intenities, since we will set them to zero anyway.
						//We would just like to avoid multiple positions to occur in the transform
						for (UInt ii=0; ii<=noi2interpol; ++ii) 
							//it must be "<=noi..." !!! not "<", since in this case we do not want to keep the last multiple	
							//the same holds for "ii=0"
						{
							transforms[c][multiple_s[c]+ii].setMZ(transforms[c][multiple_s[c]-1].getMZ() + (ii+1)*av_MZ_spacing);		
						};

						last_max_position_scan[c] = max_position_scan; //Reset
						multiple_s[c]=-1; //Reset
						repair=false;
						continue;
					}

					PeakType& c_peak2 (transforms[c][multiple_s[c]]);
					DoubleReal x1 = c_peak2.getMZ(); //transforms[c][multiple_s[c]].getMZ();
					DoubleReal y1 = c_peak2.getMZ(); transforms[c][multiple_s[c]].getIntensity();					
					DoubleReal x2 = c_peak1.getMZ();//transforms[c][i].getMZ();
					//if (i >= scan_size) //Is still possible and ugly => reset x2
					//{
					//	x2 = 2*c_peak1.getMZ()-transforms[c][i-1].getMZ(); //(transforms[c][i].getMZ() - transforms[c][i-1].getMZ()) + transforms[c][i].getMZ();
					//};
					DoubleReal y2 = c_peak1.getIntensity();//transforms[c][i].getIntensity();
					//if (i >= scan_size)
					//{
					//	y2 = transforms[c][i].getIntensity(); //Do just anything, does not matter what, since we are at the boundary
					//};
					DoubleReal dx = (x2-x1)/(DoubleReal)(noi2interpol);
					for (UInt ii=1; ii<noi2interpol; ++ii) //ii=1, not 0, since we want to keep the first of the multiples
					{	
						//transforms[c][multiple_s[c]+ii].setMZ(transforms[c][multiple_s[c]].getMZ()+ii*dx);
						transforms[c][multiple_s[c]+ii].setMZ(c_peak2.getMZ()+ii*dx);
						//transforms[c][multiple_s[c]+ii].setIntensity(y1 + (y2-y1)/(x2-x1)*(transforms[c][multiple_s[c]].getMZ()+ii*dx-x1));
						transforms[c][multiple_s[c]+ii].setIntensity(y1 + (y2-y1)/(x2-x1)*(c_peak2.getMZ()+ii*dx-x1));
					};

					last_max_position_scan[c] = max_position_scan; //Reset
					multiple_s[c]=-1; //Reset
					repair=false;
				}				
			}
		}
		return;
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates, 
		const UInt scan_index, const DoubleReal ampl_cutoff) throw ()
	{
		DoubleReal av_MZ_spacing = getAvMZSpacing_(candidates[0]);
		UInt peak_cutoff = IsotopeWavelet::getPeakCutoff();
		UInt wavelet_length = (UInt) trunc(peak_cutoff/av_MZ_spacing);	
		UInt cands_size = candidates.size();
		UInt signal_size=candidates[0].size(), c_index, i_iter, end_index; 
		Int start_index; //Do not change this to UInt
		typename MSSpectrum<PeakType>::iterator iter;
		DoubleReal seed_mz, c_av_intens, c_score, c_sd_intens, threshold;
		std::pair<DoubleReal, DoubleReal> c_processed;

		//For all charges do ...
		for (UInt c=0; c<cands_size; ++c)		
		{
			//Indicates wheter for some specific region of the scan an isotopic pattern has already been identified
			//i.e.: In the moment, we do not care about overlapping signals! 
			std::vector<bool> processed = std::vector<bool> (signal_size, false); 
			MSSpectrum<PeakType> c_sorted_candidate (candidates[c]); 

			//The fowllong hash map allows a fast transform from m/z positions to m/z indices w.r.t. the transformed vector
			hash_multimap<int, UInt> index_hash;
			for (UInt i=0; i<signal_size; ++i)
			{
				Int hash_key = (int) trunc(candidates[c][i].getMZ()*hash_precision_); 
				index_hash.insert (std::pair<int, UInt> (hash_key, i));
			};

			//Sort the transform in descending order according to the intensities present in the transform 	
			sort (c_sorted_candidate.begin(), c_sorted_candidate.end(), comparator<PeakType>); 
			c_av_intens = getAvintens_ (candidates[c]);		
			c_sd_intens = getSdintens_ (candidates[c], c_av_intens);

			//Eliminate uninteresting regions
			//In principle we should do that in a binary fashion for efficiency reasons ...  
			for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter)
			{
				if (iter->getIntensity() <= c_av_intens)
				{	
					break;
				};
			};
			c_sorted_candidate.erase (iter, c_sorted_candidate.end());	

			if (ampl_cutoff < 0)
			{
				threshold = 0;
			}
			else
			{
				threshold=ampl_cutoff*c_sd_intens + c_av_intens;
			};

			

			i_iter=0;
			for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter, ++i_iter)
			{					
				seed_mz=iter->getMZ();
				c_index = index_hash.find((int)trunc(seed_mz*hash_precision_))->second;

				if (processed[c_index])
				{
					continue;
				};
					
				//In order to determine the start and end indices, we first need to know the width of the region one should consider 
				//to estimate the mean and the sd of the pattern candidate. 
				//That region is defined by the position of the heighst amplitude +/- wavelet_length_.			
				start_index = c_index - wavelet_length-1;
				end_index = c_index + wavelet_length+1;

				if (start_index < 0)
				{
					start_index = 0;
				};
				if (end_index >= signal_size)
				{
					end_index = signal_size-1;			
				};

				//Mark as processed
				for (UInt z=start_index; z<=end_index; ++z)
				{
					processed[z] = true;
				};	
				
				c_score = scoreThis (candidates[c], seed_mz, c, iter->getIntensity(), threshold);
		
				if (c_score <= 0)
				{
					continue;
				};

				//Push the seed into its corresponding box (or create a new one, if necessary)
				push2Box_ (seed_mz, scan_index, c, c_score, iter->getIntensity(), candidates[0].getRT());
			};	
		};
	}
		

	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::sampleTheWavelet_ (const MSSpectrum<PeakType>& scan, const UInt mz_index, 
		const DoubleReal offset, const UInt z, const DoubleReal av_MZ_spacing, std::vector<DoubleReal>& psi, const Int mode) throw ()
	{
		UInt scan_size = scan.size();
		DoubleReal c_pos, lambda;

		psi.resize (scan_size); //just to be sure; if psi is already scan_size large, this will is a simple test
		
		c_pos = scan[mz_index].getMZ();				
		lambda = IsotopeWavelet::getLambdaQ(c_pos*z-mode*z*PROTON_MASS);
		UInt peak_cutoff = IsotopeWavelet::getPeakCutoff();
		UInt wavelet_length = (UInt) trunc(peak_cutoff/av_MZ_spacing);	

		if (mz_index+wavelet_length >= scan_size)
		{
			psi = std::vector<double> (wavelet_length, 0);
			return;
		}

		DoubleReal cum_spacing=offset;
		std::vector<DoubleReal> c_mzs (wavelet_length+1), c_spacings (wavelet_length);
		c_mzs[0] = scan[mz_index].getMZ();
		for (UInt j=1; j<wavelet_length+1; ++j)
		{
			c_mzs[j] = scan[mz_index+j].getMZ();
			c_spacings[j-1] = c_mzs[j]-c_mzs[j-1];
			c_spacings[j-1] = (c_spacings[j-1] > 0) ? c_spacings[j-1] : av_MZ_spacing;
		}

		//Building up (sampling) the wavelet
		const double limit=QUARTER_NEUTRON_MASS+IsotopeWavelet::getPeakCutoff();
		double tz1;
		for (UInt j=0; j<wavelet_length && cum_spacing < limit; ++j)
		{
			tz1=cum_spacing*z+1;
			psi[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambda (lambda, tz1) : 0;
			cum_spacing += c_spacings[j];
		};

		/*if (trunc(c_mzs[0]) == 1000 || trunc(c_mzs[0]) == 1700 || trunc(c_mzs[0]) == 2000 || trunc(c_mzs[0]) == 3000)
		{
			std::stringstream stream; stream << c_mzs[0] << "_" << z << ".dat\0"; 
			std::ofstream ofile (stream.str().c_str());
			for (unsigned int zzz=0; zzz<wavelet_length; ++zzz)
				ofile << scan[mz_index+zzz].getMZ() << "\t" << psi[zzz] << std::endl;
			ofile.close();
		};*/

		//Normally, we should substract the mean by now, but since the effect is marginal on real world
		//data, we will skip this step which is mainly of theoretical interest.
	}


	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThis (const MSSpectrum<PeakType>& candidate, 
		const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) throw ()
	{	
		UInt peak_cutoff = IsotopeWavelet::getPeakCutoff();
		DoubleReal c_score=0, c_check_point=-1, c_val;
		std::pair<int, int> c_between_left, c_between_right; 				
		typename MSSpectrum<PeakType>::const_iterator c_left_iter, c_right_iter, c_iter;

		DoubleReal leftBound, rightBound;
		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*peak_cutoff -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> xvec, yvec, weights;
		std::vector<DoubleReal> xs (DEFAULT_NUM_OF_INTERPOLATION_POINTS), ys(DEFAULT_NUM_OF_INTERPOLATION_POINTS);
		UInt i=0; 
		
		for (Int v=1; v<end; ++v, ++p_h_ind)
		{		
			c_check_point = seed_mz-(peak_cutoff*NEUTRON_MASS-v*0.5*NEUTRON_MASS)/((DoubleReal)c+1);
			
			leftBound = c_check_point;
			rightBound = c_check_point;		

			c_left_iter = candidate.MZBegin (c_check_point);
			c_right_iter = c_left_iter;

			//ugly, but the only way to check it I guess
			if (c_left_iter == candidate.begin()) continue;		
			--c_left_iter;
			if (c_right_iter == candidate.end()) continue;
			if (++c_right_iter == candidate.end()) continue;
		
			for (c_iter=c_left_iter, i=0; c_iter!=c_right_iter; ++c_iter, ++i)
			{
				xs[i] = c_iter->getMZ();	
				ys[i] = c_iter->getIntensity();
			}
			xs[i] = c_iter->getMZ();	
			ys[i] = c_iter->getIntensity();

			c_val = getCubicInterpolatedValue_ (xs, c_check_point, ys);

			if (p_h_ind%2 == 1) //I.e. a whole
			{
				c_score -= c_val;
			}
			else
			{
				c_score +=c_val;
			};
		};

		if (c_score <= ampl_cutoff+intens)
		{
			return(0);
		};

		return (log(c_score)+ log(intens));
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvMZSpacing_ (const MSSpectrum<PeakType>& scan, Int start_index, Int end_index) throw ()
	{ 
		DoubleReal av_MZ_spacing=0;
		if (end_index < 0)
		{
			end_index = scan.size();
		};
		for (Int i=start_index; i<end_index-1; ++i)
		{
			 av_MZ_spacing += scan[i+1].getMZ() - scan[i].getMZ();
		};
		return (av_MZ_spacing / (DoubleReal) (end_index-1-start_index));
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvintens_ (const MSSpectrum<PeakType>& scan) throw ()
	{ 
		DoubleReal av_intens=0, count=0;
		for (UInt i=0; i<scan.size(); ++i)
		{
			if (scan[i].getIntensity() >= 0)
			{
				av_intens += scan[i].getIntensity(); 
				++count;
			}
			//av_intens += (scan[i].getIntensity() < 0) ? 0 : scan[i].getIntensity();
			//av_intens += fabs(scan[i].getIntensity());
		};
		//return (av_intens / (DoubleReal) (scan.size()));	
		return (av_intens / count);
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getSdintens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) throw ()
	{
		DoubleReal res=0, intens, count=0;
		for (UInt i=0; i<scan.size(); ++i)
		{	
			if (scan[i].getIntensity() >= 0)
			{
				intens = scan[i].getIntensity();
				res += (intens-mean)*(intens-mean);
				++count;
			};
			//intens = (scan[i].getIntensity() < 0) ? 0 : scan[i].getIntensity();
			//intens = fabs(scan[i].getIntensity());
			//res += (intens-mean)*(intens-mean);
		};

		//return (sqrt(res/(DoubleReal)(scan.size()-1)));		
		return (sqrt(res/(count-1)));
	}


	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::getCubicInterpolatedValue_ (const std::vector<DoubleReal>& x, 
		const DoubleReal xi, const std::vector<DoubleReal>& y) throw ()
	{
		gsl_spline_init (spline_, &x[0], &y[0], x.size());
		DoubleReal yi = gsl_spline_eval (spline_, xi, acc_);
		return (yi);	
	}


	template <typename PeakType>
	std::pair<int, int> IsotopeWaveletTransform<PeakType>::getNearBys_ (const MSSpectrum<PeakType>& signal, const DoubleReal mz, 
		const UInt start) const throw ()
	{
		for (UInt i=start; i<signal.size(); ++i)
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
	void IsotopeWaveletTransform<PeakType>::push2Box_ (const DoubleReal mz, const UInt scan, UInt charge, 
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt) throw ()
	{	
		typename std::map<DoubleReal, Box_>::iterator upper_iter = open_boxes_.upper_bound(mz);
		typename std::map<DoubleReal, Box_>::iterator lower_iter; 
		if (open_boxes_.empty())
		{
			lower_iter = open_boxes_.end();
		}
		else
		{
			lower_iter = open_boxes_.lower_bound(mz);
		};

		//Ugly, but necessary due to the implementation of STL lower_bound
		if (mz != open_boxes_.lower_bound(mz)->first && lower_iter != open_boxes_.begin())
		{
			lower_iter = --(open_boxes_.lower_bound(mz));
		};

		typename std::map<DoubleReal, Box_>::iterator insert_iter;
		bool createNewBox=false;
		if (lower_iter == open_boxes_.end()) //I.e. there is no open Box for that mz position
		{
			createNewBox=true;
		};

		if (upper_iter == open_boxes_.end() && fabs(lower_iter->first - mz) < 0.5*NEUTRON_MASS) //Found matching Box
		{
			insert_iter = lower_iter;
			createNewBox=false;
		}
		else
		{
			createNewBox=true;
		};

		if (upper_iter != open_boxes_.end() && lower_iter != open_boxes_.end())
		{	
			//Figure out which entry is closer to m/z
			DoubleReal dist_lower = fabs(lower_iter->first - mz);
			DoubleReal dist_upper = fabs(upper_iter->first - mz);
			dist_lower = (dist_lower < 0.5*NEUTRON_MASS) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < 0.5*NEUTRON_MASS) ? dist_upper : INT_MAX;

			if (dist_lower>=0.5*NEUTRON_MASS && dist_upper>=0.5*NEUTRON_MASS) // they are both too far away
			{
				createNewBox=true;
			}
			else
			{
				insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;	
				createNewBox=false;
			};
		}; 

		BoxElement_ element; element.c = charge; element.mz = mz; element.score = score; element.RT = rt, element.intens=intens;
		std::pair<UInt, BoxElement_> help2 (scan, element);

		if (createNewBox == false)
		{	
			insert_iter->second.insert (help2);	

			//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
			Box_ replacement (insert_iter->second);	

			//We cannot devide both m/z by 2, since we already inserted some m/z's whose weight would be lowered.
			//Also note that we alread inserted the new entry, leading to size-1.
			DoubleReal c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;	
			c_mz /= ((DoubleReal) insert_iter->second.size());			
			
			//Now let's remove the old and insert the new one
			open_boxes_.erase (insert_iter);	
			std::pair<DoubleReal, std::map<UInt, BoxElement_> > help3 (c_mz, replacement);	
			open_boxes_.insert (help3);		
		}
		else
		{
			std::map<UInt, BoxElement_> help3;
			help3.insert (help2);
			std::pair<DoubleReal, std::map<UInt, BoxElement_> > help4 (mz, help3);
			open_boxes_.insert (help4);
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::updateBoxStates (const UInt c_scan_number, const UInt RT_interleave, 
		const UInt RT_votes_cutoff) throw ()
	{
		typename std::map<DoubleReal, Box_>::iterator iter, iter2;
		for (iter=open_boxes_.begin(); iter!=open_boxes_.end(); )
		{
			//For each Box we need to figure out, if and when the last RT value has been inserted
			//If the Box his unchanged since RT_interleave_ scans, we will close the Box.
			UInt lastScan = (--(iter->second.end()))->first;
			if (c_scan_number - lastScan > RT_interleave) //I.e. close the box!
			{
				iter2 = iter;
				++iter2;
				if (iter->second.size() >= RT_votes_cutoff)
				{
					closed_boxes_.insert (*iter);
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
	FeatureMap<Feature> IsotopeWaveletTransform<PeakType>::mapSeeds2Features 
		(const UInt max_charge, const UInt RT_votes_cutoff) throw ()
	{
		FeatureMap<Feature> feature_map;
		typename std::map<DoubleReal, Box_>::iterator iter;
		typename Box_::iterator box_iter;
		UInt best_charge_index; DoubleReal best_charge_score, c_mz, c_RT; UInt c_charge; 		
		ConvexHull2D c_conv_hull;
		Feature c_feature;

		#ifdef DEBUG_MALDI
			std::ofstream ofile ("mascot.query");
		#endif

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=closed_boxes_.begin(); iter!=closed_boxes_.end(); ++iter)
		{		
			Box_& c_box = iter->second;
			std::vector<DoubleReal> charge_votes (max_charge, 0), charge_binary_votes (max_charge, 0);
		
			//Let's first determine the charge
			//Therefor, we can use two types of votes: qulitative ones (charge_binary_votes) or quantitaive ones (charge_votes)
			//Collting the votes ...
			for (box_iter=iter->second.begin(); box_iter!=iter->second.end(); ++box_iter)
			{
				charge_votes[box_iter->second.c] += box_iter->second.score;
				++charge_binary_votes[box_iter->second.c];
			};
			
			//... dertermining the best fitting charge
			best_charge_index=0; best_charge_score=0; 
			for (UInt i=0; i<max_charge; ++i)
			{
				if (charge_votes[i] > best_charge_score)
				{
					best_charge_index = i;
					best_charge_score = charge_votes[i];
				};	
			};		

			//Pattern found in too few RT scan 
			if (charge_binary_votes[best_charge_index] < RT_votes_cutoff)
			{
				continue;
			};

			c_charge = best_charge_index + 1; //that's the finally predicted charge state for the pattern

			DoubleReal av_intens=0, av_mz=0, av_RT=0;
			//Now, let's get the RT boundaries for the box
			std::vector<DPosition<2> > point_set;
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
			{
				c_mz = box_iter->second.mz;
				c_RT = box_iter->second.RT;
				Averagine::getModel (c_mz, c_charge, &c_extend);
				std::cout << box_iter->first << " (" << c_RT << ")\t" << c_extend.first << "\t" 
					<< c_extend.second << "\t" << box_iter->second.c +1 << "#" << c_charge << std::endl;
				point_set.push_back (DPosition<2> (c_RT, c_extend.first)); //mz start, RT
				point_set.push_back (DPosition<2> (c_RT, c_extend.second)); //mz start, RT
				av_intens += box_iter->second.intens;
				av_mz += c_mz;
				av_RT += c_RT;
			};
			av_mz /= (DoubleReal)c_box.size();
			av_intens /= (DoubleReal)c_box.size();
			av_RT /= (DoubleReal)c_box.size();
			std::cout << "**************************************************" << std::endl;

			c_conv_hull = point_set;
			c_feature.setCharge (c_charge);
			c_feature.setConvexHulls (std::vector<ConvexHull2D> (1, c_conv_hull));
			c_feature.setMZ (av_mz);
			c_feature.setIntensity (av_intens);
			c_feature.setRT (av_RT);
			feature_map.push_back (c_feature);		

			#ifdef DEBUG_MALDI
				ofile << av_mz << "\t" << av_intens << std::endl;
			#endif

		};		
	
		#ifdef DEBUG_MALDI
			ofile.close();
		#endif

		return (feature_map);
	}

} //namespace

#endif 
