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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_TRANSFORM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_TRANSFORM_H


#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics_double.h>
#include <vector>
#include <map>
#include <backward/hash_map.h>
#include <sstream>
#include <iomanip>

#ifndef DEFAULT_NUM_OF_INTERPOLATION_POINTS
#define DEFAULT_NUM_OF_INTERPOLATION_POINTS 3 
#endif

#ifndef EPSILON_ION_COUNTS
#define EPSILON_ION_COUNTS 3 
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
					
			struct BoxElement_
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

			typedef std::multimap<UInt, BoxElement_> Box_; ///<Key: RT index, value: BoxElement_	
			/** @brief Default Constructor. 
 				* @note Provided just for inheritence reasons. You should always use the other constructor. */
			

			IsotopeWaveletTransform () throw();

			/** @brief Constructor. */
			IsotopeWaveletTransform (const UInt max_charge, const bool create_Mascot_PMF_File) throw();

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
 				* @param mode The recording mode of the mass spectrometer (+1 or -1).
 				* @return The transformed spectra. Entry i in this vector corresponds to the "i+1"-charge-state-transform 
 				* of @p scan. */ 	
			static void getTransforms (const MSSpectrum<PeakType>& scan, 
				std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge, const Int mode) throw ();


			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
 				* process is stored internally. Important: Before calling this function, apply updateRanges() to the original map.
 				*
 				* @param candidates A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
 				* charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function getTransforms.)    
 				* @param scan_index The index of the scan (w.r.t. to some map) currently under consideration. 
 				* @param ampl_cutoff The thresholding parameter. This parameter is the only (and hence a really important)
 				* parameter of the isotope wavelet transform. On the basis of @p ampl_cutoff the program tries to distinguish between 
 				* noise and signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
 				* line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
 				* get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
 				* makes sense to start with something slighly above 0, e.g. @p ampl_cutoff=0.1 or even @p ampl_cutoff=-1, 
 				* indicating no thresholding at all. Note that also ampl_cutoff=0 triggers (a moderate) thresholding based on the 
 				* averagine intensity in the wavelet transform. */ 
			virtual void identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates,
				const MSSpectrum<PeakType>& ref, const UInt scan_index, const DoubleReal ampl_cutoff=0) throw ();
			
			
			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score 
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param start_index One of the indices enclosing the region of interest.
 				* @param end_index One of the indices enclosing the region of interest.
 				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined. 
 				* @param intens The intensity of the transform at @p seed_mz.
 				* @param ampl_cutoff The threshold. */
			virtual DoubleReal scoreThis (const MSSpectrum<PeakType>& candidate, const UInt peak_cutoff, 
				const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff=0) throw ();
			

			/** @brief A function keeping track of currently open and closed sweep line boxes. 
 				* This function is used by the isotope wavelet feature finder and must be called for each processed scan. 
 				* @param scan_index The index of the scan currently under consideration w.r.t. its MS map. 
 				* This information is necessary to sweep across the map affter each scan has been evaluated. 
 				* @param RT_interleave See the IsotopeWaveletFF class. 
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class. */	
			void updateBoxStates (const UInt scan_index, const UInt RT_interleave, 
				const UInt RT_votes_cutoff) throw ();


			/** @brief Filters the candidates further more and maps the internally used data structures to the OpenMS framework. 
 				* @map The original map containing the data set to be analyzed.
 				* @param max_charge The maximal charge state under consideration. 
 				* @param RT_votes_cutoff See the IsotopeWaveletFF class.*/
			FeatureMap<Feature> mapSeeds2Features (const MSExperiment<PeakType>& map, const UInt max_charge, const UInt RT_votes_cutoff) throw ();


			/** @brief A ugly but necessary function to handle "off-by-1-Dalton predictions" due to idiosyncrasies of the data set
 				* (in comparision to the averagine model)
 				* @param candidate The wavelet transformed spectrum containing the candidate. 
 				* @param ref The original spectrum containing the candidate.
 				* @param seed_mz The m/z position of the candidate pattern.
 				* @param c The predicted charge state of the candidate.
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map).
 				*
 				* @todo THIS FUNCTION ONLY WORKS FOR LOWER RANGE MASSES AT THE MOMENT (up to about 4000 Da). */ 
			virtual void checkPosition (const MSSpectrum<PeakType>& candidate, const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, 
				const UInt c, const UInt scan_index) throw ();


			/** @brief Returns the closed boxes. */
			virtual std::multimap<DoubleReal, Box_> getClosedBoxes () throw ()
				{ return (closed_boxes_); };	


		protected:			


			/** @brief Samples the wavelet at discrete time points, s.t. they match automatically the m/z positions provided
 				* in @p scan. Returns the discrete values of psi in @psi.  
 				*	@param mz_index The start index of @p scan for which the the wavelet should be adapted.  
 				* @param offset The offset the wavelet function needs to be aligned with a signal point.
 				* @param charge The charge (not the index c!) the wavelet function should adapt (corresponds to z in the paper).
 				* @param av_MZ_spacing The average spacing in the m/z domain.
 				* @param psi The sampled values.
 				* @param peak_cutoff The number of peaks we will consider for the isotopic pattern.
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */ 
			static void sampleTheWavelet_ (const MSSpectrum<PeakType>& scan, const UInt wavelet_length, const UInt mz_index, 
				const DoubleReal offset, const UInt charge, const DoubleReal av_MZ_spacing, std::vector<DoubleReal>& psi, const DoubleReal peak_cutoff, 
				const Int mode=+1) throw ();		

			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline DoubleReal getAvIntens_ (const MSSpectrum<PeakType>& scan) throw (); 		
			
			/** @brief Computes the standard deviation (neglecting negative values) of the intensity of @p scan. */
			inline DoubleReal getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) throw ();

			/** @brief A wrapper function to the GSL interpolation routine. */
			DoubleReal getCubicInterpolatedValue_ (const std::vector<DoubleReal>& x, const DoubleReal xi, const std::vector<DoubleReal>& y) throw ();

			/** @brief A function to map m/z values to m/z indices. In particular useful if you know already the
 				* approximate position of the corresponding entry which can be indicated by @p start. */ 		
			inline std::pair<Int, Int> getNearBys_ (const MSSpectrum<PeakType>& signal, const DoubleReal mz, 
				const UInt start=0) const throw ();

			/** @brief Inserts a potential isotopic pattern into an open box or - if no such box exists - creates a new one.
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the postprocessing (sweep lining). 
 				* @param charge The estimated charge state of the pattern. 
 				* @param score The pattern's score. 
 				* @param intens The intensity at the monoisotopic peak. 
 				* @param rt The retention time of the scan (similar to @p scan, but here: no index, but the real value). */
			virtual void push2Box_ (const DoubleReal mz, const UInt scan, UInt charge, const DoubleReal score, 
				const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) throw ();

			/** @brief Essentially the same function as push2Box_. 
 				* In contrast to push2Box this function stores its candidates only temporarily. In particular, this
 				* function is only used within a single scan transform. After the wavelet transform is computed on
 				* that scan, all candidates are pushed by this function and finally clustered together by clusterSeeds. 
 				* Afterwards, a final push by push2Box_ is performed storing the clustered candidates. */
			virtual void push2TmpBox_ (const DoubleReal mz, const UInt scan, UInt charge, const DoubleReal score, 
				const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) throw ();


			/** @brief Computes the average MZ spacing of @p scan in the range @p start_index to @p end_index. */
			static DoubleReal getAvMZSpacing_ (const MSSpectrum<PeakType>& scan, Int start_index=0, Int end_index=-1) throw ();
 
				
			/** @brief The trapezoid rule for integration. 
 				* @param a first x coordinate.
 				* @param b second x coordinate.
 				* @param fa a's corresponding function value.
 				* @param fb b's corresponding function value. */
			static DoubleReal chordTrapezoidRule (const DoubleReal a, const DoubleReal b, const DoubleReal fa, const DoubleReal fb) throw ()
			{ 
				return ((fb+fa)*0.5*(b-a)); 
			};

			/** @brief Clusters the seeds stored by push2TmpBox_.
 				* @param candidates A isotope wavelet transformed spectrum. 
 				* @param ref The corresponding original spectrum (w.r.t. @p candidates). 
 				* @param scan_index The index of the scan under consideration (w.r.t. the original map). 
 				* @param max_charge The maximal charge state we will consider. */
			void clusterSeeds (const std::vector<MSSpectrum<PeakType> >& candidates, 
				const MSSpectrum<PeakType>& ref, const UInt scan_index, const UInt max_charge) throw ();


			virtual DoubleReal correctMZ (const MSSpectrum<PeakType>& ref, const DoubleReal c_mz, const DoubleReal c) throw ();

			bool create_Mascot_PMF_File_;

			//internally used data structures for the sweep line algorithm	
			std::multimap<DoubleReal, Box_> open_boxes_, closed_boxes_;	//DoubleReal = average m/z position
			std::vector<std::multimap<DoubleReal, Box_> >* tmp_boxes_; //for each charge we need a separate container
		
			gsl_interp_accel* acc_;
			gsl_spline* spline_; 
	};







	template <typename PeakType>	
	bool comparator (const PeakType& a, const PeakType& b)
	{
		return (a.getIntensity() > b.getIntensity());
	}		


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform () throw()
	{
		create_Mascot_PMF_File_	= false;
		acc_ = gsl_interp_accel_alloc ();
		spline_ = gsl_spline_alloc (gsl_interp_cspline, DEFAULT_NUM_OF_INTERPOLATION_POINTS); 
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box_> > (1);
	}

	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::IsotopeWaveletTransform (const UInt max_charge, const bool create_Mascot_PMF_File) throw()
	{	
		create_Mascot_PMF_File_ = create_Mascot_PMF_File;
		acc_ = gsl_interp_accel_alloc ();
		spline_ = gsl_spline_alloc (gsl_interp_cspline, DEFAULT_NUM_OF_INTERPOLATION_POINTS); 
		tmp_boxes_ = new std::vector<std::multimap<DoubleReal, Box_> > (max_charge);
	}


	template <typename PeakType>
	IsotopeWaveletTransform<PeakType>::~IsotopeWaveletTransform () throw()
	{
		gsl_interp_accel_free (acc_);
		gsl_spline_free (spline_);
		delete (tmp_boxes_);
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::getTransforms (const MSSpectrum<PeakType>& scan, 
		std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge, const Int mode) throw ()
	{	
		UInt scan_size = scan.size(), wavelet_length=0, peak_cutoff=0;
		DoubleReal av_MZ_spacing = getAvMZSpacing_(scan);
		
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
		//causing the upcoming code to crash since suddenly some m/z positions will occure twice. The interval of multiple 
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
				
				IsotopeWavelet::getAveragine (scan[i].getMZ()*c_charge, &peak_cutoff);
				
				wavelet_length = peak_cutoff/av_MZ_spacing;
				if (wavelet_length > scan_size)
				{		
					return;
				};
				std::vector<DoubleReal> psi (wavelet_length, 0); //The wavelet

				//Sampling the wavelet 
				sampleTheWavelet_ (scan, wavelet_length, i, cum_spacing, (UInt) c_charge, av_MZ_spacing, psi, peak_cutoff, mode);
					
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
				transforms[c][i].setIntensity(sums);
				transforms[c][i].setMZ(max_position_scan);	

				if (repair)
				{		
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
					DoubleReal x1 = c_peak2.getMZ(); 
					DoubleReal y1 = c_peak2.getIntensity(); 				
					DoubleReal x2 = c_peak1.getMZ();
					DoubleReal y2 = c_peak1.getIntensity();
					DoubleReal dx = (x2-x1)/(DoubleReal)(noi2interpol);
					for (UInt ii=1; ii<noi2interpol; ++ii) //ii=1, not 0, since we want to keep the first of the multiples
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
				std::stringstream name; name << "trans_" << c+1 << ".dat\0"; 
				std::ofstream ofile (name.str().c_str());
				for (unsigned int i=0; i<transforms[c].size(); ++i)
					ofile << transforms[c][i].getMZ() << "\t" <<  transforms[c][i].getIntensity() << std::endl;
				ofile.close();
			};
		#endif

		return;
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates,
		const MSSpectrum<PeakType>& ref, const UInt scan_index, const DoubleReal ampl_cutoff) throw ()
	{
		UInt peak_cutoff=0; 	
		UInt cands_size=candidates.size();
		UInt signal_size=candidates[0].size(), i_iter; 
		typename MSSpectrum<PeakType>::iterator iter, iter2;
		typename MSSpectrum<PeakType>::const_iterator iter_start, iter_end, iter_p;
		DoubleReal seed_mz, c_av_intens=0, c_score=0, c_sd_intens=0, threshold=0, help_mz;
	 	UInt help_dist, MZ_start, MZ_end;

		//For all charges do ...
		for (UInt c=0; c<cands_size; ++c)		
		{
			MSSpectrum<PeakType> c_sorted_candidate (candidates[c]); 
			std::vector<DoubleReal> processed (signal_size, 0);

			//Sort the transform in descending order according to the intensities present in the transform 	
			sort (c_sorted_candidate.begin(), c_sorted_candidate.end(), comparator<PeakType>); 
			
			if (ampl_cutoff < 0)
			{
				threshold=EPSILON_ION_COUNTS;
			}
			else
			{				
				c_av_intens = getAvIntens_ (c_sorted_candidate);
				c_sd_intens = getSdIntens_ (c_sorted_candidate, c_av_intens);
				threshold=ampl_cutoff*c_sd_intens + c_av_intens;
			};		
				
			//Eliminate uninteresting regions
			for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter)
			{
				if (iter->getIntensity() < c_av_intens)
				{	
					break;
				};
			};

			i_iter=0;
			for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter, ++i_iter)
			{					
				seed_mz=iter->getMZ();
			 
				if (processed[distance(candidates[c].begin(), candidates[c].MZBegin(seed_mz))] > iter->getIntensity())
				{ 	
					continue;
				}

				IsotopeWavelet::getAveragine(seed_mz*(c+1), &peak_cutoff);
				//Mark the region as processed
				//Do not move this further down, since we have to mark this as processed in any case, 
				//even when score <=0; otherwise we would look around the maximum's position unless 
				//any significant point is found
				iter_start = candidates[c].MZBegin(seed_mz-QUARTER_NEUTRON_MASS/(c+1.));
				iter_end = candidates[c].MZEnd(seed_mz+(peak_cutoff-1)-QUARTER_NEUTRON_MASS/(c+1.));
				
			
				for (iter_p=iter_start; iter_p!=iter_end; ++iter_p)
				{
					help_dist = distance(candidates[c].begin(), iter_p);
					processed[help_dist] = processed[help_dist] >= iter->getIntensity() ?  processed[help_dist] : iter->getIntensity();
				};			


				c_score = scoreThis (candidates[c], peak_cutoff, seed_mz, c, iter->getIntensity(), threshold);

				if (c_score <= 0)
				{
					continue;
				};

				MZ_start = distance (candidates[c].begin(), iter_start);
				MZ_end = distance (candidates[c].begin(), iter_end);

				//Push the seed into its corresponding box (or create a new one, if necessary)
				//Do ***NOT*** move this further down!
				push2TmpBox_ (seed_mz, scan_index, c, c_score, iter->getIntensity(), candidates[0].getRT(), MZ_start, MZ_end);
				
				for (Int h=-2; h<=2; ++h)
				{
					if (h!=0)
					{
						help_mz = seed_mz + h*NEUTRON_MASS/(c+1.);
						iter2 = c_sorted_candidate.MZBegin (help_mz);
						if (fabs(iter2->getMZ()-help_mz) > fabs(h)*NEUTRON_MASS/(2*(c+1.)))
						{
							typename MSSpectrum<PeakType>::const_iterator iter3 (candidates[c].MZBegin(help_mz));
							if (iter3 != candidates[c].end())
							{
								push2TmpBox_ (iter3->getMZ(), scan_index, c, 0, iter3->getIntensity(), candidates[0].getRT(), MZ_start, MZ_end);
							};
						};
					};
				};
			};	
		};

		clusterSeeds(candidates, ref, scan_index, candidates.size());
	}
		

	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::sampleTheWavelet_ (const MSSpectrum<PeakType>& scan, const UInt wavelet_length, 
		const UInt mz_index, const DoubleReal offset, const UInt charge, const DoubleReal av_MZ_spacing, 
		std::vector<DoubleReal>& psi, const DoubleReal peak_cutoff, const Int mode) throw ()
	{
		UInt scan_size = scan.size();
		DoubleReal c_pos, lambda;

		psi.resize (wavelet_length); //just to be sure; if psi is already scan_size large, this will is a simple test
		
		c_pos = scan[mz_index].getMZ();				
		lambda = IsotopeWavelet::getLambdaQ(c_pos*charge-mode*charge*PROTON_MASS);

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
		DoubleReal tz1;
		for (UInt j=0; j<wavelet_length; ++j)
		{
			tz1=cum_spacing*charge+1;
			psi[j] = (cum_spacing > 0) ? IsotopeWavelet::getValueByLambda (lambda, tz1) : 0;
			cum_spacing += c_spacings[j];
		}

		DoubleReal mean=0;
		for (UInt j=0; j<wavelet_length-1; ++j)
		{
			mean += chordTrapezoidRule (scan[(mz_index+j)%scan_size].getMZ(), scan[(mz_index+j+1)%scan_size].getMZ(), psi[j], psi[j+1]);
		}

		for (UInt j=0; j<wavelet_length; ++j)
		{
			psi[j] -= mean/peak_cutoff;
		}

		#ifdef DEBUG_FEATUREFINDER
			if (trunc(c_mzs[0]) == 1000 || trunc(c_mzs[0]) == 1700 || trunc(c_mzs[0]) == 2000 || trunc(c_mzs[0]) == 3000)
			{
				std::stringstream stream; stream << "wavelet_" << c_mzs[0] << "_" << charge+1 << ".dat\0"; 
				std::ofstream ofile (stream.str().c_str());
				for (unsigned int i=0; i<wavelet_length; ++i)
				{
					ofile << scan[mz_index+i].getMZ() << "\t" << psi[i] << std::endl;
				}
				ofile.close();
			};
		#endif
	}

	
	template<typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::scoreThis (const MSSpectrum<PeakType>& candidate, 
		const UInt peak_cutoff, const DoubleReal seed_mz, const UInt c, const DoubleReal intens, const DoubleReal ampl_cutoff) throw ()
	{	
		DoubleReal c_score=0, c_check_point=-1, c_val;
		std::pair<int, int> c_between_left, c_between_right; 				
		typename MSSpectrum<PeakType>::const_iterator c_left_iter, c_right_iter, c_iter;

		DoubleReal leftBound, rightBound;
		//p_h_ind indicates if we are looking for a whole or a peak
		Int p_h_ind=1, end=4*(peak_cutoff-1) -1; //4 times and not 2 times, since we move by 0.5 m/z entities

		std::vector<DoubleReal> xvec, yvec, weights;
		std::vector<DoubleReal> xs (DEFAULT_NUM_OF_INTERPOLATION_POINTS), ys(DEFAULT_NUM_OF_INTERPOLATION_POINTS);
		UInt i=0;
		
		for (Int v=1; v<end; ++v, ++p_h_ind)
		{		
			c_check_point = seed_mz-((peak_cutoff-1)*NEUTRON_MASS-v*0.5*NEUTRON_MASS)/((DoubleReal)c+1);
			
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

		#ifdef DEBUG_FEATUREFINDER	
			std::ofstream ofile_score ("scores.dat", ios::app);
			std::ofstream ofile_check_score ("check_scores.dat", ios::app);
			ofile_score << c_check_point << "\t" << c_score << std::endl;
			ofile_check_score << c_check_point << "\t" << c_check_score << std::endl;
			ofile_score.close();
			ofile_check_score.close();
		#endif

		if (c_score <= ampl_cutoff+intens)
		{
			return(0);
		};

		return (c_score);
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
	DoubleReal IsotopeWaveletTransform<PeakType>::getAvIntens_ (const MSSpectrum<PeakType>& scan) throw ()
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
	DoubleReal IsotopeWaveletTransform<PeakType>::getSdIntens_ (const MSSpectrum<PeakType>& scan, const DoubleReal mean) throw ()
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
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) throw ()
	{	
		if (intens <= 0)
		{		
			#ifdef DEBUG_FEATUREFINDER
				std::cout << "Warning: detected candidate with zero ion couts at m/z: " << mz << std::endl;
			#endif
			return;
		};

		typename std::map<DoubleReal, Box_>::iterator upper_iter = open_boxes_.upper_bound(mz);
		typename std::map<DoubleReal, Box_>::iterator lower_iter; 
		
		lower_iter = open_boxes_.lower_bound(mz);
		if (lower_iter != open_boxes_.end())
		{
			//Ugly, but necessary due to the implementation of STL lower_bound
			if (mz != lower_iter->first)
			{
				lower_iter = --(open_boxes_.lower_bound(mz));
			};
		};
		
		typename std::map<DoubleReal, Box_>::iterator insert_iter;
		bool create_new_box=true;
		if (lower_iter == open_boxes_.end()) //I.e. there is no open Box for that mz position
		{
			//There is another special case to be considered here:
			//Assume that the current box contains only a single element that is (slightly) smaller than the new mz value, 
			//then the lower bound for the new mz value is box.end and this would usually force a new entry
			if (!open_boxes_.empty())
			{
				if (fabs((--lower_iter)->first - mz) < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) //matching box
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
			if (upper_iter == open_boxes_.end() && fabs(lower_iter->first - mz) < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) //Found matching Box
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
			dist_lower = (dist_lower < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) ? dist_upper : INT_MAX;

			if (dist_lower>=0.5*NEUTRON_MASS/(/*charge+*/1.0) && dist_upper>=0.5*NEUTRON_MASS/(/*charge+*/1.0)) // they are both too far away
			{
				create_new_box=true;
			}
			else
			{
				insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;	
				create_new_box=false;
			};
		}; 
	
		BoxElement_ element; 
		element.c = charge; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end; 


		if (create_new_box == false)
		{						
			double max=intens;
			typename Box_::iterator box_iter;
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
			
			std::pair<UInt, BoxElement_> help2 (scan, element);
			insert_iter->second.insert (help2);	

			//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
			Box_ replacement (insert_iter->second);	

			//We cannot devide both m/z by 2, since we already inserted some m/z's whose weight would be lowered.
			//Also note that we alread inserted the new entry, leading to size-1.
			DoubleReal c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;	
			c_mz /= ((DoubleReal) insert_iter->second.size());		

			//Now let's remove the old and insert the new one
			open_boxes_.erase (insert_iter);	
			std::pair<DoubleReal, std::multimap<UInt, BoxElement_> > help3 (c_mz, replacement);
			open_boxes_.insert (help3);				
		}
		else
		{			
			element.max_intens=intens;			
			std::pair<UInt, BoxElement_> help2 (scan, element);
			std::multimap<UInt, BoxElement_> help3;
			help3.insert (help2);
			std::pair<DoubleReal, std::multimap<UInt, BoxElement_> > help4 (mz, help3);
			open_boxes_.insert (help4);
		};
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::push2TmpBox_ (const DoubleReal mz, const UInt scan, UInt charge, 
		const DoubleReal score, const DoubleReal intens, const DoubleReal rt, const UInt MZ_begin, const UInt MZ_end) throw ()
	{	
		std::multimap<DoubleReal, Box_>& tmp_box (tmp_boxes_->at(charge));
		typename std::map<DoubleReal, Box_>::iterator upper_iter = tmp_box.upper_bound(mz);
		typename std::map<DoubleReal, Box_>::iterator lower_iter; 
	
		lower_iter = tmp_box.lower_bound(mz);
		if (lower_iter != tmp_box.end())
		{
			//Ugly, but necessary due to the implementation of STL lower_bound
			if (mz != lower_iter->first)
			{
				lower_iter = --(tmp_box.lower_bound(mz));
			};
		};
		
		typename std::map<DoubleReal, Box_>::iterator insert_iter;
		bool create_new_box=true;
		if (lower_iter == tmp_box.end()) //I.e. there is no tmp Box for that mz position
		{
			//There is another special case to be considered here:
			//Assume that the current box contains only a single element that is (slightly) smaller than the new mz value, 
			//then the lower bound for the new mz value is box.end and this would usually force a new entry
			if (!tmp_box.empty())
			{
				if (fabs((--lower_iter)->first - mz) < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) //matching box
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
			if (upper_iter == tmp_box.end() && fabs(lower_iter->first - mz) < 0.5*NEUTRON_MASS/(/*charge+*/1.0)) //Found matching Box
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
			dist_lower = (dist_lower < 0.5*NEUTRON_MASS/(charge+1.0)) ? dist_lower : INT_MAX;
			dist_upper = (dist_upper < 0.5*NEUTRON_MASS/(charge+1.0)) ? dist_upper : INT_MAX;

			if (dist_lower>=0.5*NEUTRON_MASS/(charge+1.0) && dist_upper>=0.5*NEUTRON_MASS/(charge+1.0)) // they are both too far away
			{
				create_new_box=true;
			}
			else
			{
				insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;	
				create_new_box=false;
			};
		}; 

		BoxElement_ element; 
		element.c = charge; element.mz = mz; element.score = score; element.RT = rt; element.intens=intens;
		element.RT_index = scan; element.MZ_begin = MZ_begin; element.MZ_end = MZ_end; 
		
		if (create_new_box == false)
		{			
			double max=intens;
			typename Box_::iterator box_iter;
			for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
				if (box_iter->second.max_intens > max)
					max=box_iter->second.max_intens;
			if (max==intens)
			for (box_iter=insert_iter->second.begin(); box_iter != insert_iter->second.end(); ++box_iter)
				box_iter->second.max_intens=intens;
			element.max_intens=max;			
			
			std::pair<UInt, BoxElement_> help2 (scan, element);
			insert_iter->second.insert (help2);	

			//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
			Box_ replacement (insert_iter->second);	

			//We cannot devide both m/z by 2, since we already inserted some m/z's whose weight would be lowered.
			//Also note that we alread inserted the new entry, leading to size-1.
			DoubleReal c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;	
			c_mz /= ((DoubleReal) insert_iter->second.size());		

			//Now let's remove the old and insert the new one
			tmp_box.erase (insert_iter);	
			std::pair<DoubleReal, std::multimap<UInt, BoxElement_> > help3 (c_mz, replacement);	
			tmp_box.insert (help3);				
		}
		else
		{			
			element.max_intens=intens;			
			std::pair<UInt, BoxElement_> help2 (scan, element);
			std::multimap<UInt, BoxElement_> help3;
			help3.insert (help2);
			std::pair<DoubleReal, std::multimap<UInt, BoxElement_> > help4 (mz, help3);
			tmp_box.insert (help4);
		};
	}



	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::updateBoxStates (const UInt scan_index, const UInt RT_interleave, 
		const UInt RT_votes_cutoff) throw ()
	{
		typename std::map<DoubleReal, Box_>::iterator iter, iter2;
		for (iter=open_boxes_.begin(); iter!=open_boxes_.end(); )
		{
			//For each Box we need to figure out, if and when the last RT value has been inserted
			//If the Box his unchanged since RT_interleave_ scans, we will close the Box.
			UInt lastScan = (--(iter->second.end()))->first;
			if (scan_index - lastScan > RT_interleave) //I.e. close the box!
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
	void IsotopeWaveletTransform<PeakType>::clusterSeeds (const std::vector<MSSpectrum<PeakType> >& candidates, 
		const MSSpectrum<PeakType>& ref,  const UInt scan_index, const UInt max_charge) throw ()
	{
		typename std::map<DoubleReal, Box_>::iterator iter;
		typename Box_::iterator box_iter;
		std::vector<BoxElement_> final_box;
	 	DoubleReal c_mz, c_RT, av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, count=0;
		UInt num_o_feature, l_mz, r_mz;

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (unsigned int c=0; c<max_charge; ++c)
		{
			std::vector<BoxElement_> final_box;
			for (iter=tmp_boxes_->at(c).begin(); iter!=tmp_boxes_->at(c).end(); ++iter)
			{	
				Box_& c_box = iter->second;
				std::vector<DoubleReal> charge_votes (max_charge, 0), charge_binary_votes (max_charge, 0);
			
				av_max_intens=0, av_score=0, av_mz=0, av_RT=0, av_intens=0, count=0, l_mz=INT_MAX, r_mz=0;
				//Now, let's get the RT boundaries for the box
				for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
				{
					c_mz = box_iter->second.mz;
					c_RT = box_iter->second.RT;
					av_score += box_iter->second.score;
					av_intens += box_iter->second.intens;
					av_mz += c_mz*box_iter->second.intens;

					if (l_mz > box_iter->second.MZ_begin) l_mz=box_iter->second.MZ_begin;
					if (r_mz < box_iter->second.MZ_end) r_mz=box_iter->second.MZ_end;

					++count;
				};

				av_max_intens = c_box.begin()->second.max_intens;
				av_intens /= count; 
				//in contrast to the key entry of tmp_box_, this mz average is weighted by intensity
				av_mz /= count*av_intens;
				av_score /= count;
				av_RT = c_box.begin()->second.RT;

				BoxElement_ c_box_element;
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

			//Computing the derivates
			std::vector<DoubleReal> bwd_diffs(num_o_feature, 0), fwd_diffs(num_o_feature, 0); //, c_diffs (num_o_feature);
			/*c_diffs[0]=0; if (num_o_feature >= 1) c_diffs[num_o_feature-1]=0; 
			for (UInt i=1; i<num_o_feature-1; ++i)
			{
				c_diffs[i] = (final_box[i+1].max_intens-final_box[i-1].max_intens)/(final_box[i+1].mz-final_box[i+1].mz);
			};*/

			bwd_diffs[0]=0;
			for (UInt i=1; i<num_o_feature; ++i)
			{
				bwd_diffs[i] = (final_box[i].max_intens-final_box[i-1].max_intens)/(final_box[i].mz-final_box[i-1].mz);
			};		

			if (num_o_feature >= 1) fwd_diffs[num_o_feature-1]=0;
			for (UInt i=0; i<num_o_feature-1; ++i)
			{					
				fwd_diffs[i] = (final_box[i+1].max_intens-final_box[i].max_intens)/(final_box[i+1].mz-final_box[i].mz);
			};

			#ifdef DEBUG_FEATUREFINDER
			std::ofstream ofile_bwd ("bwd.dat"), ofile_fwd ("fwd.dat");
			for (unsigned int i=0; i<num_o_feature; ++i)
			{
				ofile_fwd << final_box[i].mz << "\t" << fwd_diffs[i] << std::endl;
				ofile_bwd << final_box[i].mz << "\t" << bwd_diffs[i] << std::endl;
			};
			ofile_bwd.close(); ofile_fwd.close();
			#endif

			for (UInt i=0; i<num_o_feature; ++i)
			{	
				while (i<num_o_feature-1)
				{
					if(final_box[i].score>0) //this has been an helping point
						break;
					++i;
				};

				if (bwd_diffs[i]>0 && fwd_diffs[i]<0) //at the moment we will only use the forward and the backward differences
				{
					checkPosition (candidates[c], ref, final_box[i].mz, final_box[i].c, scan_index);	
					continue;
				};
			};
			tmp_boxes_->at(c).clear();
		};
	}


	template <typename PeakType>
	FeatureMap<Feature> IsotopeWaveletTransform<PeakType>::mapSeeds2Features 
		(const MSExperiment<PeakType>& map, const UInt max_charge, const UInt RT_votes_cutoff) throw ()
	{
		FeatureMap<Feature> feature_map;
		typename std::map<DoubleReal, Box_>::iterator iter;
		typename Box_::iterator box_iter;
		UInt best_charge_index; DoubleReal best_charge_score, c_mz, c_RT; UInt c_charge; 	
		DoubleReal av_intens=0, av_score=0, av_mz=0, av_RT=0, av_max_intens=0; UInt peak_cutoff;
		ConvexHull2D c_conv_hull;
		std::ofstream* ofile = NULL;

		if (create_Mascot_PMF_File_)
		{
			ofile = new std::ofstream ("mascot.query");
		};

		typename std::pair<DoubleReal, DoubleReal> c_extend;
		for (iter=closed_boxes_.begin(); iter!=closed_boxes_.end(); ++iter)
		{		
			Box_& c_box = iter->second;
			std::vector<DoubleReal> charge_votes (max_charge, 0), charge_binary_votes (max_charge, 0);
		
			//Let's first determine the charge
			//Therefor, we can use two types of votes: qulitative ones (charge_binary_votes) or quantitaive ones (charge_votes)
			//Collting the votes ...
			for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
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
        
        IsotopeWavelet::getAveragine (c_mz*c_charge, &peak_cutoff);
				point_set.push_back (DPosition<2> (c_RT, c_mz - QUARTER_NEUTRON_MASS/(DoubleReal)c_charge)); 
				point_set.push_back (DPosition<2> (c_RT, c_mz + ((peak_cutoff+0.5)*NEUTRON_MASS)/(DoubleReal)c_charge)); 
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

			//Is this weightening here really clever?
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
			if (create_Mascot_PMF_File_)
			{
				if (map.size()==1) //PMF
				{
					*ofile << ::std::setprecision(4) << std::fixed << av_mz << "\t" << av_max_intens << std::endl;
				}
				else
				{
					*ofile << ::std::setprecision(4) << std::fixed << av_RT << "\t" << av_mz << "\t" << av_max_intens << "\t" << c_charge <<  "\t" << av_score 
						<< std::endl;
				}
			};	
		};		

		if (create_Mascot_PMF_File_)
		{
			ofile->close();
			delete (ofile);
		};
		return (feature_map);
	}


	template <typename PeakType>
	void IsotopeWaveletTransform<PeakType>::checkPosition (const MSSpectrum<PeakType>& candidate,
		const MSSpectrum<PeakType>& ref, const DoubleReal seed_mz, const UInt c, const UInt scan_index) throw ()
	{
		typename MSSpectrum<PeakType>::const_iterator right_cutoff, seed, iter, pre_iter, post_iter;
		UInt peak_cutoff;
		IsotopeWavelet::getAveragine(seed_mz*(c+1), &peak_cutoff);

		right_cutoff = candidate.MZBegin(seed_mz+(peak_cutoff-1)*NEUTRON_MASS/(c+1.));
		pre_iter = candidate.MZBegin(seed_mz-NEUTRON_MASS/(c+1.));
		seed = candidate.MZBegin(seed_mz);
		post_iter = candidate.MZBegin(seed_mz+NEUTRON_MASS/(c+1.));
		iter=seed; 
		//we can ignore those cases
		if (iter==candidate.begin() || iter==candidate.end() || pre_iter==candidate.begin() || post_iter == candidate.end()) 
		{
			return;
		};

		DoubleReal normal = pre_iter->getIntensity() / post_iter->getIntensity();
		DoubleReal pre = pre_iter->getIntensity() / iter->getIntensity();

		if (fabs(1-pre) <= fabs(1-normal)) //okay, let's move this peak by 1 Da to the left ...
		{
			//... but first check if the signal might be caused by an overlapping effect
			if ((candidate.MZBegin(pre_iter->getMZ()-NEUTRON_MASS/(c+1.)))->getIntensity() < pre_iter->getIntensity())
			{	
				//okay, these hard coded values should be checked again, but they definitely cover the *first* critical range
				if (seed_mz > 1500 && seed_mz < 2500)
				{
					iter = candidate.MZBegin(iter->getMZ()-NEUTRON_MASS/(c+1.));
				};
			};
		};

		if (iter->getIntensity() < 1)
		{
			return;
		};

		DoubleReal c_score = scoreThis (candidate, peak_cutoff, iter->getMZ(), c, iter->getIntensity(), 0);
		
		//Correct the position
		DoubleReal real_MZ = correctMZ (ref, iter->getMZ(), c);
		typename MSSpectrum<PeakType>::const_iterator real_l_MZ_iter = ref.MZBegin(real_MZ-QUARTER_NEUTRON_MASS/(c+1.));		
		typename MSSpectrum<PeakType>::const_iterator real_r_MZ_iter = ref.MZBegin(real_MZ+(peak_cutoff-1)*NEUTRON_MASS/(c+1.));

		UInt real_MZ_begin = distance (ref.begin(), real_l_MZ_iter);
		UInt real_MZ_end = distance (ref.begin(), real_r_MZ_iter);

		push2Box_ (real_MZ, scan_index, c, c_score, iter->getIntensity(), ref.getRT(), real_MZ_begin, real_MZ_end);
	}

	
	template <typename PeakType>
	DoubleReal IsotopeWaveletTransform<PeakType>::correctMZ (const MSSpectrum<PeakType>& ref, const DoubleReal c_mz, const DoubleReal /*c*/) throw ()
	{
 		
		return (c_mz);

		typename MSSpectrum<PeakType>::const_iterator iter = ref.MZBegin(c_mz);
		typename MSSpectrum<PeakType>::const_iterator liter=iter, riter=iter;

		while (1)
		{
			--liter;
			if (liter->getIntensity() > (liter+1)->getIntensity())
			{
				break;
			};
		};

		while (1)
		{
			++riter;
			if (riter->getIntensity() > (riter-1)->getIntensity())
			{
				break;
			};
		};

		
		if (liter == ref.end() || riter == ref.end())
		{
			return(0);
		};

		MSSpectrum<PeakType> tofit; tofit.assign (liter, ++riter);
		DoubleReal norm=0;
		for (UInt i=0; i<tofit.size(); ++i)
			norm += tofit[i].getIntensity();			

		for (UInt i=0; i<tofit.size(); ++i)
			tofit[i].setIntensity(tofit[i].getIntensity()/norm);

		//DoubleReal width = getSdIntens_(tofit, getAvIntens_(tofit));
		DoubleReal ppm = (272.706148207636 - 0.009125417520757174*c_mz)/(-0.0453027411590229 + 
      0.013066460791702246*c_mz);
		//DoubleReal ppm = (12337.641655705787 - 3.003986551355764*c_mz)/(-811.3855279967621 + 1.663571470796717*c_mz); 
		DoubleReal real_mass = c_mz* (pow(10,6)-0.5*ppm)/(pow(10,6)+0.5*ppm); 
//		std::cout << c_mz << "\t" << width*sqrt(2) << "\t" << real_mass << std::endl;

		return (real_mass);
	}
} //namespace

#endif 
