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


#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <map>
#include <hash_map.h>

#ifndef DEFAULT_HASH_PRECISION //setting the default value
#define DEFAULT_HASH_PRECISION 1000 //i.e. float precision (w.r.t. the m/z coordindate) is 3
#endif


namespace OpenMS
{
/** @brief A class implementing the isotope wavelet transform. 
 * If you just want to find features using the isotope wavelet, take a look at the IsotopeWaveletFF class. Usually, you only
 * have to consider this class if you plan to change the basic implementation of the transform.  */ 
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
			static void getTransforms (const MSSpectrum<RawDataPoint1D>& scan, 
				std::vector<MSSpectrum<RawDataPoint1D> > &transforms, const unsigned int max_charge) throw ();


			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
 				* process is stored internally.
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
			virtual void identifyCharges (const std::vector<MSSpectrum<RawDataPoint1D> >& candidates, const unsigned int scan_index, 
				const double ampl_cutoff=0) throw ();
			
			
			/** @brief Given a candidate for an isotopic pattern, this function computes the corresponding score 
 				* @param candidate A isotope wavelet transformed spectrum.
 				* @param start_index One of the indices enclosing the interesting region.
 				* @param end_index One of the indices enclosing the interesting region.
 				* @param seed_mz The predicted position of the monoisotopic peak.
 				* @param c The charge state minus 1 (e.g. c=2 means charge state 3) for which the score should be determined. 
 				* @param intens The intensity of the transform at @p seed_mz.
 				* @param The threshold. */
			virtual double scoreThis (const MSSpectrum<RawDataPoint1D>& candidate, unsigned int start_index, 
				const double seed_mz, const unsigned int end_index, const unsigned int c, const double intens, const double ampl_cutoff=0) throw ();
			

			/** @brief Returns a floating point precision, which is internally used for a hash map. 
				*	The default value is 1000, which corresponds to a floating point precision of 3. Usually, this value can be left unchanged. 		
				*	If you are considering some very high resulted spectra (e.g. FT-ICR) you might increase this value to 4 or 5. The parameter
				*	is used to mark regions of interest for which an isotopic pattern has already been found. Hence, these regions will not 
				*	be considered during the rest of the feature finding algorithm. */
			inline double getHashPrecision () const throw ()
				{ return (hash_precision_); }

			/** @brief See getHashPrecision. */
			void setHashPrecision (const double hash_precision) throw ()
				{ hash_precision_ = hash_precision; }
			
			/** @brief A function keeping track of currently open and closed sweep line boxes. 
 				* This function is used by the isotope wavelet feature finder and must be called for each processed scan. 
 				* @param c_scan_number The index of the scan currently under consideration w.r.t. its MS map. 
 				* This information is necessary to sweep across the map affter each scan has been evaluated. 
 				* @param See the IsotopeWaveletFF class. 
 				* @param The sweep line cutoff parameter. See also IsotopeWaveletFF class.*/	
			void updateBoxStates (const unsigned int c_scan_number, const unsigned int RT_interleave, 
					const unsigned int RT_votes_cutoff) throw ();

			/** @brief Maps the internally used data structures to the OpenMS framework. 
 				* @param max_charge The maximal charge state under consideration. 
 				* @param The sweep line cutoff parameter. See also IsotopeWaveletFF class.*/
			FeatureMap<Feature> mapSeeds2Features (const unsigned int max_charge, const unsigned int RT_votes_cutoff) throw ();
	

	protected:			

			struct BoxElement
			{			
				double mz;
				unsigned int c; //Note, this is not the charge (it is charge-1!!!)
				double score;
				double intens;
				double RT; //The elution time (not the scan index)
			};

			typedef std::map<unsigned int, BoxElement> Box; //Key: RT index, value: BoxElement
			typedef DRawDataPoint<2> RawDataPoint2D; 

			/** @brief Samples the wavelet at discrete time points, s.t. they match automatically to the m/z positions provided
 				* in @p scan and returns the discrete values of psi in @psi.  
 				*	@param mz_index The start index of @p scan for which the the wavelet should be adapted.  
 				* @param offset Is the offset the wavelet function needs to be aligned with a signal point.
 				* @param z The charge the wavelet function should adapt (corresponds to z in the paper).
 				* @param av_MZ_spacing The average spacing in the m/z domain.
 				* @param psi The sampled values.
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */ 
			static void sampleTheWavelet (const MSSpectrum<RawDataPoint1D>& scan, const unsigned int mz_index, 
				const double offset, const unsigned int z, const double av_MZ_spacing, std::vector<double>& psi, 
				const unsigned int mode=+1) throw ();		

			/** @brief Computes the average intensity (neglecting negative values) of @p scan. */
			inline double getAvIntens (const MSSpectrum<RawDataPoint1D>& scan) throw (); 		
			
			/** @brief Computes the standard deviation (neglecting negative values) of the intensity of @p scan. */
			inline double getSdIntens (const MSSpectrum<RawDataPoint1D>& scan, const double mean) throw ();

			/** @brief A wrapper function to the GSL interpolation routine. */
			double getCubicInterpolatedValue (const std::vector<double>& x, const double xi, const std::vector<double>& y) throw ();

			/** @brief A function to map m/z values to m/z indices. In particular useful if you know already the
 				* approximate position of the corresponding entry which can be indicated by @p start. */ 		
			inline std::pair<int, int> getNearBys (const MSSpectrum<RawDataPoint1D>& signal, const double mz, 
				const unsigned int start=0) const throw ();

			/** @brief Inserts a potential isotopic pattern into an open box or - if no such box exists - creates a new one.
 				* @param mz The position of the pattern.
 				* @param scan The index of the scan, we are currently analyzing (w.r.t. the data map).
 				* This information is necessary for the postprocessing (sweep lining). 
 				* @param The estimated charge state of the pattern. 
 				* @param The pattern's score. 
 				* @param The intensity at the monoisotopic peak. 
 				* @param The retention time of the scan (similar to @p scan, but here: no index, but the real value). */
			virtual void push2Box (const double mz, const unsigned int scan, unsigned int charge, const double score, 
				const double intens, const double rt) throw ();

			/** @brief Computes the average MZ spacing of @p scan in the range @p start_index to @p end_index. */
			static double getAvMZSpacing (const MSSpectrum<RawDataPoint1D>& scan, int start_index=0, int end_index=-1) throw ();

		protected:

			//A floating point precision for an internally used hash map
			double hash_precision_;	

			//Internally used data structures for the sweep line algorithm	
			std::map<double, Box> open_boxes_, closed_boxes_;	//double = average m/z position
	
};

} //namespace

#endif 
