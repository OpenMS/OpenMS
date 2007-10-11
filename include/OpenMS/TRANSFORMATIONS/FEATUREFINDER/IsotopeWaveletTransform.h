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
#include <gsl/gsl_spline.h>
#include <vector>
#include <map>
#include <hash_map.h>

#ifndef DEFAULT_HASH_PRECISION //setting the default value
#define DEFAULT_HASH_PRECISION 10000 //i.e. float precision (w.r.t. the m/z coordindate) is 4
#endif


namespace OpenMS
{

class IsotopeWaveletTransform
{
	public:

			struct BoxElement
			{			
				double mz;
				unsigned int c; //Note, this is not the charge (it is charge-1!!!)
				double score;
			};

			typedef std::map<unsigned int, BoxElement> Box; //Key: RT, value: BoxElement
			typedef DRawDataPoint<2> RawDataPoint2D; 

			/** @brief Default Constructor. */
			IsotopeWaveletTransform () throw();

			/** @brief Destructor. */
			virtual ~IsotopeWaveletTransform () throw();
		
			/** @brief Computes the discrete-time continuous wavelet transform simultaneously for several charges.
 				* 
 				* The function computes the isotope wavelet transformed versions of @p scan. 
 				* The function computes the transform for several charge states (up to charge @p max_charge) at the same time.
 				* Hence, the user has to ensure that the size of @p transforms equals to @p max_charge and that each spectrum in 
 				* @p transforms has the same length as @p scan.
 				*
 				* @param scan The MS scan you wish to transform.
 				* @param transforms A vector (with indices running from 0 to max_charge-1) of MS spectra (each of the size of scan).
 				* The code will NOT check the allocated memory (the sizes) for transforms and its entries.  
 				* @param max_charge The maximal charge state that is considered.
 				* @return The transformed spectra. Entry i in this vector corresponds to the "i+1"-charge-state-transform 
 				* of @p scan. The length of each MSSpectrum contained in this vector equals either the length of scan. 	
 				* If @param scan is smaller than the internally computed wavelet length, no transform can be computed and 
				* the original scan will be returend for each charge state. */
			static void getTransforms (const MSSpectrum<RawDataPoint1D>& scan, 
				std::vector<MSSpectrum<RawDataPoint1D> > &transforms, const unsigned int max_charge) throw ();


			/** @brief Given an isotope wavelet transformed spectrum @p candidates, this function assigns to every significant
 				* pattern its corresponding charge and a score indicating the reliability of the prediction. The result of this
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
			
			
			virtual double scoreThis (const MSSpectrum<RawDataPoint1D>& candidate, unsigned int start_index, 
				const double seed_mz, const unsigned int end_index, const unsigned int c, const double intens, const double ampl_cutoff=0) throw ();
			

			/** @brief Returns a floating point precision, which is internally used for a hash map. 
				*	The default value is 10000, which corresponds to a floating point precision of 4. Usually, this value can be left unchanged. 		
				*	If you are considering some very high resulted spectra (e.g. FT-ICR) you might increase this value to 5 or 6. The parameter
				*	is used to mark regions of interest for which an isotopic pattern has already been found. Hence, these regions will not 
				*	be considered during the rest of the feature finding algorithm. */
			inline double getHashPrecision () const throw ()
				{ return (hash_precision_); }

			/** @brief See getHashPrecision. */
			void setHashPrecision (const double hash_precision) throw ()
				{ hash_precision_ = hash_precision; }

			/** @brief Computes the average MZ spacing of @p scan. */
			static double getAvMZSpacing (const MSSpectrum<RawDataPoint1D>& scan, int start_index=0, int end_index=-1) throw ();

	protected:			

			/** @brief Samples the wavelet at discrete time points, s.t. they match automatically to the m/z positions provided
 				* in @p scan and returns the discrete values of psi in @psi. Usually (unless you would like to write your own 
 				* @see FeatureFinder), you do not need to call this function.
 				*	@param mz_index The start index of @p scan for which the the wavelet should be adapted.  
 				* @param offset Is the offset the wavelet function needs to be aligned with a signal point.
 				* @param z The charge the wavelet function should adapt (corresponds to z in the paper).
 				* @param av_MZ_spacing The average spacing in the m/z domain.
 				* @param psi The sampled values.
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */ 
			static void sampleTheWavelet (const MSSpectrum<RawDataPoint1D>& scan, const unsigned int mz_index, 
				const double offset, const unsigned int z, const double av_MZ_spacing, std::vector<double>& psi, 
				const unsigned int mode=+1) throw ();		

			/** @brief Trapezoid rule */
			static double chordTrapezoidRule (const double a, const double b, const double fa, const double fb) throw ()
				{ return ((fb+fa)*0.5*(b-a)); };				

			/** @brief Computes the average intensity of @p scan. */
			inline double getAvIntens (const MSSpectrum<RawDataPoint1D>& scan) throw (); 		
			
			/** @brief Computes the standard deviation of the intensity of @p scan. */
			inline double getSdIntens (const MSSpectrum<RawDataPoint1D>& scan, const double mean) throw ();

			/** @brief A wrapper function to the GSL interpolation routine. */
			double getCubicInterpolatedValue (const std::vector<double>& x, const double xi, const std::vector<double>& y) throw ();

			/** @brief A function to map m/z values to m/z indices. In particular useful if you know already the
 				* approximate position of corresponding entry which can be indicated by @p start. */ 		
			inline std::pair<int, int> getNearBys (const MSSpectrum<RawDataPoint1D>& signal, const double mz, 
				const unsigned int start=0) const throw ();

			virtual void push2Box (const double mz, const unsigned int scan, unsigned int charge, const double score) throw ();


			public:

				void updateBoxStates (const unsigned int c_scanNumber, const unsigned int RT_interleave, 
					const unsigned int RT_votes_cutoff) throw ();

				std::list<RawDataPoint2D> mapSeeds2Features (const unsigned int max_charge, const unsigned int RT_votes_cutoff) throw ();

			protected:

			//A floating point precision for an internally used hash map
			double hash_precision_;	

			//Internally used data structures for the sweep line algorithm	
			std::map<double, Box> openBoxes_, closedBoxes_;	//double = average m/z position
	
};

} //namespace

#endif 
