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
#include <vector>


namespace OpenMS
{

class IsotopeWaveletTransform
{
	public:
		
			/** @brief Default Constructor. */
			IsotopeWaveletTransform () throw();

			/** @brief Destructor. */
			virtual ~IsotopeWaveletTransform () throw();
		
			/** @brief Computes the discrete-time continuous wavelet transform simultaneously for several charges.
 				* 
 				* The function returns a pointer to the transformed versions of @p scan. Note that after having called this 
 				* function the deletion of this pointer is up to you! The function computes the transform for several charge states
 				* at the same time.
 				* With the help of @p jiggle_bins it is possible to compute an "interpolated" transform, s.t. the maximum of
 				* wavelet is not exactly aligned with some data point. In some (few) cases this could help to find the exact 
 				* m/z position of some feature, since the wavelet even "samples" points not included in the signal. 
 				* Usually, this parameter should be left unchanged, since it involves some expert knowledge to interpret 
 				* the results. Note that the size of each MSSpectrum that the function returns increased of @p jiggle_bins 
 				* is strictly positive.
 				* This function is slightly less efficient than the non static version, since it has to do some preprocessing
 				* on the MS scan. In particular the function computes the average MZ spacing of @scan. Hence, you might get 
 				* (slightly) different results if you set up an instance of this class and provide this average
 				* from yourself. In the latter case you would prefer to compute the average spacing over the whole map instead
 				* on a single scan. 
 				*
 				* @param scan The MS scan you wish to transform.
 				* @param max_charge The maximal charge state that is considered.
 				* @param jiggle_bins The number of additional "virtual sampling points" between each pair of signal points.
 				* @return The transformed spectra. Entry i in this vector corresponds to the "i+1"-charge-state-transform 
 				* of @p scan. The length of each MSSpectrum contained in this vector equals either the length of scan. 	
 				* If @param scan is smaller than the internally computed wavelet length, no transform can be computed and 
				* the original scan will be returend for each charge state. */
			static std::vector<MSSpectrum<RawDataPoint1D> >* getTransforms (const MSSpectrum<RawDataPoint1D>& scan, 
				const unsigned int max_charge, const double jiggle_bins=0) throw ();
			

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
};

} //namespace

#endif 
