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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_H


#ifndef NEUTRON_MASS
/* According to 
 * D.M. Horn, R.A. Zubarev, F.W. MacLaﬀerty: Automated reduction and interpretation of high resolution 
 * electrospray mass spectra of large molecules. J Am Soc Mass Spectrom 11 (2000) 320–332. */
#define NEUTRON_MASS 1.00235 
// Exact mass of a neutron
#define EXACT_NEUTRON_MASS 1.00866491578
#endif

#ifndef PROTON_MASS
#define PROTON_MASS 1.00727646688
#endif

#ifndef LAMBDA_L_0
//Linear Fit (standard)
#define LAMBDA_L_0 -0.472998839574110749e-1
#define LAMBDA_L_1 0.743579753540513913e-3
#endif

#ifndef LAMBDA_Q_0
//Quadratic Fit (maybe a overkill, since the dependency looks quite linear, at least in moderate mass ranges)
#define LAMBDA_Q_0 -0.137152573151174711
#define LAMBDA_Q_1 0.851289601785403817e-3
#define LAMBDA_Q_2 -0.2834469691e-7
#endif

#ifndef OPENMS_64BIT_ARCHITECTURE	
#ifndef SHIFT_PARAMETERS
//Internal parameters used for fast computation of the power function
//Please do not modify
#define SHIFT_PARAMETERS
#define shift23 (1<<23)
#define shift23_00 (1.0/(1<<23))
#define LogBodge 0.346607f;
#define PowBodge 0.33971f;
#endif
#endif


#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include<vector>


namespace OpenMS
{

/** @brief Implements the isotope wavelet function.
 *
 * 	The IsotopeWavelet class implements the isotope wavelet as described by R.Hussong, A.Tholey, A.Hildebrandt: 
 * 	Efficient Analysis of Mass Spectrometry Data Using the Isotope Wavelet. Proceedings of the 3rd International 
 * 	Symposium in Computational Life Sciences (Complife07). American Institue of Physis (AIP) Proceedings (2007).
 *
 * 	The class mainly provides static functions in order to be easily accessible by other classes with seeding or extending
 * 	purposes. */
class IsotopeWavelet
{
	public:
		
			/** @brief Default Constructor. */
			IsotopeWavelet () throw();

			/** @brief Destructor. */
			~IsotopeWavelet () throw();

			/** @brief Returns the value of the isotope wavelet at position @p t. Usually, you do not need to call this function.
 				* Please use sampleTheWavelet instead.			 
 				* 
 				* Note that this functions returns the pure function value of psi and not the normalized (average=0)
 				* value given bei Psi. 
 				* @param t The position at which the wavelet has to be drawn (within the coordinate system of the wavelet). 
 				* @param m The m/z position within the signal (i.e. the mass not decharged) within the signal.
 				* @param z The charge z we want to detect. 
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */
			static double getValueByMass (const double t, const double m, const unsigned int z, const int mode=+1) throw ();			

			/** @brief Returns the value of the isotope wavelet at position @p t. Usually, you do not need to call this function.
 				* Please use sampleTheWavelet instead.			 
 				* 
 				* Note that this functions returns the pure function value of psi and not the normalized (average=0)
 				* value given bei Psi. 
 				* @param t The position at which the wavelet has to be drawn (within the coordinate system of the wavelet). 
 				* @param lambda The mass-parameter lambda.
 				* @param z The charge z we want to detect. 
 				* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */
			static double getValueByLambda (const double t, const double lambda, const unsigned int z) 
				throw ();

			/** @brief Returns the peak_cutoff_ parameter. */ 
			static unsigned int getPeakCutoff () throw ()
				{ return (peak_cutoff_); }			
			
			/** @brief Sets the peak_cutoff_ parameter. */
 				static unsigned int setPeakCutoff (const unsigned int peak_cutoff) throw ()
				{ return (peak_cutoff_ = peak_cutoff); }			
				

			/** @brief Returns the max_charge_ parameter. */
			static unsigned int getMaxCharge () throw ()
				{ return (max_charge_); }			
		
	  	/** @brief Sets the max_charge_ parameter. */
			static unsigned int setMaxCharge (const unsigned int max_charge) throw ()
				{ return (max_charge_ = max_charge); }	
	
		  /** @brief Returns the gamma_steps_ parameter. */
			static double getGammaSteps () throw ()
				{ return (gamma_steps_); }			
		  
			/** @brief Sets the gamma_steps_ parameter. */ 
			static double setGammaSteps (const double gamma_steps) throw ()
				{ return (gamma_steps_ = gamma_steps); }

			/** @brief Should be called once before values are drawn from the isotope wavelet function. 
 				*
				* The function precomputes the expensive gamma function. Parameters related to this function are:
 				* max_charge_ and peak_cutoff_. If both of these are set correctly @see getValue will never compute
 				* the gamma function online. Please note that in a future and more efficient version checks for precomputed
 				* values will be removed. */
			static void preComputeGammaFunction () throw ();

			/** @brief Returns the mass-parameter lambda (linear fit). 
 				* @note The only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
			static double getLambdaL (const double m) throw ();

			/** @brief Returns the mass-parameter lambda (quadratic fit). 
 				* @note The only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
			static double getLambdaQ (const double m) throw ();					


	protected:

			/** @brief Internal function using register shifts for fast computation of the power function. 
 				* @note Please, do not modify this function. */
			static float myPow (float a, float b) throw ();			
			
			#ifndef OPENMS_64BIT_ARCHITECTURE	
			/** @brief Internal function using register shifts for fast computation of the power function. 
 				* @note Please, do not modify this function. */
			static float myPow2 (float i) throw ();
	
			/** @brief Internal function using register shifts for fast computation of the power function. 
 				* @note Please, do not modify this function. */
			static float myLog2 (float i) throw ();

			/** @brief Internal union for fast computation of the power function */
			union fi
			{
				int i;
				float f;
			};
			#endif
 		
		
			/** This parameter determines the cutoff for the isotope wavelet as far as the m/z dimension is concerned.
 				* peak_cutoff_ correponds to tau in the paper version.
 				* @note The default value is 5. */
			static unsigned int peak_cutoff_; 				

			/** This parameter determines the maximal charge we will consider.
 				* @note The defualt value is 4.  				
 				* @todo At the moment each from starting from 1 to max_charge_ will be considered for a wavelet transfrom.
 				* It might be useful to pass a set of unsigned integers to fix the charges. */
			static unsigned int max_charge_; 				

			/** This parameter determines the sample rate for the precomputation of the gamma function.
 				* @note For microTOF or similar well resolved data it might be usuful to decrease this value (by powers of 10).
 				* The default is 0.0001. */ 
			static double gamma_steps_;

			/** Internal table for the precomputed values of the gamma function. */ 
			static std::vector<double> gamma_table_;
};

} //namespace

#endif 
