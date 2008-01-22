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
#define QUARTER_NEUTRON_MASS 0.25059
#define WAVELET_PERIODICITY 6.26845
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
#define SHIFT23 (1<<23)
#define SHIFT23_00 (1.0/(1<<23))
#define LOG_CONST 0.346607f;
#define POW_CONST 0.33971f;
#endif
#endif


#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include<vector>


namespace OpenMS
{
	/** @brief Implements the isotope wavelet function.
	 *
	 * 	The IsotopeWavelet class implements the isotope wavelet as described by R.Hussong, A.Tholey, A.Hildebrandt: 
	 * 	Efficient Analysis of Mass Spectrometry Data Using the Isotope Wavelet. Proceedings of the 3rd international 
	 * 	Symposium in Computational Life Sciences (Complife07). American Institue of Physis (AIP) Proceedings (2007).
	 *
	 * 	The class mainly provides static functions in order to be easily accessible by other classes with seeding or extending
	 * 	purposes. 
	 *
	 *  @ingroup FeatureFinder
	 * 	@todo Tests for negative mode. */
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
				static DoubleReal getValueByMass (const DoubleReal t, const DoubleReal m, const UInt z, const Int mode=+1) throw ();			

				/** @brief Returns the value of the isotope wavelet at position @p t. Usually, you do not need to call this function.
					* Please use sampleTheWavelet instead.			 
					* 
					* Note that this functions returns the pure function value of psi and not the normalized (average=0)
					* value given bei Psi. 
					* @param lambda The mass-parameter lambda.
					* @param tz1 t (the position) times the charge (z) plus 1. */ 
				static DoubleReal getValueByLambda (const DoubleReal lambda, const DoubleReal tz1) 
					throw ();


				/** @brief Returns the max_charge_ parameter. */
				static UInt getMaxCharge () throw ()
				{ 
					return (max_charge_); 
				}			
			
				/** @brief Sets the max_charge_ parameter. */
				static void setMaxCharge (const UInt max_charge) throw ()
				{ 
					max_charge_ = max_charge; 
				}	
		
				/** @brief Returns the table_steps_ parameter. */
				static DoubleReal getTableSteps () throw ()
				{ 
					return (table_steps_); 
				}			
				
				/** @brief Sets the table_steps_ parameter. */ 
				static void setTableSteps (const DoubleReal table_steps) throw ()
				{
					inv_table_steps_ = 1./table_steps;
					table_steps_ = table_steps; 
				}

				/** @brief Should be called once before values are drawn from the isotope wavelet function. 
					*
					* The function precomputes the expensive gamma function. Parameters related to this function are:
					* max_charge_ and peak_cutoff_. If both of these are set correctly @see getValue will never compute
					* the gamma function online. Please note that in a future and more efficient version checks for precomputed
					* values will be removed. 
					*
					* @param max_m The maximal deconvoluted mass that occures in the current data set. */
				static void preComputeExpensiveFunctions (const DoubleReal max_m) throw ();

				/** @brief Returns the mass-parameter lambda (linear fit). 
					* @note The only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
				static DoubleReal getLambdaL (const DoubleReal m) throw ();

				/** @brief Returns the mass-parameter lambda (quadratic fit). 
					* @note The only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
				static DoubleReal getLambdaQ (const DoubleReal m) throw ();					

				/** @brief Initializes the internally used averagine model; automatically called by the FeatureFinder.
 					* @param max_mz The maximal deconvoluted mass that occures in the current data set.	*/ 
				static void computeIsotopeDistributionSize (const DoubleReal max_m) throw ();

				/** @brief Computes the averagine isotopic distribution we would expect the deconvoluted mass. 
 					* @param m The deconvoluted mass m.	
 					* @param size Returns the number of significant peaks within a pattern occuring at mass @p m.
 					* @return The isotopic distribution. */ 
				static const IsotopeDistribution::ContainerType& getAveragine (const DoubleReal m, UInt* size=NULL) throw ();


		protected:

				/** @brief Internal function using register shifts for fast computation of the power function. 
					* @note Please, do not modify this function. */
				static float myPow_ (float a, float b) throw ();			
				
				#ifndef OPENMS_64BIT_ARCHITECTURE	
					/** @brief Internal function using register shifts for fast computation of the power function. 
						*	The function follows http://www.dctsystems.co.uk/Software/power.html, code by 
						* Ian Stephenson, DCT Systems, NCCA ournemouth University.
					 	* @note Please, do not modify this function. */
					static float myPow2_ (float i) throw ();
			
					/** @brief Internal function using register shifts for fast computation of the power function. 
					  * The function follows http://www.dctsystems.co.uk/Software/power.html, code by
					  * Ian Stephenson, DCT Systems, NCCA ournemouth University.
					 	* @note Please, do not modify this function. */
					static float myLog2_ (float i) throw ();

					/** @brief Internal union for fast computation of the power function. */
					union fi_
					{
						Int i;
						float f;
					};
				#endif
			
				/** This parameter determines the maximal charge we will consider.
					* @todo At the moment each from starting from 1 to max_charge_ will be considered for a wavelet transfrom.
					* It might be useful to pass a set of UIntegers to fix the charges. */
				static UInt max_charge_; 				

				/** This parameter determines the sample rate for the precomputation of the gamma function.
					* @note For microTOF or similar well resolved data it might be usuful to decrease this value (by powers of 10). */
				static DoubleReal table_steps_, inv_table_steps_;

				/** Internal table for the precomputed values of the gamma function. */ 
				static std::vector<DoubleReal> gamma_table_;
				
				/** Internal table for the precomputed values of the exponential function. */ 
				static std::vector<DoubleReal> exp_table_;

				/** Internally used averagine model. */
				static IsotopeDistribution averagine;
	};

} //namespace

#endif 
