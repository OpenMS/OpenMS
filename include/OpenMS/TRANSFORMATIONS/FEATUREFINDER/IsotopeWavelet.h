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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletConstants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

namespace OpenMS
{
	/** @brief Implements the isotope wavelet function.
	 *
	 * 	The IsotopeWavelet class implements the isotope wavelet as described by R. Hussong, A. Tholey, A. Hildebrandt: 
	 * 	Efficient Analysis of Mass Spectrometry Data Using the Isotope Wavelet. Proceedings of the 3rd international 
	 * 	Symposium in Computational Life Sciences (Complife07). American Institute of Physics (AIP) Proceedings (2007).
	 *
	 * 	@note This class features a singleton design pattern. */
	class OPENMS_DLLAPI IsotopeWavelet
	{
		public:
	
				/** The init function; creates an instance of this singleton class. */	
				static IsotopeWavelet* init (const DoubleReal max_m, const UInt max_charge)	;
				
				/** Returns an pointer to the current instance of the class. */	
				static IsotopeWavelet* getInstance ()	
				{
					return (me_);
				}

				/** Deletes the singleton instance. */
				static void destroy ();	


				/** @brief Returns the value of the isotope wavelet at position @p t. Usually, you do not need to call this function.
					* Please use @see sampleTheWavelet instead.			 
					* 
					* Note that this functions returns the pure function value of psi and not the normalized (average=0)
					* value given by Psi. 
					* @param t The position at which the wavelet has to be drawn (within the coordinate system of the wavelet). 
					* @param m The m/z position within the signal (i.e. the mass not decharged) within the signal.
					* @param z The charge @p z we want to detect. 
					* @param mode Indicates whether positive mode (+1) or negative mode (-1) has been used for ionization. */
				static DoubleReal getValueByMass (const DoubleReal t, const DoubleReal m, const UInt z, const Int mode=+1) 
				{			
					return (getValueByLambda (getLambdaL(m*z-z*mode*Constants::IW_PROTON_MASS), t*z+1));
				}
	
				/** @brief Returns the value of the isotope wavelet at position @p t via a fast table lookup. 
  				*	Usually, you do not need to call this function.
					* Please use @see sampleTheWavelet instead.			 
					* 
					* Note that this functions returns the pure function value of psi and not the normalized (average=0)
					* value given by Psi. 
					* @param lambda The mass-parameter lambda.
					* @param tz1 t (the position) times the charge (z) plus 1. */ 
				static DoubleReal getValueByLambda (const DoubleReal lambda, const DoubleReal tz1) ;

				/** @brief Returns the value of the isotope wavelet at position @p t. 
 					* This function is usually significantly slower than the table lookup performed in @see getValueByLambda.
 					* Nevertheless, it might be necessary to call this function due to extrapolating reasons caused by the 
 					* alignment of the wavelet.
 					*
 					* Usually, you do not need to call this function.
					* Please use @see sampleTheWavelet instead.			 
					* 
					* Note that this functions returns the pure function value of psi and not the normalized (average=0)
					* value given by Psi. 
					* @param lambda The mass-parameter lambda.
					* @param tz1 t (the position) times the charge (z) plus 1. */ 
				static DoubleReal getValueByLambdaExtrapol (const DoubleReal lambda, const DoubleReal tz1) ;

				static DoubleReal getValueByLambdaExact (const DoubleReal lambda, const DoubleReal tz1) ; 
				

				/** @brief Returns the largest charge state we will consider. */
				static UInt getMaxCharge () 
				{ 
					return (max_charge_); 
				}			
			
				/** @brief Sets the @p max_charge parameter. */
				static void setMaxCharge (const UInt max_charge) 
				{ 
					max_charge_ = max_charge; 
				}	
		
				/** @brief Returns the table_steps_ parameter. 
 					* 
 					* This is an internally used parameter controlling the precision of several pre-sampling steps. 
 					* Normally, this parameter can be left unchanged. */
				static DoubleReal getTableSteps () 
				{ 
					return (table_steps_); 
				}			 					

				/** @brief Returns the inv_table_steps_ parameter. 
 					* 
 					* This is an internally used parameter controlling the precision of several pre-sampling steps. 
 					* Normally, this parameter can be left unchanged. */
				static DoubleReal getInvTableSteps () 
				{ 
					return (inv_table_steps_); 
				}			
				
				/** @brief Sets the @p table_steps parameter. */ 
				static void setTableSteps (const DoubleReal table_steps) 
				{
					inv_table_steps_ = 1./table_steps;
					table_steps_ = table_steps; 
				}


				/** @brief Returns the mass-parameter lambda (linear fit). */ 
				static DoubleReal getLambdaL (const DoubleReal m) ;


				/** @brief Computes the averagine isotopic distribution we would expect at the deconvoluted mass. 
 					* @param m The deconvoluted mass m.	
 					* @param size Returns the number of significant peaks within a pattern occurring at mass @p m.
 					* @return The isotopic distribution. */ 
				static const IsotopeDistribution::ContainerType& getAveragine (const DoubleReal m, UInt* size=NULL) ;


				/** @brief Returns the largest possible index for the pre-sampled gamma table. */
				static Size getGammaTableMaxIndex ()
				{
					return (gamma_table_max_index_);
				};

				/** @brief Returns the largest possible index for the pre-sampled exp table. */
				static Size getExpTableMaxIndex ()
				{
					return (exp_table_max_index_);
				};
						

				/** @brief Internally used function; uses register shifts for fast computation of the power function. 
					* @note Please, do not modify this function. */
				static float myPow (float a, float b) ;		
			
				static UInt getMzPeakCutOffAtMonoPos (const DoubleReal mass, const UInt z);
	
				static UInt getNumPeakCutOff (const DoubleReal mass, const UInt z);
				
				static UInt getNumPeakCutOff (const DoubleReal mz);


		protected:

				/** The singleton pointer. */
				static IsotopeWavelet* me_; 

				/** @brief Default Constructor. */
				IsotopeWavelet () ;				

				/** @brief Constructor 
 					* @param max_m The maximal deconvoluted mass that occurs in the current data set. 
 					* @param max_charge The maximal charge state we would like to analyze. */
				IsotopeWavelet (const DoubleReal max_m, const UInt max_charge) ;


				/** @brief Destructor. */
				virtual ~IsotopeWavelet () ;
				

				/** @brief Should be called once before values are drawn from the isotope wavelet function.
 					* The function is automatically called by the public constructor.   
					*
					* The function pre-computes the expensive gamma function. Parameters related to this function are:
					* @see max_charge_ and @see peak_cutoff_. If both of these are set correctly @see getValue will never compute
					* the gamma function online. 
					* 
					* @param max_m The maximal deconvoluted mass that occurs in the current data set. */
				static void preComputeExpensiveFunctions_ (const DoubleReal max_m) ;


				/** @brief Initializes the internally used averagine model; automatically called by the public constructor.
 					* @param max_m The maximal deconvoluted mass that occurs in the current data set.	*/ 
				static void computeIsotopeDistributionSize_ (const DoubleReal max_m) ;

				
				/** @brief Internally used function; uses register shifts for fast computation of the power function. 
					*	The function follows http://www.dctsystems.co.uk/Software/power.html , code by 
					* Ian Stephenson, DCT Systems, NCCA Bournemouth University.
					* @note Please, do not modify this function. */
				static float myPow2_ (float i) ;
		
				/** @brief Internally used function uses register shifts for fast computation of the power function. 
					* The function follows http://www.dctsystems.co.uk/Software/power.html , code by
					* Ian Stephenson, DCT Systems, NCCA Bournemouth University.
					* @note Please, do not modify this function. */
				static float myLog2_ (float i) ;

				/** @brief Internal union for fast computation of the power function. */
				union fi_
				{
					Int i;
					float f;
				};
			
				/** This parameter determines the maximal charge state we will consider. */
				static UInt max_charge_; 				

				/** This parameter determines the sample rate for the pre-computation of the gamma function. */
				static DoubleReal table_steps_;
				static DoubleReal inv_table_steps_;

				/** Internal table for the precomputed values of the gamma function. */ 
				static std::vector<DoubleReal> gamma_table_;
				static std::vector<DoubleReal> gamma_table_new_;
				
				/** Internal table for the precomputed values of the exponential function. */ 
				static std::vector<DoubleReal> exp_table_;

				/** Internal table for the precomputed values of the exponential function. */ 
				static std::vector<DoubleReal> sine_table_;				

				/** Internally used averagine model. */
				static IsotopeDistribution averagine_;

				static Size gamma_table_max_index_;
				static Size exp_table_max_index_;

	};

} //namespace

#endif 
