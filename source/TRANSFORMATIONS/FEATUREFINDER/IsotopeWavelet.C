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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <math.h>
#include <iostream>

namespace OpenMS
{

unsigned int IsotopeWavelet::peak_cutoff_ = 5, IsotopeWavelet::max_charge_ = 4; //defaults
std::vector<double> IsotopeWavelet::gamma_table_;
double IsotopeWavelet::gamma_steps_ = 0.001;


/** Default Constructor */
IsotopeWavelet::IsotopeWavelet () throw()
{ }

/** Destructor */
IsotopeWavelet::~IsotopeWavelet () throw()
{ }
	

/** @brief Returns the value of the isotope wavelet at position t.
	* Note that this functions returns the pure function value of the psi and not the normalized (average=0)
	* value given bei Psi. 
	* @param t The position at which the wavelet has to be drawn (within the coordinate system of the wavelet). 
	* @param m The mass (decharged!) within the signal.
	* @param z The charge z we want to detect. */
double IsotopeWavelet::getValueByMass (const double t, const double m, const unsigned int z, const int mode) throw ()
{
	if (t>peak_cutoff_+NEUTRON_MASS/4.)	
		return(0);
	
	int x0, x1; double f0, f1, fi;
	x0 = (int) trunc ((t*z+1)/gamma_steps_);
	x1 = x0+1;
	if (x1 < (int) gamma_table_.size())
	{
		f0 = gamma_table_[x0];
		f1 = gamma_table_[x1];
		fi = (f0 + (f1-f0)/((x1-x0)*gamma_steps_) * ((t*z+1)-x0*gamma_steps_));
	}
	else 
		fi = tgamma(t*z+1);

	double lambda= getLambdaL(m*z-z*mode*PROTON_MASS);		
	//return (sin(2*M_PI*t*z/NEUTRON_MASS) * exp(-lambda) * myPow(lambda,t*z) / tgamma(t*z+1));
	return (sin(2*M_PI*t*z/NEUTRON_MASS) * exp(-lambda) * myPow(lambda,t*z) / fi);
}


/** @brief Returns the value of the isotope wavelet at position t.
	* Note that this functions returns the pure function value of the psi and not the normalized (average=0)
	* value given bei Psi. 
	* @param t The position at which the wavelet has to be drawn (within the coordinate system of the wavelet). 
	* @param lambda The mass-parameter lambda.
	* @param z The charge z we want to detect. */
double IsotopeWavelet::getValueByLambda (const double t, const double lambda, const unsigned int z) throw ()
{
	if (t>peak_cutoff_+NEUTRON_MASS/4.)	
		return(0);
	
	int x0, x1; double f0, f1, fi;
	x0 = (int) trunc ((t*z+1)/gamma_steps_);
	x1 = x0+1;
	if (x1 < (int) gamma_table_.size())
	{
		f0 = gamma_table_[x0];
		f1 = gamma_table_[x1];
		fi = (f0 + (f1-f0)/((x1-x0)*gamma_steps_) * ((t*z+1)-x0*gamma_steps_));
	}
	else 
		fi = tgamma(t*z+1);

	//return (sin(2*M_PI*t*z/NEUTRON_MASS) * exp(-lambda) * myPow(lambda,t*z) / tgamma(t*z+1));
	return (sin(2*M_PI*t*z/NEUTRON_MASS) * exp(-lambda) * myPow(lambda,t*z) / fi);
}



/** @brief Returns the mass-dependet parameter lambda (linear fit). 
	* Note that the only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
double IsotopeWavelet::getLambdaL (const double m) throw ()
{
	return (LAMBDA_L_0 + LAMBDA_L_1*m);
}
			
/** @brief Returns the mass-dependet parameter lambda (quadratic fit). 
	* Note that the only possibility to switch between getLambdaL and LambdaQ is pure hardcoding. */
double IsotopeWavelet::getLambdaQ (const double m) throw ()
{
	return (LAMBDA_Q_0 + LAMBDA_Q_1*m + LAMBDA_Q_2*m*m);
}
		
#ifndef OPENMS_64BIT_ARCHITECTURE	
/** @brief Internal function using register shifts for fast computation of the power function. 
	* Please, do not modify this function. */
float IsotopeWavelet::myPow (float a, float b) throw ()		
{	
	return (myPow2(b*myLog2(a))); 
}
#else
float IsotopeWavelet::myPow (float a, float b) throw ()		
{
	std::cout << "normal power function" << std::endl;	
	return (pow(a,b)); 
}
#endif


#ifndef OPENMS_64BIT_ARCHITECTURE		
/** @brief Internal function using register shifts for fast computation of the power function. 
	* Please, do not modify this function. */
float IsotopeWavelet::myPow2 (float i) throw ()		
{	
  float y=i-floorf(i);
  y=(y-y*y)*PowBodge;
  float x=i+127-y;
	x*=shift23;
	fi z;
  z.i=(int) x;
  return (z.f);
}


/** @brief Internal function using register shifts for fast computation of the power function. 
	* Please, do not modify this function. */
float IsotopeWavelet::myLog2 (float i) throw ()
{	
	fi x;
	x.f=i;
  float x2=x.i;
  x2*= shift23_00;
  x2-=127; 
  float y=x2-floorf(x2);
  y=(y-y*y)*LogBodge;
  return (x2+y);
}
#endif
	

/** @brief Should be called once before values are drawn from the isotope wavelet function. 
 	*
	* The function precomputes the expensive gamma function. Parameters related to this function are:
 	* max_charge_ and peak_cutoff_. If both of these are set correctly getValue will never compute
 	* the gamme online. Please note that in a future and more efficient version checks for precomputed
 	* values will be removed. */
void IsotopeWavelet::preComputeGammaFunction () throw ()
{
	std::cout << "Precomputing the Gamma function ...";	
	unsigned int up_to = max_charge_ * peak_cutoff_ + 1;
	gamma_table_.clear();
	double query=0; 	
	while (query <= up_to)
	{
		gamma_table_.push_back(tgamma(query));
		query += gamma_steps_;	
	};
	std::cout << " ok." << std::endl;
}
		
} //namespace
