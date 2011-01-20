// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>

#ifndef ONEOLOG2E
#define ONEOLOG2E 0.6931471806 
#endif

#ifndef TWOPI
#define TWOPI 6.283185307
#endif

#ifndef SINETRAINVECS
#define SINETRAINVECS
		double avec [] =
		{
			1.00729092576906387,
			8.47021314902288758,
			0.568589043446217191,
			0.155624756155312988,
			0.0546681591351188961,
			0.0239278712518182744,
			0.0113785215738483325,
			0.00521306547958902237,
			0.00215635823841069541,
			0.000803390368142344806
		};

		double bvec [] =
		{
			0.999413854502327892,
			0.999402216792031672,
			0.999781419604179100,
			1.00003385169685144,
			1.00019089578658660,
			1.00027323765134879,
			1.00032322673781815,
			1.00036799996074044,
			1.00041798239640700,
			1.00047277354586983
		};
		
		double cvec [] =
		{
			-1.82357924751654536e-7,
			0.0000617158396180752893,
			0.275538540246282853e-3,
			0.000340752172835340433,
			0.366877055752254583e-3,
			0.000369558360963110233,
			0.365734273818181337e-3,
			0.000363964822444512546,
			0.366472952708349505e-3,
			0.372090375255612713e-3
		};
 		
		double dvec [] =
		{
			1.00134499727445,
			0.000655514660572,
			-0.5396047396523e-1,
			-0.15108542328434,
			-0.29337450695613,
			-0.48557863368653,
			-.70449316758845,
			-.93081748650931,
			-1.15609867020231,
			-1.37929313525534
		};
#endif

namespace OpenMS
{
	//internally used variables / defaults
	IsotopeWavelet* IsotopeWavelet::me_ = NULL;
	UInt IsotopeWavelet::max_charge_ = 1;
	std::vector<DoubleReal> IsotopeWavelet::gamma_table_;
	std::vector<DoubleReal> IsotopeWavelet::gamma_table_new_;
	std::vector<DoubleReal> IsotopeWavelet::exp_table_;
	std::vector<DoubleReal> IsotopeWavelet::sine_table_;
	DoubleReal IsotopeWavelet::table_steps_ = 0.0001;
	DoubleReal IsotopeWavelet::inv_table_steps_ = 1./table_steps_;
	IsotopeDistribution IsotopeWavelet::averagine_;
	Size IsotopeWavelet::gamma_table_max_index_ = 0;
	Size IsotopeWavelet::exp_table_max_index_ = 0;


	IsotopeWavelet* IsotopeWavelet::init (const DoubleReal max_m, const UInt max_charge) 
	{
		if (me_ == NULL)
		{	
			me_ = new IsotopeWavelet (max_m, max_charge);
		};
		
		return (me_);
	}
	

	IsotopeWavelet::IsotopeWavelet () 
	{ 
	}
				
	IsotopeWavelet::IsotopeWavelet (const DoubleReal max_m, const UInt max_charge) 
	{
		max_charge_ = max_charge;
		computeIsotopeDistributionSize_ (max_m);
		preComputeExpensiveFunctions_ (max_m);
	}
	
	IsotopeWavelet::~IsotopeWavelet () 
	{ 
	}
	
	void IsotopeWavelet::destroy () 
	{
		delete (me_);
		me_ = NULL;
		max_charge_ = 1;
		gamma_table_.clear();
		exp_table_.clear();
		sine_table_.clear();
		table_steps_ = 0.0001;
		inv_table_steps_ = 1./table_steps_;
		gamma_table_max_index_ = 0;
		exp_table_max_index_ = 0;
	} 	

	DoubleReal IsotopeWavelet::getValueByLambda (const DoubleReal lambda, const DoubleReal tz1) 
	{
		DoubleReal tz (tz1-1);
		DoubleReal fi_lgamma (gamma_table_ [(Int)(tz1*inv_table_steps_)]);
		DoubleReal help (tz*Constants::WAVELET_PERIODICITY/(TWOPI));
		DoubleReal sine_index ((help-(int)(help))*TWOPI*inv_table_steps_);
		DoubleReal fac (-lambda + tz*myLog2_(lambda)*ONEOLOG2E - fi_lgamma);

		return (sine_table_[(Int)(sine_index)] * exp(fac));
	}
	

	DoubleReal IsotopeWavelet::getValueByLambdaExtrapol (const DoubleReal lambda, const DoubleReal tz1) 
	{
		DoubleReal fac (-lambda + (tz1-1)*myLog2_(lambda)*ONEOLOG2E - boost::math::lgamma(tz1));
		DoubleReal help ((tz1-1)*Constants::WAVELET_PERIODICITY/(TWOPI));
		DoubleReal sine_index ((help-(int)(help))*TWOPI*inv_table_steps_);
		
		return (sine_table_[(Int)(sine_index)] * exp(fac));
	}
	
	DoubleReal IsotopeWavelet::getValueByLambdaExact (const DoubleReal lambda, const DoubleReal tz1) 
	{
    return (sin(2*Constants::PI*(tz1-1)/Constants::IW_NEUTRON_MASS)*exp(-lambda)*pow(lambda, tz1-1)/boost::math::tgamma(tz1));//boost::math::tgamma(tz1));
	}

	DoubleReal IsotopeWavelet::getLambdaL (const DoubleReal m) 
	{
		return (Constants::LAMBDA_L_0 + Constants::LAMBDA_L_1*m);
	}
				

	UInt IsotopeWavelet::getMzPeakCutOffAtMonoPos (const DoubleReal mass, const UInt z)
	{ 	
		DoubleReal mz (mass*z);
		int res=-1;
		if (mz<Constants::CUT_LAMBDA_BREAK_0_1)
			res = ceil(Constants::CUT_LAMBDA_Q_0_A+Constants::CUT_LAMBDA_Q_0_B*mz+Constants::CUT_LAMBDA_Q_0_C*mz*mz);
		if (mz>Constants::CUT_LAMBDA_BREAK_1_2)
			res =ceil(Constants::CUT_LAMBDA_L_2_A+Constants::CUT_LAMBDA_L_2_B*mz);
		if (res < 0)
			res = ceil(Constants::CUT_LAMBDA_Q_1_A+Constants::CUT_LAMBDA_Q_1_B*mz+Constants::CUT_LAMBDA_Q_1_C*mz*mz);

		return (res);
	}

	UInt IsotopeWavelet::getNumPeakCutOff (const DoubleReal mass, const UInt z)
	{ 
		DoubleReal mz (mass*z);
		int res=-1;
		if (mz<Constants::CUT_LAMBDA_BREAK_0_1)
			res = ceil(Constants::CUT_LAMBDA_Q_0_A+Constants::CUT_LAMBDA_Q_0_B*mz+Constants::CUT_LAMBDA_Q_0_C*mz*mz-Constants::IW_QUARTER_NEUTRON_MASS);
		if (mz>Constants::CUT_LAMBDA_BREAK_1_2)
			res = ceil(Constants::CUT_LAMBDA_L_2_A+Constants::CUT_LAMBDA_L_2_B*mz-Constants::IW_QUARTER_NEUTRON_MASS);
		if (res < 0)
			res = ceil(Constants::CUT_LAMBDA_Q_1_A+Constants::CUT_LAMBDA_Q_1_B*mz+Constants::CUT_LAMBDA_Q_1_C*mz*mz-Constants::IW_QUARTER_NEUTRON_MASS);

		return (res);
	}	
	
	UInt IsotopeWavelet::getNumPeakCutOff (const DoubleReal mz)
	{ 
		int res=-1;
		if (mz<Constants::CUT_LAMBDA_BREAK_0_1)
			res = ceil(Constants::CUT_LAMBDA_Q_0_A+Constants::CUT_LAMBDA_Q_0_B*mz+Constants::CUT_LAMBDA_Q_0_C*mz*mz-Constants::IW_QUARTER_NEUTRON_MASS);
		if (mz>Constants::CUT_LAMBDA_BREAK_1_2)
			res = ceil(Constants::CUT_LAMBDA_L_2_A+Constants::CUT_LAMBDA_L_2_B*mz-Constants::IW_QUARTER_NEUTRON_MASS);
		if (res < 0)
			res = ceil(Constants::CUT_LAMBDA_Q_1_A+Constants::CUT_LAMBDA_Q_1_B*mz+Constants::CUT_LAMBDA_Q_1_C*mz*mz-Constants::IW_QUARTER_NEUTRON_MASS);

		return (res);
	}

		
	float IsotopeWavelet::myPow (float a, float b) 		
	{	
		float help (b*myLog2_(a));
		return ( (help>0 && help<127) ? myPow2_(help) : pow(2, help) ); 
	}


	/** The upcoming code follows the ideas from Ian Stephenson, DCT Systems, NCCA Bournemouth University.
		* See also: http://www.dctsystems.co.uk/Software/power.html */ 
	float IsotopeWavelet::myPow2_ (float i) 		
	{	
		float y=i-(int)i;
		y=(y-y*y)*Constants::POW_CONST;
		float x=i+127-y;
		x*=Constants::SHIFT23;
		fi_ z;
		z.i=(Int) x;
		return (z.f);
	}


	float IsotopeWavelet::myLog2_ (float i) 
	{	
		fi_ x;
		x.f=i;
		float x2 =x.i;
		x2*= Constants::SHIFT23_00;
		x2-=127; 
		float y=x2-(int)x2;
		y=(y-y*y)*Constants::LOG_CONST;
		return (x2+y);
	}

	void IsotopeWavelet::preComputeExpensiveFunctions_ (const DoubleReal max_m) 
	{
		UInt peak_cutoff = getNumPeakCutOff(max_m, max_charge_);
		UInt up_to = peak_cutoff*max_charge_+1;
		gamma_table_.clear();
		gamma_table_new_.clear();
		exp_table_.clear();
		DoubleReal query=0;
		gamma_table_.push_back (std::numeric_limits<int>::max());
		gamma_table_new_.push_back (std::numeric_limits<int>::max());
		query += table_steps_; 
		while (query <= up_to)
		{
			gamma_table_.push_back (boost::math::lgamma(query));
			gamma_table_new_.push_back (boost::math::tgamma(query));

			query += table_steps_;	
		};	
		gamma_table_max_index_ = gamma_table_.size();

		DoubleReal up_to2 = getLambdaL(max_m*max_charge_);
		query=0;
		while (query <= up_to2)
		{				
			exp_table_.push_back(exp(-query));
			query += table_steps_;	
		};
		exp_table_max_index_ = exp_table_.size();

		query=0;
		while (query < 2*Constants::PI)
		{
			sine_table_.push_back (sin(query));
			query += table_steps_;
		};
	}
											

	const IsotopeDistribution::ContainerType& IsotopeWavelet::getAveragine (const DoubleReal mass, UInt* size) 
	{
	
		averagine_.estimateFromPeptideWeight (mass);
		IsotopeDistribution::ContainerType help (averagine_.getContainer());	
		IsotopeDistribution::ContainerType::iterator iter;
		
		if (size != NULL)
		{
			*size = getNumPeakCutOff(mass);
		}; 

		return (averagine_.getContainer());
	}


	void IsotopeWavelet::computeIsotopeDistributionSize_ (const DoubleReal max_m) 
	{
		DoubleReal max_deconv_mz = max_m*max_charge_;
		averagine_.setMaxIsotope(UInt(max_deconv_mz/100.+10.)); // expect less than 10 extra Da for heavy isotopes per 1000 Da mono mass.
		// averagine_.setMaxIsotope (INT_MAX); // old version INT_MAX is C not C++, should use #include <limits> anyway
		averagine_.estimateFromPeptideWeight (max_deconv_mz);
		Int max_isotope = getNumPeakCutOff(max_deconv_mz);
		averagine_.setMaxIsotope (max_isotope-1);
	}


} //namespace
