// -*- mode: C++; tab-width: 2; -*-
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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <math.h>
#include <iostream>

#ifndef ONEOLOG2E
#define ONEOLOG2E 0.6931471806 
#endif

namespace OpenMS
{
	//internally used variables / defaults
	IsotopeWavelet* IsotopeWavelet::me_ = NULL;
	UInt IsotopeWavelet::max_charge_ = 1;
	std::vector<DoubleReal> IsotopeWavelet::gamma_table_;
	std::vector<DoubleReal> IsotopeWavelet::exp_table_;
	DoubleReal IsotopeWavelet::table_steps_ = 0.001;
	DoubleReal IsotopeWavelet::inv_table_steps_ = 1./table_steps_;
	IsotopeDistribution IsotopeWavelet::averagine_;
	Int IsotopeWavelet::gamma_table_max_index_ = -1;
	Int IsotopeWavelet::exp_table_max_index_ = -1;


	IsotopeWavelet* IsotopeWavelet::init (const DoubleReal max_m, const UInt max_charge) 
	{
		delete (me_); //either me_ is NULL or is already instantiated
		me_ = new IsotopeWavelet (max_m, max_charge);
		
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
		max_charge_ = 1;
		table_steps_ = 0.001;
		inv_table_steps_ = 1./table_steps_;
		me_ = NULL;
	}
		

	DoubleReal IsotopeWavelet::getValueByLambda (const DoubleReal lambda, const DoubleReal tz1) 
	{
		DoubleReal tz = tz1-1;
		DoubleReal fi_lgamma (gamma_table_ [(int)(tz1*inv_table_steps_)]);
		
		DoubleReal fac (-lambda + tz*myLog2_(lambda)*ONEOLOG2E - fi_lgamma);

		return (sin(tz*WAVELET_PERIODICITY) * exp(fac));
	}
	

	DoubleReal IsotopeWavelet::getValueByLambdaExtrapol (const DoubleReal lambda, const DoubleReal tz1) 
	{
		DoubleReal fac (-lambda + (tz1-1)*myLog2_(lambda)*ONEOLOG2E - lgamma(tz1));
		
		return (sin((tz1-1)*WAVELET_PERIODICITY) * exp(fac));
	}


	DoubleReal IsotopeWavelet::getLambdaL (const DoubleReal m) 
	{
		return (LAMBDA_L_0 + LAMBDA_L_1*m);
	}
				
	DoubleReal IsotopeWavelet::getLambdaQ (const DoubleReal m) 
	{
		return (LAMBDA_Q_0 + LAMBDA_Q_1*m + LAMBDA_Q_2*m*m);
	}
			
	float IsotopeWavelet::myPow (float a, float b) 		
	{	
		float help (b*myLog2_(a));
		return ( (help<127) ? myPow2_(help) : pow(2, help) ); 
	}


	/** The upcoming code follows the ideas from Ian Stephenson, DCT Systems, NCCA Bournemouth University.
		* See also: http://www.dctsystems.co.uk/Software/power.html */ 
	float IsotopeWavelet::myPow2_ (float i) 		
	{	
		float y=i-(int)i;
		y=(y-y*y)*POW_CONST;
		float x=i+127-y;
		x*=SHIFT23;
		fi_ z;
		z.i=(Int) x;
		return (z.f);
	}


	float IsotopeWavelet::myLog2_ (float i) 
	{	
		fi_ x;
		x.f=i;
		float x2 =x.i;
		x2*= SHIFT23_00;
		x2-=127; 
		float y=x2-(int)x2;
		y=(y-y*y)*LOG_CONST;
		return (x2+y);
	}

	void IsotopeWavelet::preComputeExpensiveFunctions_ (const DoubleReal max_m) 
	{
		UInt peak_cutoff;
		IsotopeWavelet::getAveragine (max_m*max_charge_, &peak_cutoff);
		++peak_cutoff; //just to be sure, since getPeakCutOff (see IsotopeWaveletTransform.h) can return slightly different values 
		//This would be the theoretically justified way to estimate the boundary ...
		//UInt up_to = (UInt) ceil(max_charge_ * (peak_cutoff+QUARTER_NEUTRON_MASS) + 1);
		//... but in practise, it pays off to sample some points more.
		UInt up_to=2*(peak_cutoff*max_charge_+1);
		gamma_table_.clear();
		exp_table_.clear();
		DoubleReal query=0; 
		while (query <= up_to)
		{
			//std::cout << log(1./tgamma(query)) << "\t" << -lgamma(query) << std::endl;

			//gamma_table_.push_back(1./tgamma(query));
			gamma_table_.push_back (lgamma(query));
			query += table_steps_;	
		};	
		gamma_table_max_index_ = gamma_table_.size();

		DoubleReal up_to2 = getLambdaQ(max_m*max_charge_);
		query=0;
		while (query <= up_to2)
		{				
			exp_table_.push_back(exp(-query));
			query += table_steps_;	
		};
		exp_table_max_index_ = exp_table_.size();
	}
											

	const IsotopeDistribution::ContainerType& IsotopeWavelet::getAveragine (const DoubleReal mass, UInt* size) 
	{
		averagine_.estimateFromPeptideWeight (mass);
		IsotopeDistribution::ContainerType help (averagine_.getContainer());	
		IsotopeDistribution::ContainerType::iterator iter;
		
		if (size != NULL)
		{
			UInt count=help.size();
			for (iter=help.end()-1; iter!=help.begin(); --iter, --count)
			{	
				//maybe we should provide some interface to that constant, although its range is rather limited and
				//its influence within this range is negligible.
				if (iter->second >= 0.05)
					break;
			};
			*size=count;
		}; 

		return (averagine_.getContainer());
	}


	void IsotopeWavelet::computeIsotopeDistributionSize_ (const DoubleReal max_m) 
	{
		DoubleReal max_deconv_mz = max_m*max_charge_;
		averagine_.setMaxIsotope(UInt(max_deconv_mz/100.+10.)); // expect less than 10 extra Da for heavy isotopes per 1000 Da mono mass.
		// averagine_.setMaxIsotope (INT_MAX); // old version INT_MAX is C not C++, should use #include <limits> anyway
		averagine_.estimateFromPeptideWeight (max_deconv_mz);
		//maybe we should provide some interface to that constant, although its range is rather limited and
		//its influence within this range is negligible.
		averagine_.trimRight (0.05); 
		averagine_.setMaxIsotope (averagine_.getContainer().size());
	}


} //namespace
