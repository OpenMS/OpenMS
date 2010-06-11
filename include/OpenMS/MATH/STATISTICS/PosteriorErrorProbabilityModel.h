// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_MATH_STATISTICS_POSTERIORERRORPROBABILITYMODEL_H
#define OPENMS_MATH_STATISTICS_POSTERIORERRORPROBABILITYMODEL_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <vector>

namespace OpenMS
{
	namespace Math
	{
	  
	  
	  /** 
	  	@brief Implements a mixture model for the inverse gumbel and the gauss distribution.
	
	    This class fits a Gumbel distribution and a Gauss distribution to a set of data points.
	    The computed probabilities can be 
			GumbelDistributionFitResult and GaussFitResult store the parameter for the distributions.
		
	    The formula with the fitted parameters can be transformed into a
	    gnuplot formula using getGumbelGnuplotFormula() and getGaussGnuplotFormula() after fitting.
	
			
			@ingroup Math
		*/
		class OPENMS_DLLAPI PosteriorErrorProbabilityModel: public DefaultParamHandler
		{
			public:
	
		
				///Default constructor
				PosteriorErrorProbabilityModel();

				///Destructor
				virtual ~PosteriorErrorProbabilityModel(); 
				
				/**
					@brief fits the distributions to the data points(x_scores) and writes the computed probabilites into the given vector (the second one).
					@param x_scores a vector which holds the data points
					@param probabilities a vector which holds the probability for each data point after running this function. If it has some content it will be overwritten.
					@exception Exception::UnableToFit is thrown if fitting cannot be performed
				*/
				void fit(std::vector<double>& x_scores, std::vector<double>& probabilities);		
				
				
				///returns estimated gauss parameters. Fit should be used before.
				GaussFitter::GaussFitResult getGaussFitResult() const
				{
					return gauss_fit_param_;
				}
				
				///returns estimated gumbel parameters. Fit should be used before.
				GumbelDistributionFitter::GumbelDistributionFitResult getGumbelFitResult() const
				{
					return gumbel_fit_param_;
				}
				
				///returns the estimated negative prior.
				DoubleReal getNegativePrior() const
				{
					return negative_prior_;
				}	
				
				/// returns the gnuplot formula of the fitted gumbel distribution.
				const String getGumbelGnuplotFormula() const;
				
				/// returns the gnuplot formula of the fitted gauss distribution.
				const String getGaussGnuplotFormula() const;
				
				/// returns the gnuplot formula of the fitted mixture distribution.
				const String getBothGnuplotFormula() const;
				
				/**
					Returns the computed posterior error probability for a given score.
					@note: fit has to be used before this function. Otherwise this function will compute nonsense.
				*/
				DoubleReal computeProbability(DoubleReal score);
				
			private:
				/// assignment operator (not implemented)
				PosteriorErrorProbabilityModel& operator = (const PosteriorErrorProbabilityModel& rhs);
				///Copy constructor (not implemented)
				PosteriorErrorProbabilityModel(const PosteriorErrorProbabilityModel& rhs);
				///stores gumbel parameters
				GumbelDistributionFitter::GumbelDistributionFitResult gumbel_fit_param_;
				///stores gauss parameters
				GaussFitter::GaussFitResult	gauss_fit_param_;
				///stores final prior probability for negative peptides
				DoubleReal negative_prior_;
				///peak of the gumbel distribution
				DoubleReal max_gumbel_;
				///peak of the gauss distribution
				DoubleReal max_gauss_;			
			  ///smallest score which was used for fitting the model
				DoubleReal smallest_score_;
		};
	}
}

#endif

