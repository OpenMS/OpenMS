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
	class String;
	class TextFile;
	namespace Math
	{
	  
	  
	  /** 
	  	@brief Implements a mixture model of the inverse gumbel and the gauss distribution or a gaussian mixture.
	
	    This class fits either a Gumbel distribution and a Gauss distribution to a set of data points or two Gaussian distributions using the EM algorithm.
	    One can output the fit as a gnuplot formula using getGumbelGnuplotFormula() and getGaussGnuplotFormula() after fitting.
			@Note All paaremters are stored in GaussFitResult. In the case of the gumbel distribution x0 and sigma represent the local parameter alpha and the scale parameter beta, respectively.
			
			@ingroup Math
		*/
		class OPENMS_DLLAPI PosteriorErrorProbabilityModel: public DefaultParamHandler
		{
			public:
			
				///default constructor
				PosteriorErrorProbabilityModel();
				
				///Destructor
				virtual ~PosteriorErrorProbabilityModel(); 

				/**
					@brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables. computeProbability can be used afterwards.
					@param search_engine_scores a vector which holds the data points
					@exception Exception::UnableToFit is thrown if fitting cannot be performed
					@note the vector is sorted from smallest to biggest value!
				*/
				void fit( std::vector<double>& search_engine_scores);	
				
				/**
					@brief fits the distributions to the data points(search_engine_scores) and writes the computed probabilites into the given vector (the second one).
					@param search_engine_scores a vector which holds the data points
					@param probabilities a vector which holds the probability for each data point after running this function. If it has some content it will be overwritten.
					@exception Exception::UnableToFit is thrown if fitting cannot be performed
					@note the vectors are sorted from smallest to biggest value!
				*/
				void fit( std::vector<double>& search_engine_scores, std::vector<double>& probabilities);				
				
				///Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorreclty assigned seqeuences.
				void fillDensities(std::vector<double>& x_scores,std::vector<DoubleReal>& incorrect_density,std::vector<DoubleReal>& correct_density);
				///computes the Maximum Likelihood with a log-likelihood funciotn.
				DoubleReal computeMaxLikelihood(std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density);
				///sums (1 - posterior porbabilities)
				DoubleReal one_minus_sum_post(std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density);
				///sums  posterior porbabilities
				DoubleReal sum_post(std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density);
				///helper function for the EM algorithm (for fitting)
				DoubleReal sum_pos_x0(std::vector<double>& x_scores, std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density);
				///helper function for the EM algorithm (for fitting)
				DoubleReal sum_neg_x0(std::vector<double>& x_scores, std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density);
				///helper function for the EM algorithm (for fitting)
				DoubleReal sum_pos_sigma(std::vector<double>& x_scores, std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density, DoubleReal positive_mean);
				///helper function for the EM algorithm (for fitting)
				DoubleReal sum_neg_sigma(std::vector<double>& x_scores, std::vector<DoubleReal>& incorrect_density, std::vector<DoubleReal>& correct_density, DoubleReal positive_mean);

				
				///returns estimated parameters for correctly assigned sequences. Fit should be used before.
				GaussFitter::GaussFitResult getCorrectlyAssignedFitResult() const
				{
					return correctly_assigned_fit_param_;
				}
				
				///returns estimated parameters for correctly assigned sequences. Fit should be used before.
				GaussFitter::GaussFitResult getIncorrectlyAssignedFitResult() const
				{
					return incorrectly_assigned_fit_param_;
				}
				
				///returns the estimated negative prior probability.
				DoubleReal getNegativePrior() const
				{
					return negative_prior_;
				}	
				///computes the gaussian density at position x with parameters params.
				DoubleReal getGauss(DoubleReal x, const GaussFitter::GaussFitResult& params)
				{
					return params.A * exp(-1.0 * pow(x - params.x0, 2) / (2 * pow(params.sigma, 2)));
				}
				///computes the gumbel density at position x with parameters params.
				DoubleReal getGumbel(DoubleReal x, const GaussFitter::GaussFitResult& params)
				{
					DoubleReal z = exp((params.x0 - x)/params.sigma);
					return (z*exp(-1* z))/params.sigma;		
				}				
				
				/**
					Returns the computed posterior error probability for a given score.
					@note: fit has to be used before using this function. Otherwise this function will compute nonsense.
				*/
				DoubleReal computeProbability(DoubleReal score);
				
				///initializes the plots
				TextFile* InitPlots(std::vector<double> & x_scores);
				
				/// returns the gnuplot formula of the fitted gumbel distribution. Only x0 and sigma are used as local parameter alpha and scale parameter beta, respectively. 
				const String getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const;
				
				/// returns the gnuplot formula of the fitted gauss distribution.
				const String getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const;

				/// returns the gnuplot formula of the fitted mixture distribution.
				const String getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect, const GaussFitter::GaussFitResult& correct) const;				
			private:
				/// assignment operator (not implemented)
				PosteriorErrorProbabilityModel& operator = (const PosteriorErrorProbabilityModel& rhs);
				///Copy constructor (not implemented)
				PosteriorErrorProbabilityModel(const PosteriorErrorProbabilityModel& rhs);
				///stores parameters for incorrectly assigned sequences. If gumbel fit was used, A can be ignored. Furthermore, in this case, x0 and sigma are the local parameter alpha and scale parameter beta, respectively.
				GaussFitter::GaussFitResult	incorrectly_assigned_fit_param_;
				///stores gauss parameters
				GaussFitter::GaussFitResult	correctly_assigned_fit_param_;
				///stores final prior probability for negative peptides
				DoubleReal negative_prior_;
				///peak of the incorrectly assigned sequences distribution
				DoubleReal max_incorrectly_;
				///peak of the gauss distribution (correctly assigned sequences)
				DoubleReal max_correctly_;			
			  ///smallest score which was used for fitting the model
				DoubleReal smallest_score_;
				///points to getGauss
				DoubleReal (PosteriorErrorProbabilityModel::*calc_incorrect_)(DoubleReal x,const GaussFitter::GaussFitResult& params);
				///points either to getGumbel or getGauss depending on whether on uses the gumbel or th gausian distribution for incorrectly assigned sequences.
				DoubleReal (PosteriorErrorProbabilityModel::*calc_correct_)(DoubleReal x,const GaussFitter::GaussFitResult& params);
				///points either to getGumbelGnuplotFormula or getGaussGnuplotFormula depending on whether on uses the gumbel or th gausian distribution for incorrectly assigned sequences.
				const String (PosteriorErrorProbabilityModel::*getNegativeGnuplotFormula_)(const GaussFitter::GaussFitResult& params) const;
				///points to getGumbelGnuplotFormula
				const String (PosteriorErrorProbabilityModel::*getPositiveGnuplotFormula_)(const GaussFitter::GaussFitResult& params) const;
				
		};
	}
}

#endif

