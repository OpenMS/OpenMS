// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

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
      @note All parameters are stored in GaussFitResult. In the case of the Gumbel distribution x0 and sigma represent the local parameter alpha and the scale parameter beta, respectively.

      @htmlinclude OpenMS_Math::PosteriorErrorProbabilityModel.parameters

      @ingroup Math
    */
    class OPENMS_DLLAPI PosteriorErrorProbabilityModel :
      public DefaultParamHandler
    {
public:

      ///default constructor
      PosteriorErrorProbabilityModel();

      ///Destructor
      ~PosteriorErrorProbabilityModel() override;

      /**
          @brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables. computeProbability can be used afterwards.
          @param search_engine_scores a vector which holds the data points
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vector is sorted from smallest to biggest value!
      */
      bool fit(std::vector<double> & search_engine_scores);

      /**
          @brief fits the distributions to the data points(search_engine_scores) and writes the computed probabilities into the given vector (the second one).
          @param search_engine_scores a vector which holds the data points
          @param probabilities a vector which holds the probability for each data point after running this function. If it has some content it will be overwritten.
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vectors are sorted from smallest to biggest value!
      */
      bool fit(std::vector<double> & search_engine_scores, std::vector<double> & probabilities);

      ///Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillDensities(std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///computes the Maximum Likelihood with a log-likelihood function.
      double computeMaxLikelihood(std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///sums (1 - posterior probabilities)
      double one_minus_sum_post(std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///sums  posterior probabilities
      double sum_post(std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///helper function for the EM algorithm (for fitting)
      double sum_pos_x0(std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///helper function for the EM algorithm (for fitting)
      double sum_neg_x0(std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///helper function for the EM algorithm (for fitting)
      double sum_pos_sigma(std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density, double positive_mean);
      ///helper function for the EM algorithm (for fitting)
      double sum_neg_sigma(std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density, double positive_mean);


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
      double getNegativePrior() const
      {
        return negative_prior_;
      }

      ///computes the gaussian density at position x with parameters params.
      double getGauss(double x, const GaussFitter::GaussFitResult & params)
      {
        return params.A * exp(-1.0 * pow(x - params.x0, 2) / (2 * pow(params.sigma, 2)));
      }

      ///computes the gumbel density at position x with parameters params.
      double getGumbel(double x, const GaussFitter::GaussFitResult & params)
      {
        double z = exp((params.x0 - x) / params.sigma);
        return (z * exp(-1 * z)) / params.sigma;
      }

      /**
          Returns the computed posterior error probability for a given score.
          @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
      */
      double computeProbability(double score);

      ///initializes the plots
      TextFile initPlots(std::vector<double> & x_scores);

      /// returns the gnuplot formula of the fitted gumbel distribution. Only x0 and sigma are used as local parameter alpha and scale parameter beta, respectively.
      const String getGumbelGnuplotFormula(const GaussFitter::GaussFitResult & params) const;

      /// returns the gnuplot formula of the fitted gauss distribution.
      const String getGaussGnuplotFormula(const GaussFitter::GaussFitResult & params) const;

      /// returns the gnuplot formula of the fitted mixture distribution.
      const String getBothGnuplotFormula(const GaussFitter::GaussFitResult & incorrect, const GaussFitter::GaussFitResult & correct) const;

      ///plots the estimated distribution against target and decoy hits
      void plotTargetDecoyEstimation(std::vector<double> & target, std::vector<double> & decoy);

      /// returns the smallest score used in the last fit
      inline double getSmallestScore()
      {
        return smallest_score_;
      }

      /// try to invoke 'gnuplot' on the file to create PDF automatically
      void tryGnuplot(const String& gp_file);

private:
      /// assignment operator (not implemented)
      PosteriorErrorProbabilityModel & operator=(const PosteriorErrorProbabilityModel & rhs);
      ///Copy constructor (not implemented)
      PosteriorErrorProbabilityModel(const PosteriorErrorProbabilityModel & rhs);
      ///stores parameters for incorrectly assigned sequences. If gumbel fit was used, A can be ignored. Furthermore, in this case, x0 and sigma are the local parameter alpha and scale parameter beta, respectively.
      GaussFitter::GaussFitResult incorrectly_assigned_fit_param_;
      ///stores gauss parameters
      GaussFitter::GaussFitResult correctly_assigned_fit_param_;
      ///stores final prior probability for negative peptides
      double negative_prior_;
      ///peak of the incorrectly assigned sequences distribution
      double max_incorrectly_;
      ///peak of the gauss distribution (correctly assigned sequences)
      double max_correctly_;
      ///smallest score which was used for fitting the model
      double smallest_score_;
      ///points to getGauss
      double (PosteriorErrorProbabilityModel::* calc_incorrect_)(double x, const GaussFitter::GaussFitResult & params);
      ///points either to getGumbel or getGauss depending on whether one uses the gumbel or the gaussian distribution for incorrectly assigned sequences.
      double (PosteriorErrorProbabilityModel::* calc_correct_)(double x, const GaussFitter::GaussFitResult & params);
      ///points either to getGumbelGnuplotFormula or getGaussGnuplotFormula depending on whether one uses the gumbel or the gaussian distribution for incorrectly assigned sequences.
      const String (PosteriorErrorProbabilityModel::* getNegativeGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;
      ///points to getGumbelGnuplotFormula
      const String (PosteriorErrorProbabilityModel::* getPositiveGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;

    };
  }
}

#endif // OPENMS_MATH_STATISTICS_POSTERIORERRORPROBABILITYMODEL_H
