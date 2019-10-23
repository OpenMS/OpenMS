// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <vector>
#include <map>

namespace OpenMS
{
  class String;
  class TextFile;
  class PeptideIdentification;
  class ProteinIdentification;
  class PeptideHit;
  namespace Math
  {


    /**
      @brief Implements a mixture model of the inverse gumbel and the gauss distribution or a gaussian mixture.

      This class fits either a Gumbel distribution and a Gauss distribution to a set of data points or two Gaussian distributions using the EM algorithm.
      One can output the fit as a gnuplot formula using getGumbelGnuplotFormula() and getGaussGnuplotFormula() after fitting.
      @note All parameters are stored in GaussFitResult. In the case of the Gumbel distribution x0 and sigma represent the local parameter alpha and the scale parameter beta, respectively.

      @todo test performance and make fitGumbelGauss available via parameters.
      @todo allow charge state based fitting
      @todo allow semi-supervised by using decoy annotations
      @todo allow non-parametric via kernel density estimation

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
       * @brief extract and transform score types to a range and score orientation that the PEP model can handle
       * @param protein_ids the protein identifications
       * @param peptide_ids the peptide identifications
       * @param split_charge whether different charge states should be treated separately
       * @param top_hits_only only consider rank 1
       * @param target_decoy_available whether target decoy information is stored as meta value
       * @param fdr_for_targets_smaller fdr threshold for targets
       * @return engine (and optional charge state) id -> vector of triplets (score, target, decoy)
       * @note supported engines are: XTandem,OMSSA,MASCOT,SpectraST,MyriMatch,SimTandem,MSGFPlus,MS-GF+,Comet
       */
      static std::map<String, std::vector<std::vector<double>>> extractAndTransformScores(
        const std::vector<ProteinIdentification> & protein_ids,
        const std::vector<PeptideIdentification> & peptide_ids,
        const bool split_charge,
        const bool top_hits_only,
        const bool target_decoy_available,
        const double fdr_for_targets_smaller);

      /**
       * @brief update score entries with PEP (or 1-PEP) estimates
       * @param PEP_model the PEP model used to update the scores
       * @param search_engine the score of search_engine will be updated
       * @param charge identifications with the given charge will be updated
       * @param prob_correct report 1-PEP
       * @param split_charge if charge states have been treated separately
       * @param protein_ids the protein identifications
       * @param peptide_ids the peptide identifications
       * @param unable_to_fit_data there was a problem fitting the data (probabilities are all smaller 0 or larger 1)
       * @param data_might_not_be_well_fit fit was successful but of bad quality (probabilities are all smaller 0.8 and larger 0.2)
       * @note supported engines are: XTandem,OMSSA,MASCOT,SpectraST,MyriMatch,SimTandem,MSGFPlus,MS-GF+,Comet
       */
      static void updateScores(
        const PosteriorErrorProbabilityModel & PEP_model,
        const String & search_engine,
        const Int charge,
        const bool prob_correct,
        const bool split_charge,
        std::vector<ProteinIdentification> & protein_ids,
        std::vector<PeptideIdentification> & peptide_ids,
        bool & unable_to_fit_data,
        bool & data_might_not_be_well_fit);

      /**
          @brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables.
          computeProbability can be used afterwards.
          Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities.
          @param search_engine_scores a vector which holds the data points
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vector is sorted from smallest to biggest value!
      */
      bool fit(std::vector<double> & search_engine_scores);

      /**
          @brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables.
          computeProbability can be used afterwards.
          Uses Gumbel+Gauss for everything. Fits Gumbel by maximizing log likelihood.
          @param search_engine_scores a vector which holds the data points
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vector is sorted from smallest to biggest value!
      */
      bool fitGumbelGauss(std::vector<double>& search_engine_scores);

      /**
          @brief fits the distributions to the data points(search_engine_scores) and writes the computed probabilities into the given vector (the second one).
          @param search_engine_scores a vector which holds the data points
          @param probabilities a vector which holds the probability for each data point after running this function. If it has some content it will be overwritten.
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vectors are sorted from smallest to biggest value!
      */
      bool fit(std::vector<double> & search_engine_scores, std::vector<double> & probabilities);

      ///Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillDensities(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillLogDensities(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///Writes the log distributions of gumbel and gauss densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillLogDensitiesGumbel(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///computes the Likelihood with a log-likelihood function.
      double computeLogLikelihood(const std::vector<double> & incorrect_density, const std::vector<double> & correct_density);
      
      /**computes the posteriors for the datapoints to belong to the incorrect distribution
       * @param incorrect_posterior resulting posteriors
       * @return the loglikelihood of the model
       */
      double computeLLAndIncorrectPosteriorsFromLogDensities(
          const std::vector<double>& incorrect_log_density,
          const std::vector<double>& correct_log_density,
          std::vector<double>& incorrect_posterior);

      /**
       * @param x_scores Scores observed "on the x-axis"
       * @param incorrect_posteriors Posteriors/responsibilities of belonging to the incorrect component
       * @return New estimate for the mean of the correct (pair.first) and incorrect (pair.second) component
       * @note only for Gaussian estimates
       */
      std::pair<double, double> pos_neg_mean_weighted_posteriors(const std::vector<double> &x_scores,
                                                                 const std::vector<double> &incorrect_posteriors);

      /**
       * @param x_scores Scores observed "on the x-axis"
       * @param incorrect_posteriors Posteriors/responsibilities of belonging to the incorrect component
       * @return New estimate for the std. deviation of the correct (pair.first) and incorrect (pair.second) component
       * @note only for Gaussian estimates
       */
      std::pair<double, double> pos_neg_sigma_weighted_posteriors(const std::vector<double> &x_scores,
                                                                 const std::vector<double> &incorrect_posteriors,
                                                                 const std::pair<double, double>& means);

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

      ///returns estimated parameters for correctly assigned sequences. Fit should be used before.
      GumbelMaxLikelihoodFitter::GumbelDistributionFitResult getIncorrectlyAssignedGumbelFitResult() const
      {
        return incorrectly_assigned_fit_gumbel_param_;
      }

      ///returns the estimated negative prior probability.
      double getNegativePrior() const
      {
        return negative_prior_;
      }

      ///computes the gumbel density at position x with parameters params.
      static double getGumbel_(double x, const GaussFitter::GaussFitResult & params)
      {
        double z = exp((params.x0 - x) / params.sigma);
        return (z * exp(-1 * z)) / params.sigma;
      }

      /**
          Returns the computed posterior error probability for a given score.
          @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
      */
      double computeProbability(double score) const;

      /// initializes the plots
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
      /// transform different score types to a range and score orientation that the model can handle (engine string is assumed in upper-case)
      static double transformScore_(const String & engine, const PeptideHit & hit);

      /// assignment operator (not implemented)
      PosteriorErrorProbabilityModel & operator=(const PosteriorErrorProbabilityModel & rhs);
      ///Copy constructor (not implemented)
      PosteriorErrorProbabilityModel(const PosteriorErrorProbabilityModel & rhs);
      ///stores parameters for incorrectly assigned sequences. If gumbel fit was used, A can be ignored. Furthermore, in this case, x0 and sigma are the local parameter alpha and scale parameter beta, respectively.
      GaussFitter::GaussFitResult incorrectly_assigned_fit_param_;
      GumbelMaxLikelihoodFitter::GumbelDistributionFitResult incorrectly_assigned_fit_gumbel_param_;
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
      ///points either to getGumbelGnuplotFormula or getGaussGnuplotFormula depending on whether one uses the gumbel or the gaussian distribution for incorrectly assigned sequences.
      const String (PosteriorErrorProbabilityModel::* getNegativeGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;
      ///points to getGumbelGnuplotFormula
      const String (PosteriorErrorProbabilityModel::* getPositiveGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;
    };
  }
}

