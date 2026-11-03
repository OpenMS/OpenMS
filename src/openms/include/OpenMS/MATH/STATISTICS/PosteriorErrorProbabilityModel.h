// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  class TextFile;
  class PeptideIdentification;
  class ProteinIdentification;
  class PeptideHit;
  class PeptideIdentificationList;
  namespace Math
  {


    /**
      @brief Estimates posterior error probabilities (PEPs) of peptide-spectrum matches from search engine scores.

      Two estimators are available:

      - <b>Unsupervised mixture model</b> (fit(), and fitTargetDecoy() without decoys): a two-component mixture of a Gumbel
        (or Gauss) distribution for incorrect and a Gauss distribution for correct assignments is fitted to all scores with the
        EM algorithm. The very same densities are used for fitting and for computing probabilities. The Gumbel component is
        fitted by weighted maximum likelihood.
      - <b>Non-parametric target-decoy estimate</b> (fitTargetDecoy() with at least one decoy):
        Under target-decoy competition (concatenated target+decoy search, top hit per spectrum), each decoy winner stands
        in for (about) one incorrect target winner. The probability p(s) that a PSM with score s is a decoy is estimated by
        weighted isotonic regression (pool-adjacent-violators, tied scores pooled) and the PEP of a target PSM is
        <tt>PEP(s) = p(s) / (1 - p(s)) / decoy_ratio</tt> (cf. Percolator/qvality, Käll et al. 2008).
        The blocks and their PEPs only depend on the ordering of the scores, so they are invariant to monotone score
        transformations (e-values, p-values, raw scores, discrete or censored scores). Block decoy fractions are shrunk with a
        Jeffreys prior, and the PEPs of the targets scoring above all decoys sum up to 1 / decoy_ratio (the expected number of
        incorrect targets above the best decoy), so that no PEP is exactly 0. Even a single decoy (or a single target) is informative.
        Within a block, the PEP decreases linearly with the score around the block's PEP (keeping the mean PEP of the block),
        such that PSMs are ordered by score; between blocks, the PEP is interpolated linearly. These ramps use the score values,
        i.e. they (but not the block PEPs) depend on the score scale.
        As in qvality, Percolator and mokapot, the PEP is not forced to 1 at the lowest scores, where it follows the local ratio
        of decoys and targets (see hasCalibrationWarning()).

      @note Without decoys, the (unsupervised) mixture model assumes that the data contain correct assignments.
      For data with few or no correct assignments, it might assign low PEPs to a few PSMs in the extreme upper tail of the scores.
      Use target/decoy information whenever possible.

      One can output the fit as a gnuplot formula using getGumbelGnuplotFormula() and getGaussGnuplotFormula() after fitting.
      @note All parameters are stored in GaussFitResult. In the case of the Gumbel distribution x0 and sigma represent the local parameter alpha and the scale parameter beta, respectively.

      @todo allow non-parametric via kernel density estimation

      @htmlinclude OpenMS_PosteriorErrorProbabilityModel.parameters

      @ingroup Math
    */
    class OPENMS_DLLAPI PosteriorErrorProbabilityModel :
      public DefaultParamHandler
    {
public:

      /// The kind of model that was estimated by the last call to one of the fit functions
      enum class ModelType
      {
        NONE,                    ///< nothing was fitted (yet), or the last fit failed
        MIXTURE,                 ///< unsupervised two-component mixture model
        TARGET_DECOY_ISOTONIC    ///< non-parametric target-decoy estimate (isotonic regression of the decoy probability)
      };

      /// human-readable name of a ModelType
      static std::string modelTypeToString(ModelType type);

      ///default constructor
      PosteriorErrorProbabilityModel();

      ///Destructor
      ~PosteriorErrorProbabilityModel() override;

      /**
       * @brief extract and transform score types to a range and score orientation that the PEP model can handle
       *
       * Scores are transformed such that higher values are better (e.g. -log10 of e-values).
       *
       * @param[in] protein_ids the protein identifications
       * @param[in] peptide_ids the peptide identifications
       * @param[in] split_charge whether different charge states should be treated separately
       * @param[in] top_hits_only only consider rank 1 (affects the first vector of the result only)
       * @return engine (and optional charge state) id -> vector of three score vectors:
       *         [0] scores of all considered hits (decoys included, no labels used),
       *         [1] scores of the top hits of all spectra that are annotated as target ("target" or "target+decoy"),
       *         [2] scores of the top hits of all spectra that are annotated as decoy.
       *         [1] and [2] are empty if no hit is annotated as decoy.
       *         Groups with very few scores are included as well (unlike in the deprecated overload), although they might be too small to fit a model.
       * @note Supported engines (i.e. with a dedicated score selection) are: XTandem, OMSSA, Mascot, SpectraST, MyriMatch, SimTandem,
       *       MSGFPlus/MS-GF+, Comet, MSFragger, Sage, SimpleSearchEngine and ProSE. For other engines, the main score
       *       of the identifications is used, if it is not on a probability scale (e.g. q-values or PEPs).
       */
      static std::map<std::string, std::vector<std::vector<double>>> extractAndTransformScores(
        const std::vector<ProteinIdentification> & protein_ids,
        const PeptideIdentificationList & peptide_ids,
        const bool split_charge,
        const bool top_hits_only);

      /**
       * @brief extract and transform score types to a range and score orientation that the PEP model can handle
       *
       * @deprecated Use the overload without @p target_decoy_available and @p fdr_for_targets_smaller.
       *             Target and decoy scores are now determined from the target/decoy annotation of the hits,
       *             @p fdr_for_targets_smaller is ignored.
       *
       * @param[in] protein_ids the protein identifications
       * @param[in] peptide_ids the peptide identifications
       * @param[in] split_charge whether different charge states should be treated separately
       * @param[in] top_hits_only only consider rank 1
       * @param[in] target_decoy_available if false, the target and decoy score vectors are left empty
       * @param[in] fdr_for_targets_smaller ignored
       * @return engine (and optional charge state) id -> vector of triplets (score, target, decoy); only groups with more than two scores
       */
      static std::map<std::string, std::vector<std::vector<double>>> extractAndTransformScores(
        const std::vector<ProteinIdentification> & protein_ids,
        const PeptideIdentificationList & peptide_ids,
        const bool split_charge,
        const bool top_hits_only,
        const bool target_decoy_available,
        const double fdr_for_targets_smaller);

      /**
       * @brief update score entries with PEP (or 1-PEP) estimates
       * @param[in] PEP_model the PEP model used to update the scores
       * @param[in] search_engine the score of search_engine will be updated (key as returned by extractAndTransformScores())
       * @param[in] charge identifications with the given charge will be updated
       * @param[in] prob_correct report 1-PEP
       * @param[in] split_charge if charge states have been treated separately
       * @param[in,out] protein_ids the protein identifications
       * @param[in,out] peptide_ids the peptide identifications
       * @param[out] unable_to_fit_data there was a problem fitting the data (probabilities are all smaller 0 or larger 1)
       * @param[out] data_might_not_be_well_fit fit was successful but of bad quality (probabilities are all smaller 0.8 and larger 0.2)
       */
      static void updateScores(
        const PosteriorErrorProbabilityModel & PEP_model,
        const std::string & search_engine,
        const Int charge,
        const bool prob_correct,
        const bool split_charge,
        std::vector<ProteinIdentification> & protein_ids,
        PeptideIdentificationList & peptide_ids,
        bool & unable_to_fit_data,
        bool & data_might_not_be_well_fit);

      /**
       * @brief update score entries with PEP (or 1-PEP) estimates and report which hits were updated
       *
       * Same as the overload without @p updated_hits, except for hits without a usable score: these are left unchanged
       * (and not added to @p updated_hits), so that the caller can treat them consistently with the other hits of their run
       * (the overload without @p updated_hits assigns them a PEP of 1).
       * @param[in] PEP_model the PEP model used to update the scores
       * @param[in] search_engine the score of search_engine will be updated (key as returned by extractAndTransformScores())
       * @param[in] charge identifications with the given charge will be updated
       * @param[in] prob_correct report 1-PEP
       * @param[in] split_charge if charge states have been treated separately
       * @param[in,out] protein_ids the protein identifications
       * @param[in,out] peptide_ids the peptide identifications
       * @param[out] unable_to_fit_data there was a problem fitting the data (probabilities are all smaller 0 or larger 1)
       * @param[out] data_might_not_be_well_fit fit was successful but of bad quality (probabilities are all smaller 0.8 and larger 0.2)
       * @param[in,out] updated_hits pairs (index into @p peptide_ids, index of the hit) of the updated hits are added
       */
      static void updateScores(
        const PosteriorErrorProbabilityModel & PEP_model,
        const std::string & search_engine,
        const Int charge,
        const bool prob_correct,
        const bool split_charge,
        std::vector<ProteinIdentification> & protein_ids,
        PeptideIdentificationList & peptide_ids,
        bool & unable_to_fit_data,
        bool & data_might_not_be_well_fit,
        std::set<std::pair<Size, Size>> & updated_hits);

      /**
       * @brief the key under which the scores of a search run are reported by extractAndTransformScores()
       *
       * This is the upper-case name of the search engine without non-alphanumeric characters, with aliases resolved
       * (e.g. "X! Tandem" -> "XTANDEM", "MS-GF+" -> "MSGFPLUS").
       */
      static std::string getSearchEngineKey(const ProteinIdentification & protein_id);

      /**
          @brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables.
          computeProbability can be used afterwards.

          Fits a Gumbel (parameter @p incorrectly_assigned = Gumbel) or a Gauss distribution for the incorrectly assigned and a Gauss distribution
          for the correctly assigned sequences with the EM algorithm. The same distributions are used to calculate the final probabilities.

          @param[in,out] search_engine_scores a vector which holds the data points
          @param[in] outlier_handling Valid values are in the Param 'outlier_handling'
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vector is sorted from smallest to biggest value!
          @note non-finite scores are ignored for fitting
      */
      bool fit(std::vector<double>& search_engine_scores, const std::string& outlier_handling);

      /**
          @brief fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables.
          computeProbability can be used afterwards.
          Uses Gumbel+Gauss (independent of the parameter @p incorrectly_assigned). Fits the Gumbel distribution by maximizing the (weighted) likelihood.
          @param[in,out] search_engine_scores a vector which holds the data points
          @param[in] outlier_handling Valid values are in the Param 'outlier_handling'
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vector is sorted from smallest to biggest value!
      */
      bool fitGumbelGauss(std::vector<double>& search_engine_scores, const std::string& outlier_handling);

      /**
          @brief fits the distributions to the data points(search_engine_scores) and writes the computed probabilities into the given vector (the second one).
          @param[in,out] search_engine_scores a vector which holds the data points
          @param[out] probabilities Probability for each data point after running this function. If it has some content it will be overwritten.
          @param[in] outlier_handling Valid values are in the Param 'outlier_handling'
          @return true if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated.
          @note the vectors are sorted from smallest to biggest value!
      */
      bool fit(std::vector<double>& search_engine_scores, std::vector<double>& probabilities, const std::string& outlier_handling);

      /**
          @brief estimates the PEP model from scores of target and decoy PSMs. computeProbability can be used afterwards.

          The scores must be oriented such that higher scores are better (as returned by extractAndTransformScores()).
          They should stem from a target-decoy competition, i.e. a search against a concatenated target+decoy database,
          with one (the top) hit per spectrum. With at least one decoy, the non-parametric target-decoy estimate is used;
          without decoys (or with @p decoy_handling = ignore), the unsupervised mixture model (see class documentation).
          Use getModelType() to find out which estimator was used.

          @param[in] target_scores scores of target PSMs
          @param[in] decoy_scores scores of decoy PSMs
          @param[in] outlier_handling Valid values are in the Param 'outlier_handling' (only used for the mixture model)
          @return true if a model could be estimated
      */
      bool fitTargetDecoy(const std::vector<double>& target_scores, const std::vector<double>& decoy_scores, const std::string& outlier_handling);

      /// returns the kind of model estimated by the last fit
      ModelType getModelType() const
      {
        return model_type_;
      }

      /**
          @brief returns the target-decoy estimate of the fraction of incorrect target PSMs, i.e.
          (number of decoys) / (number of targets) / decoy_ratio (capped at 1), as computed by the last call to fitTargetDecoy().
          Returns -1 if it is not available.
      */
      double getTargetDecoyIncorrectFraction() const
      {
        return td_incorrect_fraction_;
      }

      /**
          @brief returns true if the last call to fitTargetDecoy() found the target-decoy assumptions to be violated at low scores

          A warning is raised if even the lowest-scoring PSMs contain significantly (binomial test) and clearly too few decoys
          for incorrect PSMs (decoy fraction below decoy_ratio / (1 + decoy_ratio), PEP below 0.8), so that the PEPs stay well below 1
          down to the lowest scores. Then either correct PSMs extend to the lowest scores, or decoys are underrepresented among the
          incorrect PSMs (e.g. a wrong decoy_ratio), which makes all target-decoy PEPs (and FDRs) too low.
          The set of lowest-scoring PSMs with the highest decoy fraction among those with at least a quarter of the estimated number of
          incorrect PSMs (and at least 100 PSMs) is tested.
      */
      bool hasCalibrationWarning() const
      {
        return calibration_warning_;
      }

      /**
          @brief returns the PEP curve of the non-parametric target-decoy estimate as pairs (score, PEP), sorted by score.

          Probabilities of scores in between are interpolated linearly, probabilities of scores outside are constant.
          Empty unless the last fit was of type ModelType::TARGET_DECOY_ISOTONIC.
      */
      std::vector<std::pair<double, double>> getIsotonicPEPCurve() const;

      ///Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillDensities(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillLogDensities(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///Writes the log distributions of gumbel and gauss densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
      void fillLogDensitiesGumbel(const std::vector<double> & x_scores, std::vector<double> & incorrect_density, std::vector<double> & correct_density);
      ///computes the Likelihood with a log-likelihood function.
      double computeLogLikelihood(const std::vector<double> & incorrect_density, const std::vector<double> & correct_density) const;

      /**computes the posteriors for the datapoints to belong to the incorrect distribution
        @param[in] incorrect_log_density ...
        @param[in] correct_log_density ...
        @param[out] incorrect_posterior resulting posteriors; same size as @p incorrect_log_density
        @return the log-likelihood of the model
      */
      double computeLLAndIncorrectPosteriorsFromLogDensities(
          const std::vector<double>& incorrect_log_density,
          const std::vector<double>& correct_log_density,
          std::vector<double>& incorrect_posterior) const;

      /**
       * @param[in] x_scores Scores observed "on the x-axis"
       * @param[in] incorrect_posteriors Posteriors/responsibilities of belonging to the incorrect component
       * @return New estimate for the mean of the correct (pair.first) and incorrect (pair.second) component
       * @note only for Gaussian estimates
       */
      std::pair<double, double> pos_neg_mean_weighted_posteriors(const std::vector<double> &x_scores,
                                                                 const std::vector<double> &incorrect_posteriors);

      /**
       * @param[in] x_scores Scores observed "on the x-axis"
       * @param[in] incorrect_posteriors Posteriors/responsibilities of belonging to the incorrect component
       * @param[in] pos_neg_mean Positive(correct) and negative(incorrect) means, respectively
       * @return New estimate for the std. deviation of the correct (pair.first) and incorrect (pair.second) component
       * @note only for Gaussian estimates
       */
      std::pair<double, double> pos_neg_sigma_weighted_posteriors(const std::vector<double> &x_scores,
                                                                 const std::vector<double> &incorrect_posteriors,
                                                                 const std::pair<double, double>& pos_neg_mean);

      ///returns estimated parameters for correctly assigned sequences. Fit should be used before.
      GaussFitter::GaussFitResult getCorrectlyAssignedFitResult() const
      {
        return correctly_assigned_fit_param_;
      }

      ///returns estimated parameters for incorrectly assigned sequences (for the Gumbel distribution, x0 and sigma are the location and scale parameter). Fit should be used before.
      GaussFitter::GaussFitResult getIncorrectlyAssignedFitResult() const
      {
        return incorrectly_assigned_fit_param_;
      }

      ///returns estimated Gumbel parameters for incorrectly assigned sequences. Fit should be used before.
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

          The probability is non-increasing with the score.
          @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
      */
      double computeProbability(double score) const;

      /// initializes the plots
      TextFile initPlots(std::vector<double> & x_scores);

      /// returns the gnuplot formula of the fitted gumbel distribution. Only x0 and sigma are used as local parameter alpha and scale parameter beta, respectively.
      const std::string getGumbelGnuplotFormula(const GaussFitter::GaussFitResult & params) const;

      /// returns the gnuplot formula of the fitted gauss distribution.
      const std::string getGaussGnuplotFormula(const GaussFitter::GaussFitResult & params) const;

      /// returns the gnuplot formula of the fitted mixture distribution.
      const std::string getBothGnuplotFormula(const GaussFitter::GaussFitResult & incorrect, const GaussFitter::GaussFitResult & correct) const;

      /**
          @brief plots the score distributions of target and decoy PSMs together with the estimated model

          The decoy histogram is scaled such that, under target-decoy competition, it estimates the score distribution of the incorrect target PSMs.
          @param[in,out] target scores of target PSMs (will be sorted)
          @param[in,out] decoy scores of decoy PSMs (will be sorted)
      */
      void plotTargetDecoyEstimation(std::vector<double> & target, std::vector<double> & decoy);

      /// returns the smallest score used in the last fit
      inline double getSmallestScore() const
      {
        return smallest_score_;
      }

      /// try to invoke 'gnuplot' on the file to create PDF automatically
      void tryGnuplot(const std::string& gp_file);

private:
      /// resets the state of the last fit (model type, target-decoy diagnostics, non-parametric curve)
      void resetFit_();

      /// remove (or clamp) outliers for fitting
      void processOutliers_(std::vector<double>& x_scores, const std::string& outlier_handling) const;

      /**
          @brief fits a two-component mixture model with the EM algorithm

          @param[in] scores (finite) scores of the PSMs
          @param[in] use_gumbel use a Gumbel (instead of a Gauss) distribution for the incorrect component
          @param[in] outlier_handling see parameter 'outlier_handling'
          @return true on success
      */
      bool fitMixture_(std::vector<double> scores, bool use_gumbel, const std::string& outlier_handling);

      /// fits the target-decoy model (or, without decoys or with decoy_handling 'ignore', the mixture model) to finite target and decoy scores
      bool fitTargetDecoyModel_(const std::vector<double>& targets, const std::vector<double>& decoys, const std::string& outlier_handling);

      /// estimates the non-parametric target-decoy model; returns false if it cannot be estimated
      bool fitIsotonic_(const std::vector<double>& target_scores, const std::vector<double>& decoy_scores);

      /// log density of the incorrect component at @p x (in the shifted score space)
      double logIncorrectDensity_(double x) const;

      /// log density of the correct component at @p x (in the shifted score space)
      double logCorrectDensity_(double x) const;

      /// mode of the incorrect component (in the shifted score space)
      double incorrectMode_() const;

      /// transform different score types to a range and score orientation that the model can handle (engine string is assumed in upper-case)
      /// @param[in] engine the search engine key as returned by getSearchEngineKey()
      /// @param[in] hit the PeptideHit to extract transformed scores from
      /// @param[in] current_score_type the current score type of the PeptideIdentification to take precedence
      /// @param[in] higher_score_better the score orientation of the PeptideIdentification (used for the main score of unknown engines)
      /// @param[in] score_source the score type to use (see chooseScoreSources_()); empty: the main score of the identification
      /// @return transformed score (higher is better), or NaN if no usable score is available
      static double transformScore_(const std::string& engine, const PeptideHit& hit, const std::string& current_score_type, bool higher_score_better, const std::string& score_source);

      /// For each engine (key) of @p run_to_engine (run identifier -> engine key): the score type to use for all of its hits.
      /// For engines with dedicated support, the most preferred of their score types that any hit has; otherwise (or if no hit
      /// has any of them) an empty string, i.e. the main score of the identifications. One score type per engine keeps the scores
      /// of a model on the same scale.
      /// implementation of updateScores(); hits without a usable score get a PEP of 1 if @p unusable_to_one is set, otherwise they are left unchanged
      static void updateScores_(
        const PosteriorErrorProbabilityModel & PEP_model,
        const std::string & search_engine,
        const Int charge,
        const bool prob_correct,
        const bool split_charge,
        std::vector<ProteinIdentification> & protein_ids,
        PeptideIdentificationList & peptide_ids,
        bool & unable_to_fit_data,
        bool & data_might_not_be_well_fit,
        bool unusable_to_one,
        std::set<std::pair<Size, Size>> & updated_hits);

      static std::map<std::string, std::string> chooseScoreSources_(const std::map<std::string, std::string>& run_to_engine, const PeptideIdentificationList& peptide_ids);

      /// gets a specific score (either main score [preferred] or metavalue)
      /// @param[in] requested_score_types the requested score_types in order of preference (will be tested with a "_score" suffix as well)
      /// @param[in] hit the PeptideHit to extract from
      /// @param[in] actual_score_type the current score type to take preference if matching
      /// @param[out] value the score
      /// @return true if one of the requested scores was found
      static bool getScore_(const std::vector<std::string>& requested_score_types, const PeptideHit & hit, const std::string& actual_score_type, double& value);

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
      ///smallest score which was used for fitting the model
      double smallest_score_;
      ///points either to getGumbelGnuplotFormula or getGaussGnuplotFormula depending on whether one uses the gumbel or the gaussian distribution for incorrectly assigned sequences.
      const std::string (PosteriorErrorProbabilityModel::* getNegativeGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;
      ///points to getGumbelGnuplotFormula
      const std::string (PosteriorErrorProbabilityModel::* getPositiveGnuplotFormula_)(const GaussFitter::GaussFitResult & params) const;
      /// kind of the last estimated model
      ModelType model_type_;
      /// whether the incorrect component of the (last) mixture model is a Gumbel distribution (otherwise Gauss)
      bool incorrect_is_gumbel_;
      /// knots of the non-parametric PEP curve: scores (ascending) and PEPs (non-increasing)
      std::vector<double> isotonic_scores_;
      std::vector<double> isotonic_peps_;
      /// target-decoy estimate of the fraction of incorrect targets (-1 if not available)
      double td_incorrect_fraction_;
      /// whether the last target-decoy fit found too few decoys among the lowest-scoring PSMs (see hasCalibrationWarning())
      bool calibration_warning_;
    };
  }
}
