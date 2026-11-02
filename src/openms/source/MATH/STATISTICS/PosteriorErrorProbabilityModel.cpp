// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>

#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <filesystem>

#include <boost/math/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>



using namespace std;

namespace OpenMS::Math
{
    namespace
    {
      /// Euler-Mascheroni constant (mean of the standard Gumbel distribution)
      constexpr double EULER_GAMMA = 0.57721566490153286;

      double logGaussDensity(double x, double mu, double sigma)
      {
        const double z = (x - mu) / sigma;
        return -std::log(sigma) - 0.5 * std::log(2.0 * Constants::PI) - 0.5 * z * z;
      }

      /// weighted mean and (biased) standard deviation
      std::pair<double, double> weightedMeanSD(const std::vector<double>& x, const std::vector<double>& w)
      {
        double sum_w = 0.0, sum_wx = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
        {
          sum_w += w[i];
          sum_wx += w[i] * x[i];
        }
        if (!(sum_w > 0.0))
        {
          return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }
        const double mean = sum_wx / sum_w;
        double sum_wvar = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
        {
          sum_wvar += w[i] * (x[i] - mean) * (x[i] - mean);
        }
        return {mean, std::sqrt(sum_wvar / sum_w)};
      }

      /// mean and (biased) standard deviation of a range
      std::pair<double, double> meanSD(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end)
      {
        const std::vector<double> x(begin, end);
        return weightedMeanSD(x, std::vector<double>(x.size(), 1.0));
      }

      enum class TDLabel
      {
        TARGET,
        DECOY,
        UNKNOWN
      };

      /// target/decoy status of a hit, parsed defensively (missing, empty or unexpected annotations yield UNKNOWN)
      TDLabel getTDLabel(const PeptideHit& hit)
      {
        if (!hit.metaValueExists(Constants::UserParam::TARGET_DECOY))
        {
          return TDLabel::UNKNOWN;
        }
        std::string value = hit.getMetaValue(Constants::UserParam::TARGET_DECOY).toString();
        StringUtils::trim(value);
        StringUtils::toLower(value);
        if (value == "target" || value == "target+decoy")
        {
          return TDLabel::TARGET;
        }
        if (value == "decoy")
        {
          return TDLabel::DECOY;
        }
        return TDLabel::UNKNOWN;
      }

      /// upper-case and remove all non-alphanumeric characters, e.g. "X! Tandem" -> "XTANDEM"
      std::string compactEngineName(const std::string& name)
      {
        std::string result;
        for (char c : name)
        {
          if (std::isalnum(static_cast<unsigned char>(c)))
          {
            result += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
          }
        }
        return result;
      }

      /// resolve aliases of search engine names
      std::string normalizeEngineKey(const std::string& name)
      {
        std::string key = compactEngineName(name);
        if (key == "MSGF") // "MS-GF+"
        {
          key = "MSGFPLUS";
        }
        if (key.empty())
        {
          key = "UNKNOWN";
        }
        return key;
      }

      enum class ScoreScale
      {
        RAW,         ///< raw score (orientation given by the identification)
        EVALUE,      ///< e-value or p-value: transformed with -log10
        PROBABILITY  ///< probability-like (PEP, posterior probability, q-value, FDR): not usable as input
      };

      /// classify the main score of an identification of an engine without dedicated support
      ScoreScale classifyScore(const std::string& score_type)
      {
        Scores::IDType type;
        if (Scores::findIDTypeByName(score_type, type))
        {
          switch (type)
          {
            case Scores::IDType::RAW_EVAL:
              return ScoreScale::EVALUE;
            case Scores::IDType::PP:
            case Scores::IDType::PEP:
            case Scores::IDType::FDR:
            case Scores::IDType::QVAL:
              return ScoreScale::PROBABILITY;
            default:
              break;
          }
        }
        std::string name = Scores::normalizeScoreName(score_type);
        StringUtils::toLower(name);
        for (const char* probability_name : {"posterior", "q-value", "qvalue", "q_value", "fdr", "probability"})
        {
          if (StringUtils::hasSubstring(name, probability_name))
          {
            return ScoreScale::PROBABILITY;
          }
        }
        // short abbreviations only as separate words (e.g. "Andromeda:PEP", but not "PeptideScore")
        std::string token;
        for (char c : name + " ")
        {
          if (std::isalnum(static_cast<unsigned char>(c)))
          {
            token += c;
            continue;
          }
          if (token == "pep" || token == "pp")
          {
            return ScoreScale::PROBABILITY;
          }
          token.clear();
        }
        for (const char* evalue_name : {"evalue", "e-value", "e_value", "expect", "p-value", "pvalue", "p_value"})
        {
          if (StringUtils::hasSubstring(name, evalue_name))
          {
            return ScoreScale::EVALUE;
          }
        }
        return ScoreScale::RAW;
      }

      /// strict weak ordering that sorts NaN values to the end
      bool nanLastLess(double a, double b)
      {
        return std::isnan(b) ? !std::isnan(a) : a < b;
      }

      /// -log10 of an e-value or p-value; 0 is mapped to the smallest positive double, negative or NaN values yield NaN
      double negLog10(double value)
      {
        if (!(value >= 0.0) || std::isinf(value))
        {
          return std::numeric_limits<double>::quiet_NaN();
        }
        return -std::log10(std::max(value, std::numeric_limits<double>::min()));
      }

      /// block of consecutive PSMs (sorted by score) for the pool-adjacent-violators algorithm
      struct PAVABlock
      {
        double n;          ///< number of PSMs
        double decoys;     ///< number of decoys
        double min_score;  ///< smallest score in the block
        double max_score;  ///< largest score in the block
        double score_sum;  ///< sum of the scores of the PSMs in the block

        /// fraction of decoys, optionally with a Jeffreys prior (pseudo-count of 1/2 decoy and 1/2 target)
        double decoyFraction(bool jeffreys) const
        {
          return jeffreys ? (decoys + 0.5) / (n + 1.0) : decoys / n;
        }

        void merge(const PAVABlock& other)
        {
          n += other.n;
          decoys += other.decoys;
          min_score = std::min(min_score, other.min_score);
          max_score = std::max(max_score, other.max_score);
          score_sum += other.score_sum;
        }
      };

      /**
        @brief adds @p block to @p blocks (sorted from best to worst score) and merges blocks until their decoy fractions are non-decreasing

        @param merge_equal also merge neighbours with equal decoy fractions (yields maximal blocks; does not change the fitted values)
        @param jeffreys use the decoy fractions with a Jeffreys prior (one pseudo-count per block)
      */
      void addNonDecreasing(std::vector<PAVABlock>& blocks, PAVABlock block, bool merge_equal, bool jeffreys)
      {
        blocks.push_back(block);
        while (blocks.size() > 1)
        {
          const double previous = blocks[blocks.size() - 2].decoyFraction(jeffreys);
          const double last = blocks.back().decoyFraction(jeffreys);
          if (previous < last || (previous == last && !merge_equal))
          {
            break;
          }
          PAVABlock last_block = blocks.back();
          blocks.pop_back();
          blocks.back().merge(last_block);
        }
      }

      /// area under the ROC curve for separating targets (positives) from decoys by score (ties count half)
      double targetDecoyAUC(const std::vector<double>& targets, const std::vector<double>& decoys)
      {
        std::vector<std::pair<double, bool>> all; // (score, is_target)
        all.reserve(targets.size() + decoys.size());
        for (double t : targets) { all.emplace_back(t, true); }
        for (double d : decoys) { all.emplace_back(d, false); }
        std::sort(all.begin(), all.end(), [](const std::pair<double, bool>& a, const std::pair<double, bool>& b) { return a.first < b.first; });
        // Mann-Whitney U statistic with average ranks for ties
        double rank_sum_targets = 0.0;
        for (Size i = 0; i < all.size(); )
        {
          Size j = i;
          Size n_targets_tied = 0;
          for (; j < all.size() && all[j].first == all[i].first; ++j)
          {
            n_targets_tied += all[j].second ? 1 : 0;
          }
          const double average_rank = 0.5 * double(i + 1 + j); // ranks i+1 .. j
          rank_sum_targets += average_rank * double(n_targets_tied);
          i = j;
        }
        const double n_t = double(targets.size()), n_d = double(decoys.size());
        return (rank_sum_targets - n_t * (n_t + 1.0) / 2.0) / (n_t * n_d);
      }
    }

    std::string PosteriorErrorProbabilityModel::modelTypeToString(ModelType type)
    {
      switch (type)
      {
        case ModelType::NONE: return "none";
        case ModelType::MIXTURE: return "unsupervised mixture model";
        case ModelType::SEMI_SUPERVISED_MIXTURE: return "semi-supervised (decoy-anchored) mixture model";
        case ModelType::TARGET_DECOY_ISOTONIC: return "non-parametric target-decoy estimate (isotonic regression)";
      }
      return "unknown";
    }

    PosteriorErrorProbabilityModel::PosteriorErrorProbabilityModel() :
      DefaultParamHandler("PosteriorErrorProbabilityModel"),
      incorrectly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      incorrectly_assigned_fit_gumbel_param_(GumbelMaxLikelihoodFitter::GumbelDistributionFitResult(-1,-1)),
      correctly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      negative_prior_(0.5), max_incorrectly_(0), max_correctly_(0), smallest_score_(0),
      model_type_(ModelType::NONE), incorrect_is_gumbel_(true), td_incorrect_fraction_(-1), calibration_warning_(false)
    {
      defaults_.setValue("out_plot", "", "If given, the some output files will be saved in the following manner: <out_plot>_scores.txt for the scores and <out_plot> which contains the fitted values for each step of the EM-algorithm, e.g., out_plot = /usr/home/OMSSA123 leads to /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123 will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123' just OMSSA123, the files will be written into the working directory.", {"advanced","output file"});
      defaults_.setValue("number_of_bins", 100, "Number of bins used for visualization. Only needed if each iteration step of the EM-Algorithm will be visualized", {"advanced"});
      defaults_.setValue("incorrectly_assigned", "Gumbel", "Distribution of the scores of incorrectly assigned sequences in the mixture models. For 'Gumbel', the Gumbel distribution is used (fitted by maximum likelihood). For 'Gauss', the Gauss distribution is used.", {"advanced"});
      defaults_.setValue("max_nr_iterations", 1000, "Bounds the number of iterations for the EM algorithm when convergence is slow.", {"advanced"});
      defaults_.setValidStrings("incorrectly_assigned", {"Gumbel","Gauss"});
      defaults_.setValue("neg_log_delta",6, "The negative logarithm of the convergence threshold for the likelihood increase.");
      defaults_.setValue("outlier_handling","ignore_iqr_outliers", "What to do with outliers (mixture models only):\n"
                                                                   "- ignore_iqr_outliers: ignore outliers outside of 3*IQR from Q1/Q3 for fitting\n"
                                                                   "- set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting\n"
                                                                   "- ignore_extreme_percentiles: ignore everything outside 99.9th and 1st percentile (also removes equal values like potential censored max values in XTandem)\n"
                                                                   "- none: do nothing");
      defaults_.setValidStrings("outlier_handling", {"ignore_iqr_outliers","set_iqr_to_closest_valid","ignore_extreme_percentiles","none"});
      defaults_.setValue("decoy_handling", "auto", "How target/decoy annotations are used (if decoy PSMs are present):\n"
                                                   "- auto: 'isotonic' if there are at least 'min_decoys_isotonic' decoys, 'semi_supervised' if there are at least 'min_decoys' decoys, otherwise an unsupervised mixture model is fitted to the target scores\n"
                                                   "- isotonic: non-parametric estimate from the local ratio of decoys and targets (isotonic regression). Invariant to monotone transformations of the score, i.e. suitable for any score type.\n"
                                                   "- semi_supervised: mixture model in which decoys belong to the incorrect component\n"
                                                   "- ignore: ignore the annotations, i.e. targets and decoys are pooled and modeled without supervision (behavior of OpenMS 3.5 and earlier)\n"
                                                   "Decoy-based estimates assume a search against a concatenated target+decoy database (target-decoy competition). Only top hits are used for fitting.");
      defaults_.setValidStrings("decoy_handling", {"auto", "isotonic", "semi_supervised", "ignore"});
      defaults_.setValue("decoy_ratio", 1.0, "Number of decoy sequences per target sequence in the searched database (cf. DecoyDatabase's 'shuffle_decoy_ratio'). Under target-decoy competition, the number of incorrect targets is estimated as the number of decoys divided by this ratio.", {"advanced"});
      defaults_.setMinFloat("decoy_ratio", 1e-6);
      defaults_.setValue("min_decoys", 10, "Minimum number of decoy PSMs (per model, i.e. per charge state with 'split_charge') to use decoy information. Otherwise an unsupervised mixture model is fitted to the target scores.", {"advanced"});
      defaults_.setMinInt("min_decoys", 1);
      defaults_.setValue("min_decoys_isotonic", 50, "Minimum number of decoy PSMs (per model) for the 'isotonic' estimate with decoy_handling 'auto'.", {"advanced"});
      defaults_.setMinInt("min_decoys_isotonic", 1);
      defaultsToParam_();
      getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
    }

    PosteriorErrorProbabilityModel::~PosteriorErrorProbabilityModel() = default;

    double PosteriorErrorProbabilityModel::logIncorrectDensity_(double x) const
    {
      if (incorrect_is_gumbel_)
      {
        return incorrectly_assigned_fit_gumbel_param_.log_eval_no_normalize(x);
      }
      return logGaussDensity(x, incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_.sigma);
    }

    double PosteriorErrorProbabilityModel::logCorrectDensity_(double x) const
    {
      return logGaussDensity(x, correctly_assigned_fit_param_.x0, correctly_assigned_fit_param_.sigma);
    }

    double PosteriorErrorProbabilityModel::incorrectMode_() const
    {
      return incorrect_is_gumbel_ ? incorrectly_assigned_fit_gumbel_param_.mode() : incorrectly_assigned_fit_param_.x0;
    }

    void PosteriorErrorProbabilityModel::resetFit_()
    {
      model_type_ = ModelType::NONE;
      td_incorrect_fraction_ = -1;
      calibration_warning_ = false;
      isotonic_scores_.clear();
      isotonic_peps_.clear();
    }

    bool PosteriorErrorProbabilityModel::fitGumbelGauss(std::vector<double>& search_engine_scores, const std::string& outlier_handling)
    {
      resetFit_();
      std::sort(search_engine_scores.begin(), search_engine_scores.end(), nanLastLess);
      std::vector<double> finite_scores;
      std::copy_if(search_engine_scores.begin(), search_engine_scores.end(), std::back_inserter(finite_scores), [](double x) { return std::isfinite(x); });
      if (finite_scores.size() < 2) { return false; }

      if (!fitMixture_(finite_scores, {}, true, outlier_handling)) { return false; }
      model_type_ = ModelType::MIXTURE;
      return true;
    }

    bool PosteriorErrorProbabilityModel::fit(std::vector<double>& search_engine_scores, const std::string& outlier_handling)
    {
      resetFit_();
      // nothing to fit?
      if (search_engine_scores.empty()) { return false; }

      std::sort(search_engine_scores.begin(), search_engine_scores.end(), nanLastLess);
      std::vector<double> finite_scores;
      std::copy_if(search_engine_scores.begin(), search_engine_scores.end(), std::back_inserter(finite_scores), [](double x) { return std::isfinite(x); });
      if (finite_scores.size() < 2) { return false; }

      const bool use_gumbel = (param_.getValue("incorrectly_assigned") == "Gumbel");
      if (!fitMixture_(finite_scores, {}, use_gumbel, outlier_handling)) { return false; }
      model_type_ = ModelType::MIXTURE;
      return true;
    }

    bool PosteriorErrorProbabilityModel::fitMixture_(std::vector<double> scores, std::vector<double> decoy_scores, bool use_gumbel, const std::string& outlier_handling)
    {
      if (scores.size() < 2) { return false; }

      //-------------------------------------------------------------
      // Initializing Parameters
      //-------------------------------------------------------------
      // transform to a positive range (the shift is irrelevant for the location-scale families, but kept for compatibility)
      smallest_score_ = *std::min_element(scores.begin(), scores.end());
      if (!decoy_scores.empty())
      {
        smallest_score_ = std::min(smallest_score_, *std::min_element(decoy_scores.begin(), decoy_scores.end()));
      }
      const double shift = fabs(smallest_score_) + 0.001;
      for (double& d : scores) { d += shift; }
      for (double& d : decoy_scores) { d += shift; }
      sort(scores.begin(), scores.end());
      sort(decoy_scores.begin(), decoy_scores.end());

      processOutliers_(scores, outlier_handling);
      if (!decoy_scores.empty()) { processOutliers_(decoy_scores, outlier_handling); }

      std::vector<double>& x = scores;
      std::vector<double>& xd = decoy_scores;
      const Size n = x.size();
      const Size nd = xd.size();
      if (n < 2) { return false; }

      incorrect_is_gumbel_ = use_gumbel;
      getNegativeGnuplotFormula_ = use_gumbel ? &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula : &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;

      const double overall_sd = meanSD(x.cbegin(), x.cend()).second;
      const double min_sd = 1e-6 * std::max(x.back() - x.front(), 1e-6);
      if (!(overall_sd > min_sd))
      {
        OPENMS_LOG_WARN << "Warning: scores have (almost) no spread. Unable to fit the mixture model." << std::endl;
        return false;
      }

      // incorrect component: from decoys if available, otherwise from the lower half of the scores
      std::pair<double, double> incorrect_ms = (nd >= 2) ? meanSD(xd.cbegin(), xd.cend()) : meanSD(x.cbegin(), x.cbegin() + std::max<Size>(2, n / 2));
      if (!(incorrect_ms.second > min_sd)) { incorrect_ms.second = overall_sd; }
      // correct component: from the upper 30% of the scores
      const Size correct_start = std::min(n - 2, static_cast<Size>(std::ceil(0.7 * n)));
      std::pair<double, double> correct_ms = meanSD(x.cbegin() + correct_start, x.cend());
      if (!(correct_ms.second > min_sd)) { correct_ms.second = overall_sd; }

      GumbelMaxLikelihoodFitter gumbel_fitter;
      if (use_gumbel)
      {
        // method of moments (Gumbel: sd = pi / sqrt(6) * b, mean = a + gamma * b)
        const double b = incorrect_ms.second * std::sqrt(6.0) / Constants::PI;
        incorrectly_assigned_fit_gumbel_param_ = GumbelMaxLikelihoodFitter::GumbelDistributionFitResult(incorrect_ms.first - EULER_GAMMA * b, b);
        if (nd >= 2)
        {
          try
          {
            gumbel_fitter.setInitialParameters(incorrectly_assigned_fit_gumbel_param_);
            incorrectly_assigned_fit_gumbel_param_ = gumbel_fitter.fitWeighted(xd, std::vector<double>(nd, 1.0));
          }
          catch (Exception::UnableToFit&)
          {
            // keep the moment estimates
          }
        }
        gumbel_fitter.setInitialParameters(incorrectly_assigned_fit_gumbel_param_);
      }
      else
      {
        incorrectly_assigned_fit_param_.x0 = incorrect_ms.first;
        incorrectly_assigned_fit_param_.sigma = incorrect_ms.second;
      }
      correctly_assigned_fit_param_.x0 = correct_ms.first;
      correctly_assigned_fit_param_.sigma = correct_ms.second;
      correctly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2.0));

      // prior probability of an incorrect assignment (of a target)
      if (nd > 0)
      {
        const double decoy_ratio = param_.getValue("decoy_ratio");
        negative_prior_ = std::clamp(double(nd) / (decoy_ratio * double(n)), 0.05, 0.95);
      }
      else
      {
        negative_prior_ = 0.7;
      }

      // keep the GaussFitResult representation in sync (used for plotting and by the getters)
      auto sync_incorrect_params = [&]()
      {
        if (use_gumbel)
        {
          incorrectly_assigned_fit_param_.x0 = incorrectly_assigned_fit_gumbel_param_.a;
          incorrectly_assigned_fit_param_.sigma = incorrectly_assigned_fit_gumbel_param_.b;
          incorrectly_assigned_fit_param_.A = std::exp(-1.0) / incorrectly_assigned_fit_gumbel_param_.b; // density at the mode
        }
        else
        {
          incorrectly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2.0));
        }
      };
      sync_incorrect_params();

      //-------------------------------------------------------------
      // create files for output
      //-------------------------------------------------------------
      bool output_plots  = (!StringUtils::trimmed(param_.getValue("out_plot").toString()).empty());
      TextFile file;
      if (output_plots)
      {
        // create output directory (if not already present)
        namespace fs = std::filesystem;
        auto plot_path = to_path(std::string(param_.getValue("out_plot").toString()));
        auto parent_dir = plot_path.parent_path();
        if (parent_dir.empty())
        {
          OPENMS_LOG_ERROR << "Could not navigate to output directory for plots from '" << plot_path.filename().string() << "'." << std::endl;
          return false;
        }
        if (!fs::exists(parent_dir) && !File::makeDir(parent_dir.generic_string()))
        {
          OPENMS_LOG_ERROR << "Could not create output directory for plots '" << parent_dir.generic_string() << "'." << std::endl;
          return false;
        }
        //
        file = initPlots(x);
      }

      //-------------------------------------------------------------
      // Estimate Parameters - EM algorithm
      //-------------------------------------------------------------
      const Int max_itns = param_.getValue("max_nr_iterations");
      const int delta = param_.getValue("neg_log_delta");

      // E-step: responsibilities of the incorrect component (targets/unlabeled PSMs only) and the log-likelihood
      std::vector<double> incorrect_posteriors(n);
      auto e_step = [&]() -> double
      {
        double loglikelihood = 0.0;
        const double log_prior_neg = std::log(negative_prior_);
        const double log_prior_pos = std::log(1.0 - negative_prior_);
        for (Size i = 0; i < n; ++i)
        {
          const double log_resp_incorrect = log_prior_neg + logIncorrectDensity_(x[i]);
          const double log_resp_correct = log_prior_pos + logCorrectDensity_(x[i]);
          const double max_log_resp = std::max(log_resp_correct, log_resp_incorrect);
          const double resp_incorrect = std::exp(log_resp_incorrect - max_log_resp);
          const double resp_correct = std::exp(log_resp_correct - max_log_resp);
          incorrect_posteriors[i] = resp_incorrect / (resp_incorrect + resp_correct);
          loglikelihood += max_log_resp + std::log(resp_incorrect + resp_correct);
        }
        for (double d : xd)
        {
          loglikelihood += logIncorrectDensity_(d); // decoys are known to be incorrect
        }
        return loglikelihood;
      };

      // all PSMs that (may) contribute to the incorrect component: targets/unlabeled PSMs, then decoys
      std::vector<double> x_incorrect(x);
      x_incorrect.insert(x_incorrect.end(), xd.begin(), xd.end());
      std::vector<double> w_incorrect(n + nd, 1.0);

      double maxlike = e_step();
      if (std::isnan(maxlike))
      {
        OPENMS_LOG_WARN << "Numerical instabilities. Aborting." << endl;
        return false;
      }
      int itns = 0;
      while (true)
      {
        //-------------------------------------------------------------
        // M-STEP
        double sum_incorrect_posteriors = std::accumulate(incorrect_posteriors.begin(), incorrect_posteriors.end(), 0.0);
        const double sum_correct_posteriors = double(n) - sum_incorrect_posteriors;
        if (!(sum_correct_posteriors > 1e-6) || !(sum_incorrect_posteriors + double(nd) > 1e-6))
        {
          OPENMS_LOG_WARN << "Warning: all scores were assigned to a single component. Unable to fit the mixture model." << std::endl;
          return false;
        }
        negative_prior_ = std::clamp(sum_incorrect_posteriors / double(n), 1e-6, 1.0 - 1e-6);

        // correct component (Gauss)
        std::vector<double> correct_weights(n);
        for (Size i = 0; i < n; ++i) { correct_weights[i] = 1.0 - incorrect_posteriors[i]; }
        const std::pair<double, double> correct_new = weightedMeanSD(x, correct_weights);

        // incorrect component (Gumbel or Gauss)
        std::copy(incorrect_posteriors.begin(), incorrect_posteriors.end(), w_incorrect.begin());
        GumbelMaxLikelihoodFitter::GumbelDistributionFitResult gumbel_new = incorrectly_assigned_fit_gumbel_param_;
        std::pair<double, double> incorrect_new(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_.sigma);
        if (use_gumbel)
        {
          try
          {
            gumbel_new = gumbel_fitter.fitWeighted(x_incorrect, w_incorrect);
          }
          catch (Exception::UnableToFit&)
          {
            OPENMS_LOG_WARN << "Warning: unable to fit the Gumbel distribution. Aborting fit." << std::endl;
            return false;
          }
        }
        else
        {
          incorrect_new = weightedMeanSD(x_incorrect, w_incorrect);
        }

        const double incorrect_scale = use_gumbel ? gumbel_new.b : incorrect_new.second;
        if (!(correct_new.second > min_sd) || !(incorrect_scale > min_sd) || std::isnan(correct_new.first) || std::isnan(incorrect_new.first))
        {
          OPENMS_LOG_WARN << "Warning: encountered impossible standard deviations. Aborting fit." << std::endl;
          return false;
        }

        // update parameters
        correctly_assigned_fit_param_.x0 = correct_new.first;
        correctly_assigned_fit_param_.sigma = correct_new.second;
        correctly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(correct_new.second, 2.0));
        if (use_gumbel)
        {
          incorrectly_assigned_fit_gumbel_param_ = gumbel_new;
        }
        else
        {
          incorrectly_assigned_fit_param_.x0 = incorrect_new.first;
          incorrectly_assigned_fit_param_.sigma = incorrect_new.second;
        }
        sync_incorrect_params();

        //-------------------------------------------------------------
        // E-STEP
        const double new_maxlike = e_step();
        if (std::isnan(new_maxlike))
        {
          OPENMS_LOG_WARN << "Numerical instabilities. Aborting." << endl;
          return false;
        }
        ++itns;

        if (output_plots)
        {
          std::string formula1, formula2, formula3;
          formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "* " + StringUtils::toStr(negative_prior_);
          formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + StringUtils::toStr(negative_prior_) + ")";
          formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
          // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
          file.addLine("plot '" + StringUtils::toStr(param_.getValue("out_plot")) + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        }

        // check termination criterion
        const double increase = new_maxlike - maxlike;
        maxlike = new_maxlike;
        if (increase < pow(10.0, -delta))
        {
          if (increase < -1e-8 * std::fabs(new_maxlike))
          {
            OPENMS_LOG_DEBUG << "Log likelihood of the fit decreased by " << -increase << " (numerical inaccuracy). Stopping." << endl;
          }
          break;
        }
        if (itns >= max_itns)
        {
          OPENMS_LOG_WARN << "Number of iterations exceeded. Convergence criterion not met. Last log likelihood increase: " << increase << endl;
          OPENMS_LOG_WARN << "Algorithm returns probabilities for suboptimal fit. You might want to try raising the max. number of iterations and have a look at the distribution." << endl;
          break;
        }
      }

      //-------------------------------------------------------------
      // Finished fitting
      //-------------------------------------------------------------
      if (!(incorrectMode_() < correctly_assigned_fit_param_.x0))
      {
        OPENMS_LOG_WARN << "Warning: the fitted distribution of incorrect assignments is not located below the one of correct assignments. "
                           "The data might not contain enough correct assignments. Unable to fit the data." << std::endl;
        return false;
      }
      max_incorrectly_ = std::exp(logIncorrectDensity_(incorrectMode_()));
      max_correctly_ = std::exp(logCorrectDensity_(correctly_assigned_fit_param_.x0));

      if (output_plots)
      {
        std::string formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "*" + StringUtils::toStr(negative_prior_);
        std::string formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + StringUtils::toStr(negative_prior_) + ")";
        std::string formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
        // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
        file.addLine("plot '" + StringUtils::toStr(param_.getValue("out_plot")) + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        file.store(StringUtils::toStr(param_.getValue("out_plot")));
        tryGnuplot(StringUtils::toStr(param_.getValue("out_plot")));
      }
      return true;
    }

    bool PosteriorErrorProbabilityModel::fit(std::vector<double>& search_engine_scores, vector<double>& probabilities, const std::string& outlier_handling)
    {
      bool return_value = fit(search_engine_scores, outlier_handling);

      if (!return_value)
      {
        return false;
      }
      probabilities = std::vector<double>(search_engine_scores);
      for (double & p : probabilities)
      {
        p = computeProbability(p);
      }

      return true;
    }

    bool PosteriorErrorProbabilityModel::fitTargetDecoy(const std::vector<double>& target_scores, const std::vector<double>& decoy_scores, const std::string& outlier_handling)
    {
      resetFit_();

      std::vector<double> targets, decoys;
      std::copy_if(target_scores.begin(), target_scores.end(), std::back_inserter(targets), [](double x) { return std::isfinite(x); });
      std::copy_if(decoy_scores.begin(), decoy_scores.end(), std::back_inserter(decoys), [](double x) { return std::isfinite(x); });
      if (targets.size() < 2) { return false; }

      const double decoy_ratio = param_.getValue("decoy_ratio");
      td_incorrect_fraction_ = std::min(1.0, double(decoys.size()) / (decoy_ratio * double(targets.size())));

      const std::string decoy_handling = param_.getValue("decoy_handling").toString();
      if (decoy_handling != "ignore" && !decoys.empty())
      {
        if (double(decoys.size()) >= 0.95 * decoy_ratio * double(targets.size()))
        {
          OPENMS_LOG_WARN << "Warning: there are about as many (or more) decoy as target PSMs (" << decoys.size() << " decoys, " << targets.size()
                          << " targets, decoy_ratio " << decoy_ratio << "). Either the data contain (almost) no correct identifications, or targets and decoys"
                          << " did not compete (e.g. separate target and decoy searches that cannot be matched by spectrum). PEPs will be (close to) 1." << std::endl;
        }
        const double auc = targetDecoyAUC(targets, decoys);
        if (auc <= 0.5)
        {
          OPENMS_LOG_WARN << "Warning: target PSMs do not score better than decoy PSMs (AUC = " << auc
                          << "). Please check the score type and its orientation (higher_score_better)." << std::endl;
        }
      }

      if (!fitTargetDecoyModel_(targets, decoys, outlier_handling))
      {
        return false;
      }

      // Calibration check: under target-decoy competition, the number of decoys scoring at least as well as the targets accepted
      // at a score threshold is (decoy_ratio times) the number of incorrect targets among them, which the model estimates by the
      // sum of their PEPs. Compare observed and expected decoys at the thresholds of 1%, 5% and 10% target-decoy FDR:
      // flag deviations that are both statistically significant (two-sided Poisson test) and larger than a factor of 2.
      if (decoys.size() >= 10 && decoy_handling != "ignore")
      {
        std::vector<double> t_sorted(targets), d_sorted(decoys);
        std::sort(t_sorted.begin(), t_sorted.end(), std::greater<double>());
        std::sort(d_sorted.begin(), d_sorted.end(), std::greater<double>());
        const std::vector<double> fdr_levels = {0.01, 0.05, 0.1};
        std::vector<double> expected_at_level(fdr_levels.size(), 0.0), observed_at_level(fdr_levels.size(), 0.0), size_at_level(fdr_levels.size(), 0.0);
        double sum_pep = 0.0;
        Size d_index = 0;
        for (Size i = 0; i < t_sorted.size(); )
        {
          // process groups of tied target scores together
          Size j = i;
          for (; j < t_sorted.size() && t_sorted[j] == t_sorted[i]; ++j)
          {
            sum_pep += computeProbability(t_sorted[j]);
          }
          while (d_index < d_sorted.size() && d_sorted[d_index] >= t_sorted[i]) { ++d_index; }
          const double n_targets = double(j);
          const double fdr = double(d_index) / decoy_ratio / n_targets;
          for (Size l = 0; l < fdr_levels.size(); ++l)
          {
            if (fdr <= fdr_levels[l]) // keep the largest set below the FDR level
            {
              expected_at_level[l] = sum_pep * decoy_ratio;
              observed_at_level[l] = double(d_index);
              size_at_level[l] = n_targets;
            }
          }
          i = j;
        }
        std::string report;
        bool miscalibrated = false;
        for (Size l = 0; l < fdr_levels.size(); ++l)
        {
          if (size_at_level[l] == 0.0) { continue; }
          const double observed = observed_at_level[l], expected = std::max(expected_at_level[l], 1e-12);
          report += " FDR " + StringUtils::toStr(fdr_levels[l]) + ": " + StringUtils::toStr(int(size_at_level[l])) + " targets, "
                    + StringUtils::toStr(int(observed)) + " decoys, " + StringUtils::toStr(expected) + " expected from PEPs;";
          // two-sided Poisson test (each tail at 0.5e-3)
          const double p_upper = observed > 0.0 ? boost::math::gamma_p(observed, expected) : 1.0; // P(X >= observed)
          const double p_lower = boost::math::gamma_q(observed + 1.0, expected);                // P(X <= observed)
          if ((p_upper < 5e-4 && observed > 2.0 * expected) || (p_lower < 5e-4 && 2.0 * observed < expected))
          {
            miscalibrated = true;
          }
        }
        if (!report.empty())
        {
          OPENMS_LOG_INFO << "Calibration check (decoys above the target-decoy FDR thresholds vs. expected from the PEPs of the accepted targets):" << report << std::endl;
        }
        if (miscalibrated)
        {
          calibration_warning_ = true;
          OPENMS_LOG_WARN << "Warning: the estimated PEPs disagree significantly with the target-decoy estimates:"
                          << report << " The model (" << modelTypeToString(model_type_) << ") might not fit the data well." << std::endl;
        }
      }
      return true;
    }

    bool PosteriorErrorProbabilityModel::fitTargetDecoyModel_(const std::vector<double>& targets, const std::vector<double>& decoys, const std::string& outlier_handling)
    {
      const std::string decoy_handling = param_.getValue("decoy_handling").toString();
      const bool use_gumbel = (param_.getValue("incorrectly_assigned") == "Gumbel");

      if (decoy_handling == "ignore")
      {
        std::vector<double> pooled(targets);
        pooled.insert(pooled.end(), decoys.begin(), decoys.end());
        if (!fitMixture_(pooled, {}, use_gumbel, outlier_handling)) { return false; }
        model_type_ = ModelType::MIXTURE;
        return true;
      }

      const Size min_decoys = static_cast<Size>(int(param_.getValue("min_decoys")));
      const Size min_decoys_isotonic = static_cast<Size>(int(param_.getValue("min_decoys_isotonic")));
      const Size n_decoys = decoys.size();
      const bool enough_decoys = n_decoys >= min_decoys;
      if (!enough_decoys && n_decoys > 0)
      {
        OPENMS_LOG_WARN << "Only " << n_decoys << " decoy PSMs available (min_decoys: " << min_decoys << "). Fitting an unsupervised mixture model to the target scores." << std::endl;
      }

      bool try_isotonic = enough_decoys && (decoy_handling == "isotonic" || (decoy_handling == "auto" && n_decoys >= min_decoys_isotonic));
      if (try_isotonic)
      {
        if (fitIsotonic_(targets, decoys))
        {
          model_type_ = ModelType::TARGET_DECOY_ISOTONIC;
          return true;
        }
        OPENMS_LOG_WARN << "Unable to estimate the non-parametric target-decoy model. Trying the semi-supervised mixture model." << std::endl;
      }

      if (enough_decoys)
      {
        if (fitMixture_(targets, decoys, use_gumbel, outlier_handling))
        {
          model_type_ = ModelType::SEMI_SUPERVISED_MIXTURE;
          if (decoysContradictModel_(targets, decoys))
          {
            if (decoy_handling == "auto" && fitIsotonic_(targets, decoys))
            {
              OPENMS_LOG_WARN << "Using the non-parametric target-decoy estimate instead of the semi-supervised mixture model." << std::endl;
              model_type_ = ModelType::TARGET_DECOY_ISOTONIC;
              return true;
            }
            calibration_warning_ = true;
          }
          const double ratio = negative_prior_ / std::max(td_incorrect_fraction_, 1e-12);
          if (ratio > 1.5 || ratio < 1.0 / 1.5)
          {
            OPENMS_LOG_WARN << "The fraction of incorrect target PSMs estimated by the mixture model (" << negative_prior_
                            << ") differs from the target-decoy estimate (" << td_incorrect_fraction_
                            << "). The model might not fit the data well, or the decoys might not represent incorrect targets (e.g. separate target and decoy searches)." << std::endl;
          }
          return true;
        }
        OPENMS_LOG_WARN << "Unable to fit the semi-supervised mixture model. Fitting an unsupervised mixture model to the target scores." << std::endl;
      }

      if (!fitMixture_(targets, {}, use_gumbel, outlier_handling)) { return false; }
      model_type_ = ModelType::MIXTURE;
      if (!decoys.empty() && decoysContradictModel_(targets, decoys))
      {
        calibration_warning_ = true;
      }
      return true;
    }

    bool PosteriorErrorProbabilityModel::decoysContradictModel_(const std::vector<double>& targets, const std::vector<double>& decoys) const
    {
      // Under target-decoy competition, the number of decoys scoring at least as well as a set of targets
      // is (decoy_ratio times) the number of incorrect targets in the set, which the model estimates by the sum of their PEPs.
      // Test this for the targets that the model considers (fairly) confident (PEP < 0.1): many more decoys than expected
      // indicate a spurious "correct" component, e.g. fitted to the upper tail of data without (many) correct assignments.
      const double decoy_ratio = param_.getValue("decoy_ratio");
      double threshold = std::numeric_limits<double>::infinity();
      double n_targets = 0.0, expected_decoys = 0.0;
      for (double t : targets)
      {
        const double pep = computeProbability(t);
        if (pep < 0.1)
        {
          threshold = std::min(threshold, t);
          n_targets += 1.0;
          expected_decoys += pep;
        }
      }
      if (n_targets == 0.0)
      {
        return false;
      }
      expected_decoys *= decoy_ratio;
      const double observed_decoys = double(std::count_if(decoys.begin(), decoys.end(), [threshold](double d) { return d >= threshold; }));
      if (observed_decoys < 2.0 || observed_decoys <= 2.0 * expected_decoys)
      {
        return false;
      }
      // Poisson upper tail: P(X >= observed | expected) = P(observed, expected) (regularized lower incomplete gamma function)
      const double p_value = boost::math::gamma_p(observed_decoys, std::max(expected_decoys, 1e-12));
      if (p_value >= 1e-3)
      {
        return false;
      }
      OPENMS_LOG_WARN << "Warning: " << observed_decoys << " decoys score at least as well as the " << n_targets
                      << " target PSMs with PEP < 0.1, but the PEPs of these targets imply only " << expected_decoys
                      << " (p = " << p_value << "). The model (" << modelTypeToString(model_type_) << ") is overconfident;"
                      << " the data might contain few or no correct assignments." << std::endl;
      return true;
    }

    bool PosteriorErrorProbabilityModel::fitIsotonic_(const std::vector<double>& target_scores, const std::vector<double>& decoy_scores)
    {
      if (target_scores.empty() || decoy_scores.empty()) { return false; }

      // all PSMs, best scores first
      std::vector<std::pair<double, bool>> psms; // (score, is_decoy)
      psms.reserve(target_scores.size() + decoy_scores.size());
      for (double s : target_scores) { psms.emplace_back(s, false); }
      for (double s : decoy_scores) { psms.emplace_back(s, true); }
      std::sort(psms.begin(), psms.end(), [](const std::pair<double, bool>& a, const std::pair<double, bool>& b) { return a.first > b.first; });

      // 1) maximum likelihood estimate of a decoy probability that does not decrease with decreasing score:
      //    pool-adjacent-violators on groups of tied scores (tied scores always get the same estimate).
      //    Neighbours with equal decoy fractions are merged as well, which yields maximal blocks.
      std::vector<PAVABlock> blocks;
      for (Size i = 0; i < psms.size(); )
      {
        const double score = psms[i].first;
        PAVABlock group{0.0, 0.0, score, score, 0.0};
        for (; i < psms.size() && psms[i].first == score; ++i)
        {
          group.n += 1.0;
          group.decoys += psms[i].second ? 1.0 : 0.0;
          group.score_sum += score;
        }
        addNonDecreasing(blocks, group, true, false);
      }

      // 2) Jeffreys prior on the decoy fraction of each block (avoids PEPs of exactly 0 for blocks without decoys,
      //    e.g. above the best scoring decoy, in particular if such a block is small), followed by pool-adjacent-violators
      //    to restore monotonicity (one pseudo-count per final block)
      std::vector<PAVABlock> shrunk_blocks;
      for (const PAVABlock& b : blocks)
      {
        addNonDecreasing(shrunk_blocks, b, true, true);
      }

      // 3) PEP of a target under target-decoy competition: PEP = (#decoys / decoy_ratio) / #targets = p / (1 - p) / decoy_ratio
      const double decoy_ratio = param_.getValue("decoy_ratio");
      std::vector<double> block_peps; // best scores first
      for (const PAVABlock& b : shrunk_blocks)
      {
        const double p = std::min(b.decoyFraction(true), 1.0 - 1e-12);
        block_peps.push_back(std::min(1.0, p / (1.0 - p) / decoy_ratio));
      }
      // The targets above the best decoy contain on average 1 / decoy_ratio incorrect targets (the number of incorrect targets
      // before the first decoy is geometrically distributed); cf. the "+1" correction of target-decoy FDR estimates.
      if (shrunk_blocks.front().decoys == 0.0)
      {
        block_peps.front() = std::max(block_peps.front(), std::min(1.0, 1.0 / (decoy_ratio * shrunk_blocks.front().n)));
      }
      for (Size b = 1; b < block_peps.size(); ++b)
      {
        block_peps[b] = std::max(block_peps[b], block_peps[b - 1]);
      }

      // Curve: within a block, the PEP decreases linearly with the score around the block's PEP at the mean score of the block,
      // which keeps the mean PEP of the block (calibration) and orders the PSMs of a block by score (no ties of different scores,
      // e.g. for FDR estimates based on the PEPs). The slope is limited such that the PEP stays strictly between the PEPs
      // of the neighbouring blocks. Between blocks, the PEP is interpolated linearly.
      isotonic_scores_.clear();
      isotonic_peps_.clear();
      for (Size b = shrunk_blocks.size(); b-- > 0; ) // ascending scores
      {
        const PAVABlock& block = shrunk_blocks[b];
        const double pep = block_peps[b];
        const double gap_better = (b == 0) ? pep : pep - block_peps[b - 1]; // towards the better scores (PEP must stay >= 0)
        const double gap_worse = (b + 1 == shrunk_blocks.size()) ? 1.0 - pep : block_peps[b + 1] - pep; // towards the worse scores (PEP <= 1)
        const double half_range = 0.45 * std::max(0.0, std::min(gap_better, gap_worse));
        const double center = block.score_sum / block.n;
        const double half_width = std::max(block.max_score - center, center - block.min_score);
        if (half_width > 0.0)
        {
          isotonic_scores_.push_back(block.min_score);
          isotonic_peps_.push_back(pep + half_range * (center - block.min_score) / half_width);
          isotonic_scores_.push_back(block.max_score);
          isotonic_peps_.push_back(pep - half_range * (block.max_score - center) / half_width);
        }
        else
        {
          isotonic_scores_.push_back(block.min_score);
          isotonic_peps_.push_back(pep);
        }
      }

      smallest_score_ = psms.back().first;
      negative_prior_ = td_incorrect_fraction_ >= 0 ? td_incorrect_fraction_ : std::min(1.0, double(decoy_scores.size()) / (decoy_ratio * double(target_scores.size())));
      OPENMS_LOG_INFO << "Non-parametric target-decoy estimate: " << target_scores.size() << " targets, " << decoy_scores.size() << " decoys, "
                      << shrunk_blocks.size() << " blocks, PEP range [" << isotonic_peps_.back() << ", " << isotonic_peps_.front() << "]." << std::endl;
      return !isotonic_scores_.empty();
    }

    std::vector<std::pair<double, double>> PosteriorErrorProbabilityModel::getIsotonicPEPCurve() const
    {
      std::vector<std::pair<double, double>> curve;
      for (Size i = 0; i < isotonic_scores_.size(); ++i)
      {
        curve.emplace_back(isotonic_scores_[i], isotonic_peps_[i]);
      }
      return curve;
    }

    void PosteriorErrorProbabilityModel::fillDensities(const vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      incorrect_density.resize(x_scores.size());
      correct_density.resize(x_scores.size());
      for (Size i = 0; i < x_scores.size(); ++i)
      {
        incorrect_density[i] = std::exp(logIncorrectDensity_(x_scores[i]));
        correct_density[i] = std::exp(logCorrectDensity_(x_scores[i]));
      }
    }

    void PosteriorErrorProbabilityModel::fillLogDensitiesGumbel(const vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      if (incorrect_density.size() != x_scores.size())
      {
        incorrect_density.resize(x_scores.size());
        correct_density.resize(x_scores.size());
      }
      auto incorrect(incorrect_density.begin());
      auto correct(correct_density.begin());
      for (double const & score : x_scores)
      {
        *incorrect = incorrectly_assigned_fit_gumbel_param_.log_eval_no_normalize(score);
        *correct = correctly_assigned_fit_param_.log_eval_no_normalize(score);
        ++incorrect;
        ++correct;
      }
    }

    void PosteriorErrorProbabilityModel::fillLogDensities(const vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      incorrect_density.resize(x_scores.size());
      correct_density.resize(x_scores.size());
      for (Size i = 0; i < x_scores.size(); ++i)
      {
        incorrect_density[i] = logIncorrectDensity_(x_scores[i]);
        correct_density[i] = logCorrectDensity_(x_scores[i]);
      }
    }

    double PosteriorErrorProbabilityModel::computeLogLikelihood(const vector<double>& incorrect_density, const vector<double>& correct_density) const
    {
      double maxlike(0);
      auto incorrect = incorrect_density.cbegin();
      for (auto correct = correct_density.cbegin(); correct < correct_density.cend(); ++correct, ++incorrect)
      {
        maxlike += log10(negative_prior_ * (*incorrect) + (1 - negative_prior_) * (*correct));
      }
      return maxlike;
    }

    double PosteriorErrorProbabilityModel::computeLLAndIncorrectPosteriorsFromLogDensities(
        const vector<double>& incorrect_log_density, const vector<double>& correct_log_density,
        vector<double>& incorrect_posterior) const
    {
      double loglikelihood = 0.0;
      double log_prior_pos = log(1. - negative_prior_);
      double log_prior_neg = log(negative_prior_);
      auto incorrect = incorrect_log_density.cbegin();
      if (incorrect_posterior.size() != incorrect_log_density.size())
      {
        incorrect_posterior.resize(incorrect_log_density.size());
      }
      auto incorrect_posterior_it = incorrect_posterior.begin();

      for (auto correct = correct_log_density.cbegin(); correct < correct_log_density.cend(); ++correct, ++incorrect, ++incorrect_posterior_it)
      {
        double log_resp_correct = log_prior_pos + *correct;
        double log_resp_incorrect = log_prior_neg + *incorrect;
        double max_log_resp = std::max(log_resp_correct,log_resp_incorrect);
        log_resp_correct -= max_log_resp;
        log_resp_incorrect -= max_log_resp;
        double resp_correct = exp(log_resp_correct);
        double resp_incorrect = exp(log_resp_incorrect);
        double sum = resp_correct + resp_incorrect;
        // normalize
        *incorrect_posterior_it = resp_incorrect / sum;
        loglikelihood += max_log_resp + log(sum);
      }
      return loglikelihood;
    }

    std::pair<double,double> PosteriorErrorProbabilityModel::pos_neg_mean_weighted_posteriors(const vector<double>& x_scores, const vector<double>& incorrect_posteriors)
    {
      double pos_x0(0);
      double neg_x0(0);
      auto the_x = x_scores.cbegin();

      for (auto incorrect = incorrect_posteriors.cbegin(); incorrect < incorrect_posteriors.end(); ++incorrect, ++the_x)
      {
        pos_x0 += (1. - *incorrect) * (*the_x);
        neg_x0 += (*incorrect) * (*the_x);
      }
      return {pos_x0,neg_x0};
    }



    std::pair<double,double> PosteriorErrorProbabilityModel::pos_neg_sigma_weighted_posteriors(
        const vector<double>& x_scores,
        const vector<double>& incorrect_posteriors,
        const std::pair<double,double>& pos_neg_mean)
    {
      double pos_sigma(0);
      double neg_sigma(0);
      auto the_x = x_scores.cbegin();

      for (auto incorrect = incorrect_posteriors.cbegin(); incorrect < incorrect_posteriors.end(); ++incorrect, ++the_x)
      {
        pos_sigma += (1. - *incorrect) * pow((*the_x) - pos_neg_mean.first, 2);
        neg_sigma += (*incorrect) * pow((*the_x) - pos_neg_mean.second, 2);
      }
      return {pos_sigma, neg_sigma};
    }

    double PosteriorErrorProbabilityModel::computeProbability(double score) const
    {
      if (model_type_ == ModelType::TARGET_DECOY_ISOTONIC)
      {
        if (std::isnan(score)) { return score; }
        if (score <= isotonic_scores_.front()) { return isotonic_peps_.front(); }
        if (score >= isotonic_scores_.back()) { return isotonic_peps_.back(); }
        const Size upper = std::upper_bound(isotonic_scores_.begin(), isotonic_scores_.end(), score) - isotonic_scores_.begin();
        const double x0 = isotonic_scores_[upper - 1], x1 = isotonic_scores_[upper];
        const double t = (score - x0) / (x1 - x0);
        return isotonic_peps_[upper - 1] + t * (isotonic_peps_[upper] - isotonic_peps_[upper - 1]);
      }

      // apply the same transformation that was applied before fitting
      const double x = score + fabs(smallest_score_) + 0.001;

      // Below the mode of the incorrect component, its density is kept at its maximum,
      // and above the mean of the correct component, its density is kept at its maximum.
      // Together with the mode of the incorrect component being below the mean of the correct one,
      // this ensures that the probability does not increase with the score.
      // (The tails of the Gumbel/Gauss densities would otherwise lead to pathological probabilities.)
      const double log_incorrect = (x < incorrectMode_()) ? logIncorrectDensity_(incorrectMode_()) : logIncorrectDensity_(x);
      const double log_correct = (x > correctly_assigned_fit_param_.x0) ? logCorrectDensity_(correctly_assigned_fit_param_.x0) : logCorrectDensity_(x);
      const double log_odds_correct = std::log(1.0 - negative_prior_) + log_correct - std::log(negative_prior_) - log_incorrect;
      if (std::isnan(log_odds_correct)) { return log_odds_correct; }
      return 1.0 / (1.0 + std::exp(log_odds_correct));
    }

    TextFile PosteriorErrorProbabilityModel::initPlots(vector<double>& x_scores)
    {
      std::vector<DPosition<2> > points;
      Int number_of_bins = param_.getValue("number_of_bins");
      points.resize(number_of_bins);
      DPosition<2> temp;
      double dividing_score = (x_scores.back() - x_scores[0]) / number_of_bins;

      temp.setX(dividing_score / 2);
      temp.setY(0);
      Int bin = 0;
      points[bin] = temp;
      double temp_divider = dividing_score;
      for (std::vector<double>::iterator it = x_scores.begin(); it < x_scores.end(); ++it)
      {
        if (temp_divider - *it >= 0 && bin < number_of_bins - 1)
        {
          points[bin].setY(points[bin].getY() + 1);
        }
        else if (bin  == number_of_bins - 1)
        {
          points[bin].setY(points[bin].getY() + 1);
        }
        else
        {
          temp.setX((temp_divider + temp_divider + dividing_score) / 2);
          temp.setY(1);
          ++bin;
          points[bin] = temp;
          temp_divider += dividing_score;
        }
      }

      TextFile data_points;
      for (DPosition<2>& dp : points)
      {
        dp.setY(dp.getY() / (x_scores.size()  * dividing_score));
        data_points << (StringUtils::toStr(dp.getX()) + "\t" + dp.getY());
      }
      data_points.store(StringUtils::toStr(param_.getValue("out_plot")) + "_scores.txt");

      TextFile file;
      file << "set terminal pdf color solid linewidth 2.0 rounded";
      //file<<"set style empty solid 0.5 border -1";
      //file<<"set style function lines";
      file << "set xlabel \"discriminant score\"";
      file << "set ylabel \"density\"";
      //TODO: file<<"set title ";
      file << "set key off";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file <<  "set output '" + StringUtils::toStr(param_.getValue("out_plot")) + ".pdf'";
      std::string formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "* " + StringUtils::toStr(negative_prior_); //StringUtils::toStr(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + StringUtils::toStr(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + StringUtils::toStr(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + StringUtils::toStr(negative_prior_);
      std::string formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + StringUtils::toStr(negative_prior_) + ")"; //StringUtils::toStr(correctly_assigned_fit_param_.A) +" * exp(-(x - " + StringUtils::toStr(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + StringUtils::toStr(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + StringUtils::toStr(negative_prior_) + ")";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file << ("plot '" + StringUtils::toStr(param_.getValue("out_plot")) + "_scores.txt' with boxes, " + formula1 + " , " + formula2);
      return file;
    }


    //TODO those functions should be members of the Fitter/Function classes!
    const std::string PosteriorErrorProbabilityModel::getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const
    {
      // build a formula with the fitted parameters for gnuplot
      stringstream formula;
      formula << "(1/" << params.sigma << ") * " << "exp(( " << params.x0 << "- x)/" << params.sigma << ") * exp(-exp((" << params.x0 << " - x)/" << params.sigma << "))";
      return formula.str();
    }

    const std::string PosteriorErrorProbabilityModel::getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const
    {
      stringstream formula;
      formula << params.A << " * exp(-(x - " << params.x0 << ") ** 2 / 2 / (" << params.sigma << ") ** 2)";
      return formula.str();
    }

    const std::string PosteriorErrorProbabilityModel::getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect, const GaussFitter::GaussFitResult& correct) const
    {
      stringstream formula;
      formula << negative_prior_ << "*" <<  ((this)->*(getNegativeGnuplotFormula_))(incorrect) << " + (1-" << negative_prior_ << ")*" << ((this)->*(getPositiveGnuplotFormula_))(correct);
      return formula.str();
    }

    void PosteriorErrorProbabilityModel::plotTargetDecoyEstimation(vector<double>& target, vector<double>& decoy)
    {
      if (target.empty() || decoy.empty())
      {
        StringList empty;
        if (target.empty())
        {
          empty.push_back("target");
        }
        if (decoy.empty())
        {
          empty.push_back("decoy");
        }
        OPENMS_LOG_WARN << "Target-Decoy plot was called, but '" << ListUtils::concatenate(empty, "' and '") << "' has no data! Unable to create a target-decoy plot." << std::endl;
        return;
      }
      sort(target.begin(), target.end());
      sort(decoy.begin(), decoy.end());

      // plot in the same (shifted) score space as the mixture model
      const bool mixture = (model_type_ == ModelType::MIXTURE || model_type_ == ModelType::SEMI_SUPERVISED_MIXTURE);
      const double shift = mixture ? fabs(smallest_score_) + 0.001 : 0.0;
      const double lo = std::min(target.front(), decoy.front()) + shift;
      const double hi = std::max(target.back(), decoy.back()) + shift;
      const Int number_of_bins = std::max<Int>(1, param_.getValue("number_of_bins"));
      const double bin_width = (hi > lo) ? (hi - lo) / number_of_bins : 1.0;
      std::vector<double> target_counts(number_of_bins, 0.0), decoy_counts(number_of_bins, 0.0);
      auto bin_of = [&](double x) -> Size
      {
        return std::min<Size>(number_of_bins - 1, static_cast<Size>(std::max(0.0, (x + shift - lo) / bin_width)));
      };
      for (double t : target) { target_counts[bin_of(t)] += 1.0; }
      for (double d : decoy) { decoy_counts[bin_of(d)] += 1.0; }

      // targets: density of the target scores; decoys: scaled to estimate the density of the incorrect target scores
      // (under target-decoy competition, #incorrect targets = #decoys / decoy_ratio)
      const double decoy_ratio = param_.getValue("decoy_ratio");
      const double normalization = double(target.size()) * bin_width;
      TextFile data_points;
      for (Int b = 0; b < number_of_bins; ++b)
      {
        const double center = lo + (b + 0.5) * bin_width;
        data_points << (StringUtils::toStr(center) + "\t" + StringUtils::toStr(target_counts[b] / normalization) + "\t"
                        + StringUtils::toStr(decoy_counts[b] / decoy_ratio / normalization) + "\t"
                        + StringUtils::toStr(computeProbability(center - shift)));
      }
      const std::string out_plot = StringUtils::toStr(param_.getValue("out_plot"));
      data_points.store(out_plot + "_target_decoy_scores.txt");
      TextFile file;
      file << "set terminal pdf color solid linewidth 2.0 rounded";
      file << "set xlabel \"discriminant score\"";
      file << "set ylabel \"density\"";
      file << "set y2label \"posterior error probability\"";
      file << "set y2range [0:1]";
      file << "set y2tics";
      file << "set key off";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file << "set output '" + out_plot + "_target_decoy.pdf'";
      std::string plot = "plot '" + out_plot + "_target_decoy_scores.txt' using 1:2 with boxes, '"
                         + out_plot + "_target_decoy_scores.txt' using 1:3 with boxes fill solid 0.8 noborder, '"
                         + out_plot + "_target_decoy_scores.txt' using 1:4 with lines axes x1y2";
      if (mixture)
      {
        plot += ", " + ((this)->*(getNegativeGnuplotFormula_))(getIncorrectlyAssignedFitResult()) + "* " + StringUtils::toStr(getNegativePrior());
        plot += ", " + getGaussGnuplotFormula(getCorrectlyAssignedFitResult()) + "* (1 - " + StringUtils::toStr(getNegativePrior()) + ")";
      }
      file << plot;
      file.store(out_plot + "_target_decoy");
      if (!mixture)
      {
        file.store(out_plot); // the non-parametric estimate has no other plot
      }
      tryGnuplot(out_plot + "_target_decoy");
    }

    void PosteriorErrorProbabilityModel::tryGnuplot(const std::string& gp_file)
    {
      OPENMS_LOG_INFO << "Attempting to call 'gnuplot' ...";
      std::string cmd =std::string("gnuplot \"") + gp_file + "\"";
      if (system(cmd.c_str()))  // 0 is success!
      {
        OPENMS_LOG_WARN << "Calling 'gnuplot' on '" << gp_file << "' failed. Please create plots manually." << std::endl;
      }
      else OPENMS_LOG_INFO << " success!" << std::endl;

    }

    void PosteriorErrorProbabilityModel::processOutliers_(vector<double>& x_scores, const std::string& outlier_handling) const
    {
      if (x_scores.size() < 3)
      {
        return; // too few values to determine outliers
      }
      if (outlier_handling != "none")
      {
        const std::vector<double> original(x_scores);
        Size nr_outliers = 0;
        Size before = x_scores.size();
        auto q1 = Math::quantile1st(x_scores.begin(),x_scores.end(), true);
        auto q3 = Math::quantile3rd(x_scores.begin(),x_scores.end(), true);
        double iqr = q3 - q1;
        if (outlier_handling == "ignore_iqr_outliers")
        {
          x_scores.erase(
              std::remove_if(x_scores.begin(),x_scores.end(),
                  [&q1,&q3,&iqr](double x){ return x < q1 - 3 * iqr || x > q3 + 3 * iqr;}),
                  x_scores.end()
          );
          nr_outliers = before - x_scores.size();
        }
        else if (outlier_handling == "set_iqr_to_closest_valid")
        {
          auto closest_lower = std::lower_bound(x_scores.begin(), x_scores.end(), q1 - 3 * iqr);
          auto closest_upper = --std::upper_bound(x_scores.begin(), x_scores.end(), q3 + 3 * iqr);

          for (auto it = x_scores.begin(); it != closest_lower; ++it)
          {
            nr_outliers++;
            *it = *closest_lower;
          }
          auto it = closest_upper;
          it++;
          for (; it != x_scores.end(); ++it)
          {
            nr_outliers++;
            *it = *closest_upper;
          }
        }
        else //"ignore_extreme_percentiles"
        {
          Size ninetyninth_idx = std::min<Size>(x_scores.size() - 1, x_scores.size() * 99.9 / 100.);
          double ninetyninth_value = x_scores[ninetyninth_idx];
          Size first_idx = std::min<Size>(x_scores.size() - 1, (x_scores.size() * 1. / 100.) + 1);
          double first_value = x_scores[first_idx];
          x_scores.erase(
              std::remove_if(x_scores.begin(),x_scores.end(),
                             [&first_value,&ninetyninth_value](double x){ return x <= first_value || x >= ninetyninth_value;}),
              x_scores.end()
          );
          nr_outliers = before - x_scores.size();
        }

        if (x_scores.size() < 3)
        {
          OPENMS_LOG_WARN << "Warning: outlier handling '" << outlier_handling << "' would leave fewer than 3 scores (many tied or censored values?). Outliers are kept." << std::endl;
          x_scores = original;
          return;
        }

        double outlier_percent = nr_outliers * 100. / before;
        if (outlier_percent > 2.1)
        {
          OPENMS_LOG_WARN << "Warning: " << outlier_percent << "% outliers detected and corrected. Please double check"
                                                               " the score distribution.\n";
        }
        else
        {
          OPENMS_LOG_INFO << nr_outliers << " outliers detected.\n";
        }
      }
    }

    bool PosteriorErrorProbabilityModel::getScore_(const std::vector<std::string>& requested_score_types, const PeptideHit & hit, const std::string& actual_score_type, double& value)
    {
      for (const auto& requested_score_type : requested_score_types)
      {
        if (actual_score_type == requested_score_type)
        {
          value = hit.getScore();
          return true;
        }
        for (const std::string& key : {requested_score_type, requested_score_type + "_score"})
        {
          if (hit.metaValueExists(key))
          {
            const DataValue& dv = hit.getMetaValue(key);
            if (dv.valueType() == DataValue::DOUBLE_VALUE || dv.valueType() == DataValue::INT_VALUE)
            {
              value = static_cast<double>(dv);
              return true;
            }
            if (dv.valueType() == DataValue::STRING_VALUE)
            {
              try
              {
                value = StringUtils::toDouble(dv.toString());
                return true;
              }
              catch (Exception::ConversionError&)
              {
                // not a number: try the next score type
              }
            }
          }
        }
      }
      return false;
    }

    namespace
    {
      /// dedicated score selection of supported search engines
      struct DedicatedScores
      {
        std::vector<std::string> names; ///< requested score types in order of preference (empty: engine without dedicated support)
        bool neg_log = true;            ///< transform with -log10 (e-values)
        double factor = 1.0;            ///< otherwise: multiply with this factor
      };

      DedicatedScores dedicatedScores(const std::string& engine)
      {
        DedicatedScores d;
        if (engine == "OMSSA") { d.names = {"OMSSA", "expect", "MS:1001328"}; }
        else if (engine == "MYRIMATCH") { d.names = {"mvh", "MS:1001589"}; d.neg_log = false; }
        else if (engine == "XTANDEM") { d.names = {"E-Value", "expect", "MS:1001330"}; }
        else if (engine == "MASCOT") { d.names = {"EValue", "expect", "MS:1001172"}; }
        else if (engine == "SPECTRAST") { d.names = {"f-val", "fval", "MS:1001419"}; d.neg_log = false; d.factor = 100.0; }
        else if (engine == "SIMTANDEM") { d.names = {"E-Value"}; }
        else if (engine == "MSGFPLUS") { d.names = {"MS:1002053", "expect", "SpecEValue", "MS:1002052"}; }
        else if (engine == "COMET") { d.names = {"MS:1002257", "expect"}; }
        else if (engine == "SIMPLESEARCHENGINE" || engine == "PROSE" || engine == "SAGE")
        {
          // HyperScore already returns the log-space value (HyperScore.cpp: log1p(dot_product) + log-factorials)
          // and is labeled "ln(hyperscore)"; "hyperscore" is the older name for the same quantity.
          d.names = {"hyperscore", "ln(hyperscore)"};
          d.neg_log = false;
        }
        else if (engine == "MSFRAGGER") { d.names = {"expect", "MS:1001330"}; }
        return d;
      }
    }

    bool PosteriorErrorProbabilityModel::hasDedicatedScore_(const std::string& engine, const PeptideHit& hit, const std::string& current_score_type)
    {
      const DedicatedScores d = dedicatedScores(engine);
      double value;
      return !d.names.empty() && getScore_(d.names, hit, current_score_type, value);
    }

    std::set<std::string> PosteriorErrorProbabilityModel::runsWithDedicatedScores_(const std::map<std::string, std::string>& run_to_engine, const PeptideIdentificationList& peptide_ids)
    {
      std::set<std::string> runs;
      for (const PeptideIdentification& pep : peptide_ids)
      {
        auto engine_it = run_to_engine.find(pep.getIdentifier());
        if (engine_it == run_to_engine.end() || runs.count(pep.getIdentifier()) > 0)
        {
          continue;
        }
        for (const PeptideHit& hit : pep.getHits())
        {
          if (hasDedicatedScore_(engine_it->second, hit, pep.getScoreType()))
          {
            runs.insert(pep.getIdentifier());
            break;
          }
        }
      }
      return runs;
    }

    double PosteriorErrorProbabilityModel::transformScore_(const std::string & engine, const PeptideHit & hit, const std::string& current_score_type, bool higher_score_better, bool generic_fallback)
    {
      //TODO implement censoring. 1) if value is below censoring take cumulative density below it, instead of point estimate
      const double nan = std::numeric_limits<double>::quiet_NaN();

      if (!generic_fallback)
      {
        // issue #740: unable to fit data with (ion) score 0
        if (engine == "MASCOT" && higher_score_better && hit.getScore() == 0.0)
        {
          return nan;
        }
        const DedicatedScores d = dedicatedScores(engine);
        double value = 0;
        if (!d.names.empty() && getScore_(d.names, hit, current_score_type, value))
        {
          return d.neg_log ? negLog10(value) : d.factor * value;
        }
        return nan; // the score of this search engine is missing for this hit
      }

      // generic fallback: main score of the identification
      const double value = hit.getScore();
      switch (classifyScore(current_score_type))
      {
        case ScoreScale::PROBABILITY:
          return nan; // already a probability (e.g. q-value or PEP): not a meaningful input for the model
        case ScoreScale::EVALUE:
          return higher_score_better ? value : negLog10(value);
        case ScoreScale::RAW:
        default:
          return higher_score_better ? value : -value;
      }
    }

    std::string PosteriorErrorProbabilityModel::getSearchEngineKey(const ProteinIdentification & protein_id)
    {
      std::string search_engine = protein_id.getSearchEngine();
      if (StringUtils::hasPrefix(search_engine, "OpenMS/ConsensusID"))
      {
        if (protein_id.metaValueExists("ConsensusIDBaseSearch"))
        {
          search_engine = StringUtils::prefix(StringUtils::toStr(protein_id.getMetaValue("ConsensusIDBaseSearch")), ':');
        }
        else
        {
          // resolve to the original search engine only if there is exactly one (otherwise, the consensus score is used)
          std::vector<std::string> keys;
          protein_id.getSearchParameters().getKeys(keys);
          const Size n_engines = std::count_if(keys.begin(), keys.end(), [](const std::string& key)
          {
            return StringUtils::hasPrefix(key, "SE:") && !StringUtils::hasSubstring(key, "percolator");
          });
          const std::string original = protein_id.getOriginalSearchEngineName();
          if (n_engines == 1 && original != "Unknown")
          {
            search_engine = original;
          }
        }
      }
      return normalizeEngineKey(search_engine);
    }

    map<std::string, vector<vector<double>>> PosteriorErrorProbabilityModel::extractAndTransformScores(
      const vector<ProteinIdentification> & protein_ids,
      const PeptideIdentificationList & peptide_ids,
      const bool split_charge,
      const bool top_hits_only,
      const bool target_decoy_available,
      const double /* fdr_for_targets_smaller */)
    {
      map<std::string, vector<vector<double>>> all_scores = extractAndTransformScores(protein_ids, peptide_ids, split_charge, top_hits_only);
      for (auto it = all_scores.begin(); it != all_scores.end(); )
      {
        if (it->second[0].size() <= 2) // previous behavior: groups with too few scores for fitting are not reported
        {
          it = all_scores.erase(it);
          continue;
        }
        if (!target_decoy_available)
        {
          it->second[1].clear();
          it->second[2].clear();
        }
        ++it;
      }
      return all_scores;
    }

    map<std::string, vector<vector<double>>> PosteriorErrorProbabilityModel::extractAndTransformScores(
      const vector<ProteinIdentification> & protein_ids,
      const PeptideIdentificationList & peptide_ids,
      const bool split_charge,
      const bool top_hits_only)
    {
      if (split_charge)
      {  // determine different charges in data
        bool has_hits = std::any_of(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& pep) { return !pep.getHits().empty(); });
        if (!has_hits)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'split_charge' is set, but the list of charge states is empty");
        }
      }

      // search engine of each run
      std::map<std::string, std::string> run_to_engine;
      std::map<std::string, StringList> run_to_ms_files; // spectra files of each run (for per-spectrum competition)
      for (const ProteinIdentification& prot : protein_ids)
      {
        run_to_engine[prot.getIdentifier()] = getSearchEngineKey(prot);
        prot.getPrimaryMSRunPath(run_to_ms_files[prot.getIdentifier()]);
      }
      // runs of supported engines use the dedicated score of the engine (if any hit has it), otherwise the main score
      const std::set<std::string> dedicated_runs = runsWithDedicatedScores_(run_to_engine, peptide_ids);

      struct ScoreGroup
      {
        std::vector<double> scores, targets, decoys;
      };
      std::map<std::string, ScoreGroup> groups;
      auto group_key = [split_charge](const std::string& engine, Int charge)
      {
        return split_charge ? engine + "," + StringUtils::toStr(charge) : engine;
      };

      // Top hits (one per spectrum) are used for target-decoy based estimates. If a spectrum has several identifications
      // (e.g. separate target and decoy searches of the same run), the best scoring top hit wins (ties: decoy),
      // which reconstructs the target-decoy competition.
      struct TopHit
      {
        double score;
        bool is_decoy;
      };
      std::map<std::string, std::map<std::string, TopHit>> top_hits; // group key -> spectrum key -> top hit
      Size n_competitions = 0;

      for (Size pep_index = 0; pep_index < peptide_ids.size(); ++pep_index)
      {
        const PeptideIdentification& pep_id = peptide_ids[pep_index];
        auto engine_it = run_to_engine.find(pep_id.getIdentifier());
        if (engine_it == run_to_engine.end() || pep_id.getHits().empty())
        {
          continue;
        }
        const std::string& engine = engine_it->second;
        PeptideIdentification pep(pep_id);
        pep.sort();
        const vector<PeptideHit>& hits = pep.getHits();
        const std::string& score_type = pep.getScoreType();
        const bool higher_score_better = pep.isHigherScoreBetter();
        const bool generic_fallback = (dedicated_runs.count(pep.getIdentifier()) == 0);

        {
          const PeptideHit& top = hits[0];
          const TDLabel label = getTDLabel(top);
          const double score = transformScore_(engine, top, score_type, higher_score_better, generic_fallback);
          if (!std::isnan(score) && label != TDLabel::UNKNOWN)
          {
            std::string spectrum_key;
            const std::string spectrum_reference = pep.getSpectrumReference();
            if (spectrum_reference.empty())
            {
              spectrum_key = "#" + StringUtils::toStr(pep_index); // unknown spectrum: no competition
            }
            else
            {
              // identify the spectrum by its spectra file (e.g. separate target and decoy searches of the same file, stored as separate runs)
              const StringList& ms_files = run_to_ms_files[pep.getIdentifier()];
              Size ms_file_index = 0;
              if (pep.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
              {
                try
                {
                  ms_file_index = static_cast<Size>(std::max(0, StringUtils::toInt32(pep.getMetaValue(Constants::UserParam::ID_MERGE_INDEX).toString())));
                }
                catch (Exception::BaseException&)
                {
                  ms_file_index = 0;
                }
              }
              if (ms_file_index < ms_files.size() && !ms_files[ms_file_index].empty())
              {
                spectrum_key = "file:" + ms_files[ms_file_index] + "\t" + spectrum_reference;
              }
              else
              {
                // no spectra file annotated: identify the spectrum within its identification run
                spectrum_key = "run:" + pep.getIdentifier() + "\t" + pep.getMetaValue("file_origin", std::string()).toString() + "|"
                               + pep.getMetaValue(Constants::UserParam::ID_MERGE_INDEX, std::string()).toString() + "\t" + spectrum_reference;
              }
            }
            const TopHit hit{score, label == TDLabel::DECOY};
            auto [it, inserted] = top_hits[group_key(engine, top.getCharge())].emplace(spectrum_key, hit);
            if (!inserted)
            {
              ++n_competitions;
              if (hit.score > it->second.score || (hit.score == it->second.score && hit.is_decoy))
              {
                it->second = hit;
              }
            }
          }
        }

        // all considered hits: used for the unsupervised model
        const Size n_hits = top_hits_only ? 1 : hits.size();
        for (Size i = 0; i < n_hits; ++i)
        {
          const double score = transformScore_(engine, hits[i], score_type, higher_score_better, generic_fallback);
          if (!std::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
          {
            groups[group_key(engine, hits[i].getCharge())].scores.push_back(score);
          }
        }
      }

      if (n_competitions > 0)
      {
        OPENMS_LOG_INFO << "Found " << n_competitions << " additional identification(s) of already identified spectra (e.g. from separate target and decoy searches). "
                           "Only the best scoring top hit per spectrum is used for target-decoy based estimates." << std::endl;
      }

      bool any_decoy = false;
      for (const auto& [key, spectra] : top_hits)
      {
        ScoreGroup& group = groups[key];
        for (const auto& [spectrum_key, hit] : spectra)
        {
          (hit.is_decoy ? group.decoys : group.targets).push_back(hit.score);
          any_decoy = any_decoy || hit.is_decoy;
        }
      }

      map<std::string, vector<vector<double>>> all_scores;
      for (auto& [key, group] : groups)
      {
        if (!any_decoy)
        {
          group.targets.clear();
        }
        all_scores[key] = {std::move(group.scores), std::move(group.targets), std::move(group.decoys)};
      }
      return all_scores;
    }

    void PosteriorErrorProbabilityModel::updateScores(
      const PosteriorErrorProbabilityModel & PEP_model,
      const std::string & search_engine,
      const Int charge,
      const bool prob_correct,
      const bool split_charge,
      vector<ProteinIdentification> & protein_ids,
      PeptideIdentificationList & peptide_ids,
      bool & unable_to_fit_data,
      bool & data_might_not_be_well_fit)
    {
      std::set<std::pair<Size, Size>> updated_hits;
      updateScores(PEP_model, search_engine, charge, prob_correct, split_charge, protein_ids, peptide_ids, unable_to_fit_data, data_might_not_be_well_fit, updated_hits);
    }

    void PosteriorErrorProbabilityModel::updateScores(
      const PosteriorErrorProbabilityModel & PEP_model,
      const std::string & search_engine,
      const Int charge,
      const bool prob_correct,
      const bool split_charge,
      vector<ProteinIdentification> & protein_ids,
      PeptideIdentificationList & peptide_ids,
      bool & unable_to_fit_data,
      bool & data_might_not_be_well_fit,
      std::set<std::pair<Size, Size>> & updated_hits)
    {
      const std::string engine = normalizeEngineKey(search_engine);
      unable_to_fit_data = false;
      data_might_not_be_well_fit = false;

      std::map<std::string, std::string> runs; // identifier -> engine key
      for (const ProteinIdentification& prot : protein_ids)
      {
        if (getSearchEngineKey(prot) == engine)
        {
          runs[prot.getIdentifier()] = engine;
        }
      }
      // decided before any score is changed (consistent with extractAndTransformScores())
      const std::set<std::string> dedicated_runs = runsWithDedicatedScores_(runs, peptide_ids);
      Size n_unusable = 0;

      for (Size index = 0; index < peptide_ids.size(); ++index)
      {
        PeptideIdentification& pep = peptide_ids[index];
        if (runs.count(pep.getIdentifier()) == 0)
        {
          continue;
        }
        const std::string score_type = pep.getScoreType() + "_score";
        vector<PeptideHit>& hits = pep.getHits();
        for (Size hit_index = 0; hit_index < hits.size(); ++hit_index)
        {
          PeptideHit& hit = hits[hit_index];
          if (!split_charge || (hit.getCharge() == charge))
          {
            double score;
            hit.setMetaValue(score_type, hit.getScore());
            score = PosteriorErrorProbabilityModel::transformScore_(engine, hit, pep.getScoreType(), pep.isHigherScoreBetter(), dedicated_runs.count(pep.getIdentifier()) == 0);

            //TODO they should be ignored during fitting already!
            // and in this issue the -log(10^99) should actually be an acceptable value.
            if (std::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
            {
              score = 1.0;
              ++n_unusable;
            }
            else
            {
              score = PEP_model.computeProbability(score);

              // invalid score? invalid fit!
              if (!(score >= 0.0 && score <= 1.0))
              {
                unable_to_fit_data = true;
              }
              //TODO implement something to check the quality of fit and set data_might_not_be_well_fit
            }
            hit.setScore(prob_correct ? 1.0 - score : score);
            updated_hits.emplace(index, hit_index);
          }
        }
      }
      if (n_unusable > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << n_unusable << " hit(s) of search engine " << engine << (split_charge ? " (charge " + StringUtils::toStr(charge) + ")" : std::string())
                        << " have no usable score and were assigned a PEP of 1." << std::endl;
      }
    }
} // namespace OpenMS // namespace Math
