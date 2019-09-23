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
//
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <QDir>

#include <algorithm>



using namespace std;

namespace OpenMS
{
  namespace Math
  {
    PosteriorErrorProbabilityModel::PosteriorErrorProbabilityModel() :
      DefaultParamHandler("PosteriorErrorProbabilityModel"),
      incorrectly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      incorrectly_assigned_fit_gumbel_param_(GumbelMaxLikelihoodFitter::GumbelDistributionFitResult(-1,-1)),
      correctly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      negative_prior_(0.5), max_incorrectly_(0), max_correctly_(0), smallest_score_(0)
    {
      defaults_.setValue("out_plot", "", "If given, the some output files will be saved in the following manner: <out_plot>_scores.txt for the scores and <out_plot> which contains the fitted values for each step of the EM-algorithm, e.g., out_plot = /usr/home/OMSSA123 leads to /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123 will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123' just OMSSA123, the files will be written into the working directory.", ListUtils::create<String>("advanced,output file"));
      defaults_.setValue("number_of_bins", 100, "Number of bins used for visualization. Only needed if each iteration step of the EM-Algorithm will be visualized", ListUtils::create<String>("advanced"));
      defaults_.setValue("incorrectly_assigned", "Gumbel", "for 'Gumbel', the Gumbel distribution is used to plot incorrectly assigned sequences. For 'Gauss', the Gauss distribution is used.", ListUtils::create<String>("advanced"));
      defaults_.setValue("max_nr_iterations", 1000, "Bounds the number of iterations for the EM algorithm when convergence is slow.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("incorrectly_assigned", ListUtils::create<String>("Gumbel,Gauss"));
      defaults_.setValue("neg_log_delta",6, "The negative logarithm of the convergence threshold for the likelihood increase.");
      defaultsToParam_();
      getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
    }

    PosteriorErrorProbabilityModel::~PosteriorErrorProbabilityModel() = default;

    bool PosteriorErrorProbabilityModel::fitGumbelGauss(std::vector<double>& search_engine_scores)
    {
      // nothing to fit?
      if (search_engine_scores.empty()) { return false; }

      //-------------------------------------------------------------
      // Initializing Parameters
      //-------------------------------------------------------------
      sort(search_engine_scores.begin(), search_engine_scores.end());

      smallest_score_ = search_engine_scores[0];

      vector<double> x_scores{search_engine_scores};

      //transform to a positive range
      for (double & d : x_scores) { d += fabs(smallest_score_) + 0.001; }

      incorrectly_assigned_fit_gumbel_param_.a = Math::mean(x_scores.begin(), x_scores.begin() + ceil(0.5 * x_scores.size())) + x_scores[0];
      incorrectly_assigned_fit_gumbel_param_.b = Math::sd(x_scores.begin(), x_scores.end(), incorrectly_assigned_fit_gumbel_param_.a);
      negative_prior_ = 0.7;

      getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;

      Size x_score_start = std::min(x_scores.size() - 1, (Size) ceil(x_scores.size() * 0.7)); // if only one score is present, ceil(...) will yield 1, which is an invalid index
      correctly_assigned_fit_param_.x0 = Math::mean(x_scores.begin() + x_score_start, x_scores.end()) + x_scores[x_score_start]; //(gauss_scores.begin()->getX() + (gauss_scores.end()-1)->getX())/2;
      correctly_assigned_fit_param_.sigma = incorrectly_assigned_fit_gumbel_param_.b;
      correctly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2.0));

      //-------------------------------------------------------------
      // create files for output
      //-------------------------------------------------------------
      bool output_plots  = (param_.getValue("out_plot").toString().trim().length() > 0);
      TextFile file;
      if (output_plots)
      {
        // create output directory (if not already present)
        QDir dir(param_.getValue("out_plot").toString().toQString());
        if (!dir.cdUp())
        {
          OPENMS_LOG_ERROR << "Could not navigate to output directory for plots from '" << String(dir.dirName()) << "'." << std::endl;
          return false;
        }
        if (!dir.exists() && !dir.mkpath("."))
        {
          OPENMS_LOG_ERROR << "Could not create output directory for plots '" << String(dir.dirName()) << "'." << std::endl;
          return false;
        }
        //
        file = initPlots(x_scores);
      }

      //-------------------------------------------------------------
      // Estimate Parameters - EM algorithm
      //-------------------------------------------------------------
      bool stop_em_init = false;
      Int max_itns = param_.getValue("max_nr_iterations");
      int delta = param_.getValue("neg_log_delta");
      int itns = 0;

      vector<double> incorrect_log_density, correct_log_density;
      fillLogDensitiesGumbel(x_scores, incorrect_log_density, correct_log_density);
      vector<double> bins;
      vector<double> incorrect_posteriors;
      double maxlike = computeLLAndIncorrectPosteriorsFromLogDensities(incorrect_log_density, correct_log_density, incorrect_posteriors);
      double sumIncorrectPosteriors = Math::sum(incorrect_posteriors.begin(),incorrect_posteriors.end());
      double sumCorrectPosteriors = x_scores.size() - sumIncorrectPosteriors;

      OpenMS::Math::GumbelMaxLikelihoodFitter gmlf{incorrectly_assigned_fit_gumbel_param_};

      do
      {
        //-------------------------------------------------------------
        // E-STEP (gauss)
        double newGaussMean = 0.0;
        auto the_x = x_scores.cbegin();

        for (auto incorrect = incorrect_posteriors.cbegin(); incorrect != incorrect_posteriors.cend(); ++incorrect, ++the_x)
        {
          newGaussMean += (1.-*incorrect) * *the_x;
        }
        newGaussMean /= sumCorrectPosteriors;

        double newGaussSigma = 0.0;
        the_x = x_scores.cbegin();

        for (auto incorrect = incorrect_posteriors.cbegin(); incorrect != incorrect_posteriors.cend(); ++incorrect, ++the_x)
        {
          newGaussSigma += (1. - *incorrect) * pow((*the_x) - newGaussMean, 2);
        }
        newGaussSigma = sqrt(newGaussSigma/sumCorrectPosteriors);

        GumbelMaxLikelihoodFitter::GumbelDistributionFitResult newGumbelParams = gmlf.fitWeighted(x_scores, incorrect_posteriors);

        if (newGumbelParams.b <= 0 || std::isnan(newGumbelParams.b))
        {
          OPENMS_LOG_WARN << "Warning: encountered impossible standard deviations. Aborting fit." << std::endl;
          break;
        }

        // update parameters
        correctly_assigned_fit_param_.x0 = newGaussMean;
        correctly_assigned_fit_param_.sigma = newGaussSigma;
        correctly_assigned_fit_param_.A = 1 / sqrt(2 * Constants::PI * pow(newGaussSigma, 2));

        incorrectly_assigned_fit_gumbel_param_ = newGumbelParams;


        // compute new prior probabilities negative peptides
        fillLogDensitiesGumbel(x_scores, incorrect_log_density, correct_log_density);
        double new_maxlike = computeLLAndIncorrectPosteriorsFromLogDensities(incorrect_log_density, correct_log_density, incorrect_posteriors);
        sumIncorrectPosteriors = Math::sum(incorrect_posteriors.begin(),incorrect_posteriors.end());
        sumCorrectPosteriors = x_scores.size() - sumIncorrectPosteriors;
        negative_prior_ = sumIncorrectPosteriors / x_scores.size();

        if (std::isnan(new_maxlike - maxlike)
            || new_maxlike < maxlike)
        {
          return false;
        }

        // check termination criterium
        if ((new_maxlike - maxlike) < pow(10.0, -delta) || itns >= max_itns)
        {
          if (itns >= max_itns)
          {
            OPENMS_LOG_WARN << "Number of iterations exceeded. Convergence criterion not met. Last likelihood increase: " << (new_maxlike - maxlike) << endl;
            OPENMS_LOG_WARN << "Algorithm returns probabilities for suboptimal fit. You might want to try raising the max. number of iterations and have a look at the distribution." << endl;
          }
          stop_em_init = true;
        }

        if (output_plots)
        {
          String formula1, formula2, formula3;
          formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "* " + String(negative_prior_); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
          formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + String(negative_prior_) + ")"; //String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
          formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
          // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
          file.addLine("plot '" + (String)param_.getValue("out_plot") + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        }
        //update maximum likelihood
        maxlike = new_maxlike;
        ++itns;
      } while (!stop_em_init);

      //-------------------------------------------------------------
      // Finished fitting
      //-------------------------------------------------------------
      max_incorrectly_ = getGumbel_(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_);
      max_correctly_ = correctly_assigned_fit_param_.eval(correctly_assigned_fit_param_.x0);

      if (output_plots)
      {
        String formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "*" + String(negative_prior_); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
        String formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + String(negative_prior_) + ")"; // String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
        String formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
        // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
        file.addLine("plot '" + (String)param_.getValue("out_plot") + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        file.store((String)param_.getValue("out_plot"));
        tryGnuplot((String)param_.getValue("out_plot"));
      }
      return true;
    }

    bool PosteriorErrorProbabilityModel::fit(std::vector<double>& search_engine_scores)
    {
      // nothing to fit?
      if (search_engine_scores.empty()) { return false; }

      //-------------------------------------------------------------
      // Initializing Parameters
      //-------------------------------------------------------------
      sort(search_engine_scores.begin(), search_engine_scores.end());

      smallest_score_ = search_engine_scores[0];

      vector<double> x_scores{search_engine_scores};

      //transform to a positive range
      for (double & d : x_scores) { d += fabs(smallest_score_) + 0.001; }

      negative_prior_ = 0.7;
      if (param_.getValue("incorrectly_assigned") == "Gumbel")
      {
        incorrectly_assigned_fit_param_.x0 = Math::mean(x_scores.begin(), x_scores.begin() + ceil(0.5 * x_scores.size())) + x_scores[0];
        incorrectly_assigned_fit_param_.sigma = Math::sd(x_scores.begin(), x_scores.end(), incorrectly_assigned_fit_param_.x0);
        incorrectly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2.0));
        //TODO: Currently, the fit is calculated using the Gauss. 
        getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      }
      else
      {
        incorrectly_assigned_fit_param_.x0 = Math::mean(x_scores.begin(), x_scores.begin() + ceil(0.5 * x_scores.size())) + x_scores[0];
        incorrectly_assigned_fit_param_.sigma = Math::sd(x_scores.begin(), x_scores.end(), incorrectly_assigned_fit_param_.x0);
        incorrectly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2.0));
        getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
      }
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;

      Size x_score_start = std::min(x_scores.size() - 1, (Size) ceil(x_scores.size() * 0.7)); // if only one score is present, ceil(...) will yield 1, which is an invalid index
      correctly_assigned_fit_param_.x0 = Math::mean(x_scores.begin() + x_score_start, x_scores.end()) + x_scores[x_score_start]; //(gauss_scores.begin()->getX() + (gauss_scores.end()-1)->getX())/2;
      correctly_assigned_fit_param_.sigma = incorrectly_assigned_fit_param_.sigma;
      correctly_assigned_fit_param_.A = 1.0 / sqrt(2.0 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2.0));
 
      //-------------------------------------------------------------
      // create files for output
      //-------------------------------------------------------------
      bool output_plots  = (param_.getValue("out_plot").toString().trim().length() > 0);
      TextFile file;
      if (output_plots)
      {
        // create output directory (if not already present)
        QDir dir(param_.getValue("out_plot").toString().toQString());
        if (!dir.cdUp())
        {
          OPENMS_LOG_ERROR << "Could not navigate to output directory for plots from '" << String(dir.dirName()) << "'." << std::endl;
          return false;
        }
        if (!dir.exists() && !dir.mkpath("."))
        {
          OPENMS_LOG_ERROR << "Could not create output directory for plots '" << String(dir.dirName()) << "'." << std::endl;
          return false;
        }
        //
        file = initPlots(x_scores);
      }

      //-------------------------------------------------------------
      // Estimate Parameters - EM algorithm
      //-------------------------------------------------------------
      bool stop_em_init = false;
      Int max_itns = param_.getValue("max_nr_iterations");
      int delta = param_.getValue("neg_log_delta");
      int itns = 0;

      vector<double> incorrect_log_density, correct_log_density;
      fillLogDensities(x_scores, incorrect_log_density, correct_log_density);
      vector<double> incorrect_posteriors;
      double maxlike = computeLLAndIncorrectPosteriorsFromLogDensities(incorrect_log_density, correct_log_density, incorrect_posteriors);
      double sumIncorrectPosteriors = Math::sum(incorrect_posteriors.begin(),incorrect_posteriors.end());
      double sumCorrectPosteriors = x_scores.size() - sumIncorrectPosteriors;

      do
      {
        //-------------------------------------------------------------
        // E-STEP
        std::pair<double,double> newMeans = pos_neg_mean_weighted_posteriors(x_scores, incorrect_posteriors);
        newMeans.first /= sumCorrectPosteriors;
        newMeans.second /= sumIncorrectPosteriors;

        //new standard deviation
        std::pair<double,double> newSigmas = pos_neg_sigma_weighted_posteriors(x_scores, incorrect_posteriors, newMeans);
        newSigmas.first = sqrt(newSigmas.first/sumCorrectPosteriors);
        newSigmas.second = sqrt(newSigmas.second/sumIncorrectPosteriors);

        if (newSigmas.first <= 0 || newSigmas.second <= 0 || std::isnan(newSigmas.first) || std::isnan(newSigmas.second) )
        {
          OPENMS_LOG_WARN << "Warning: encountered impossible standard deviations. Aborting fit." << std::endl;
          break;
        }

        // update parameters
        correctly_assigned_fit_param_.x0 = newMeans.first;
        incorrectly_assigned_fit_param_.x0 = newMeans.second;

        correctly_assigned_fit_param_.sigma = newSigmas.first;
        correctly_assigned_fit_param_.A = 1 / sqrt(2 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2));

        incorrectly_assigned_fit_param_.sigma = newSigmas.second;
        incorrectly_assigned_fit_param_.A = 1 / sqrt(2 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2));


        // compute new prior probabilities negative peptides
        fillLogDensities(x_scores, incorrect_log_density, correct_log_density);
        double new_maxlike = computeLLAndIncorrectPosteriorsFromLogDensities(incorrect_log_density, correct_log_density, incorrect_posteriors);
        sumIncorrectPosteriors = Math::sum(incorrect_posteriors.begin(),incorrect_posteriors.end());
        sumCorrectPosteriors = x_scores.size() - sumIncorrectPosteriors;
        negative_prior_ = sumIncorrectPosteriors / x_scores.size();

        if (std::isnan(new_maxlike - maxlike) 
          || new_maxlike < maxlike)
        {
          return false;
        }

        // check termination criterium
        if ((new_maxlike - maxlike) < pow(10.0, -delta) || itns >= max_itns)
        {
          if (itns >= max_itns)
          {
            OPENMS_LOG_WARN << "Number of iterations exceeded. Convergence criterion not met. Last likelihood increase: " << (new_maxlike - maxlike) << endl;
            OPENMS_LOG_WARN << "Algorithm returns probabilities for suboptimal fit. You might want to try raising the max. number of iterations and have a look at the distribution." << endl;
          }
          stop_em_init = true;
        }

        if (output_plots)
        {
          String formula1, formula2, formula3;
          formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "* " + String(negative_prior_); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
          formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + String(negative_prior_) + ")"; //String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
          formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
          // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
          file.addLine("plot '" + (String)param_.getValue("out_plot") + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        }
        //update maximum likelihood
        maxlike = new_maxlike;
        ++itns;
      } while (!stop_em_init);

      //-------------------------------------------------------------
      // Finished fitting
      //-------------------------------------------------------------
      if (param_.getValue("incorrectly_assigned") == "Gumbel")
      {
        max_incorrectly_ = getGumbel_(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_);
      }
      else
      {
        max_incorrectly_ = incorrectly_assigned_fit_param_.eval(incorrectly_assigned_fit_param_.x0);
      }
      max_correctly_ = correctly_assigned_fit_param_.eval(correctly_assigned_fit_param_.x0);

      if (output_plots)
      {
        String formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "*" + String(negative_prior_); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
        String formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + String(negative_prior_) + ")"; // String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
        String formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_, correctly_assigned_fit_param_);
        // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
        file.addLine("plot '" + (String)param_.getValue("out_plot") + "_scores.txt' with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
        file.store((String)param_.getValue("out_plot"));
        tryGnuplot((String)param_.getValue("out_plot"));
      }
      return true;
    }

    bool PosteriorErrorProbabilityModel::fit(std::vector<double>& search_engine_scores, vector<double>& probabilities)
    {
      bool return_value = fit(search_engine_scores);

      if (!return_value) return false;

      probabilities = std::vector<double>(search_engine_scores);
      for (double & p : probabilities) { p = computeProbability(p); }

      return true;
    }

    void PosteriorErrorProbabilityModel::fillDensities(const vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
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
        // TODO: incorrect is currently filled with gauss as fitting gumble is not supported
        *incorrect = incorrectly_assigned_fit_param_.eval(score);
        *correct = correctly_assigned_fit_param_.eval(score);
        ++incorrect;
        ++correct;
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
      if (incorrect_density.size() != x_scores.size())
      {
        incorrect_density.resize(x_scores.size());
        correct_density.resize(x_scores.size());
      }
      auto incorrect(incorrect_density.begin());
      auto correct(correct_density.begin());
      for (double const & score : x_scores)
      {
        // TODO: incorrect is currently filled with gauss as fitting gumble is not supported
        *incorrect = incorrectly_assigned_fit_param_.log_eval_no_normalize(score);
        *correct = correctly_assigned_fit_param_.log_eval_no_normalize(score);
        ++incorrect;
        ++correct;
      }
    }

    double PosteriorErrorProbabilityModel::computeLogLikelihood(const vector<double>& incorrect_density, const vector<double>& correct_density)
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
        vector<double>& incorrect_posterior)
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
        *incorrect_posterior_it = resp_incorrect / sum; //TODO can we somehow stay in log space (i.e. fill as log posteriors?)
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
      // apply the same transformation that was applied before fitting
      score = score + fabs(smallest_score_) + 0.001;
      double x_neg, x_pos;

      // the score is smaller than the peak of incorrectly assigned sequences.
      // To ensure that the probabilities wont rise again use the incorrectly assigned peak for computation
      if (score < incorrectly_assigned_fit_param_.x0)
      {
        x_neg = max_incorrectly_;
        x_pos = correctly_assigned_fit_param_.eval(score);
      }
      // same as above. However, this time to ensure that probabilities wont drop again.
      //TODO this does not consider the possibility of using Gauss as negative function anymore! Confusing at best!
      else if (score > correctly_assigned_fit_param_.x0)
      {
        x_neg = getGumbel_(score, incorrectly_assigned_fit_param_);
        x_pos = max_correctly_;
      }
      // if it's in-between use the normal formula
      else
      {
        x_neg = getGumbel_(score, incorrectly_assigned_fit_param_);
        x_pos = correctly_assigned_fit_param_.eval(score);
      }
      return (negative_prior_ * x_neg) / ((negative_prior_ * x_neg) + (1 - negative_prior_) * x_pos);
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
      for (vector<DPosition<2> >::iterator it = points.begin(); it < points.end(); ++it)
      {
        it->setY(it->getY() / (x_scores.size()  * dividing_score));
        data_points << (String(it->getX()) + "\t" + it->getY());
      }
      data_points.store((String)param_.getValue("out_plot") + "_scores.txt");

      TextFile file;
      file << "set terminal pdf color solid linewidth 2.0 rounded";
      //file<<"set style empty solid 0.5 border -1";
      //file<<"set style function lines";
      file << "set xlabel \"discriminant score\"";
      file << "set ylabel \"density\"";
      //TODO: file<<"set title ";
      file << "set key off";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file <<  "set output '" + (String)param_.getValue("out_plot") + ".pdf'";
      String formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_) + "* " + String(negative_prior_); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
      String formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_) + "* (1 - " + String(negative_prior_) + ")"; //String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file << ("plot '" + (String)param_.getValue("out_plot") + "_scores.txt' with boxes, " + formula1 + " , " + formula2);
      return file;
    }


    //TODO those functions should be members of the Fitter/Function classes!
    const String PosteriorErrorProbabilityModel::getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const
    {
      // build a formula with the fitted parameters for gnuplot
      stringstream formula;
      formula << "(1/" << params.sigma << ") * " << "exp(( " << params.x0 << "- x)/" << params.sigma << ") * exp(-exp((" << params.x0 << " - x)/" << params.sigma << "))";
      return formula.str();
    }

    const String PosteriorErrorProbabilityModel::getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const
    {
      stringstream formula;
      formula << params.A << " * exp(-(x - " << params.x0 << ") ** 2 / 2 / (" << params.sigma << ") ** 2)";
      return formula.str();
    }

    const String PosteriorErrorProbabilityModel::getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect, const GaussFitter::GaussFitResult& correct) const
    {
      stringstream formula;
      formula << negative_prior_ << "*" <<  ((this)->*(getNegativeGnuplotFormula_))(incorrect) << " + (1-" << negative_prior_ << ")*" << ((this)->*(getPositiveGnuplotFormula_))(correct);
      return formula.str();
    }

    void PosteriorErrorProbabilityModel::plotTargetDecoyEstimation(vector<double>& target, vector<double>& decoy)
    {
      if (target.size() == 0 || decoy.size() == 0)
      {
        StringList empty;
        if (target.size() == 0) empty.push_back("target");
        if (decoy.size() == 0) empty.push_back("decoy");
        OPENMS_LOG_WARN << "Target-Decoy plot was called, but '" << ListUtils::concatenate(empty, "' and '") << "' has no data! Unable to create a target-decoy plot." << std::endl;
        return;
      }
      Int number_of_bins = param_.getValue("number_of_bins");
      std::vector<DPosition<3> > points(number_of_bins);
      DPosition<3> temp;

      sort(target.begin(), target.end());
      sort(decoy.begin(), decoy.end());

      double dividing_score = (max(target.back(), decoy.back()) /*scores.back()*/ - min(target[0], decoy[0]) /*scores[0]*/) / number_of_bins;

      temp[0] = (dividing_score / 2);
      temp[1] = 0;
      temp[2] = 0;
      Int bin = 0;
      points[bin] = temp;
      double temp_divider = dividing_score;
      for (std::vector<double>::iterator it = target.begin(); it < target.end(); ++it)
      {
        *it = *it + fabs(smallest_score_) + 0.001;
        if (temp_divider - *it >= 0 && bin < number_of_bins - 1)
        {
          points[bin][1] += 1;
        }
        else if (bin  == number_of_bins - 1)
        {
          points[bin][1] += 1;
        }
        else
        {
          temp[0] = ((temp_divider + temp_divider + dividing_score) / 2);
          temp[1] = 1;
          ++bin;
          points[bin] = temp;
          temp_divider += dividing_score;
        }
      }

      bin = 0;
      temp_divider = dividing_score;
      for (std::vector<double>::iterator it = decoy.begin(); it < decoy.end(); ++it)
      {
        *it = *it + fabs(smallest_score_) + 0.001;
        if (temp_divider - *it >= 0 && bin < number_of_bins - 1)
        {
          points[bin][2] += 1;
        }
        else if (bin  == number_of_bins - 1)
        {
          points[bin][2] += 1;
        }
        else
        {
          // temp[0] = ((temp_divider + temp_divider + dividing_score)/2);
          // temp[2] = 1;
          ++bin;
          points[bin][2] = 1;
          temp_divider += dividing_score;
        }
      }

      TextFile data_points;
      for (vector<DPosition<3> >::iterator it = points.begin(); it < points.end(); ++it)
      {
        (*it)[1] = ((*it)[1] / ((decoy.size() + target.size())  * dividing_score));
        (*it)[2] = ((*it)[2] / ((decoy.size() + target.size())  * dividing_score));
        String temp_ = (*it)[0];
        temp_ += "\t";
        temp_ += (*it)[1];
        temp_ += "\t";
        temp_ += (*it)[2];
        data_points << temp_;
      }
      data_points.store((String)param_.getValue("out_plot") + "_target_decoy_scores.txt");
      TextFile file;
      file << "set terminal pdf color solid linewidth 2.0 rounded";
      //file<<"set style empty solid 0.5 border -1";
      //file<<"set style function lines";
      file << "set xlabel \"discriminant score\"";
      file << "set ylabel \"density\"";
      //TODO: file<<"set title ";
      file << "set key off";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file << String("set output '") +  (String)param_.getValue("out_plot") + "_target_decoy.pdf'";
      String formula1, formula2;
      formula1 = getGumbelGnuplotFormula(getIncorrectlyAssignedFitResult()) + "* " + String(getNegativePrior()); //String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
      formula2 = getGaussGnuplotFormula(getCorrectlyAssignedFitResult()) + "* (1 - " + String(getNegativePrior()) + ")"; //String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
      // important: use single quotes for paths, since otherwise backslashes will not be accepted on Windows!
      file << ("plot '" + (String)param_.getValue("out_plot") + "_target_decoy_scores.txt'   using 1:3  with boxes fill solid 0.8 noborder, \"" + (String)param_.getValue("out_plot") + "_target_decoy_scores.txt\"  using 1:2  with boxes, " + formula1 + " , " + formula2);
      file.store((String)param_.getValue("out_plot") + "_target_decoy");
      tryGnuplot((String)param_.getValue("out_plot") + "_target_decoy");
    }

    void PosteriorErrorProbabilityModel::tryGnuplot(const String& gp_file)
    {
      OPENMS_LOG_INFO << "Attempting to call 'gnuplot' ...";
      String cmd = String("gnuplot \"") + gp_file + "\"";
      if (system(cmd.c_str()))  // 0 is success!
      {
        OPENMS_LOG_WARN << "Calling 'gnuplot' on '" << gp_file << "' failed. Please create plots manually." << std::endl;
      }
      else OPENMS_LOG_INFO << " success!" << std::endl;

    }

    double PosteriorErrorProbabilityModel::transformScore_(const String & engine, const PeptideHit & hit)
    {
      //TODO implement censoring. 1) if value is below censoring take cumulative density below it, instead of point estimate
      //Set fixed e-value threshold
      //rather take min to be safe even with optimizations
      //const double smallest_e_value_ = numeric_limits<double>::denorm_min();
      const double smallest_e_value_ = numeric_limits<double>::min();

      //TODO we don't care about score types here? What if the data was processed with
      // IDPEP or Percolator already?
      if (engine == "OMSSA")
      {
        return (-1) * log10(max(hit.getScore(), smallest_e_value_));
      }
      else if (engine == "MYRIMATCH" ) 
      {
        return hit.getScore();
      }
      else if (engine == "XTANDEM")
      {
        return (-1) * log10(max((double)hit.getMetaValue("E-Value"), smallest_e_value_));
      }
      else if (engine == "MASCOT")
      {
        // issue #740: unable to fit data with score 0
        if (hit.getScore() == 0.0) 
        {
          return numeric_limits<double>::quiet_NaN();
        }
        // end issue #740
        if (hit.metaValueExists("EValue"))
        {
          return (-1) * log10(max((double)hit.getMetaValue("EValue"), smallest_e_value_));
        }
        if (hit.metaValueExists("expect"))
        {
          return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
        }
      }
      else if (engine == "SPECTRAST")
      {
        return 100 * hit.getScore(); // f-val
      }
      else if (engine == "SIMTANDEM")
      {
        if (hit.metaValueExists("E-Value"))
        {
          return (-1) * log10(max((double)hit.getMetaValue("E-Value"), smallest_e_value_));
        }
      }
      else if ((engine == "MSGFPLUS") || (engine == "MS-GF+"))
      {
        if (hit.metaValueExists("MS:1002053"))  // name: MS-GF:EValue
        {
          return (-1) * log10(max((double)hit.getMetaValue("MS:1002053"), smallest_e_value_));
        }
        else if (hit.metaValueExists("expect"))
        {
          return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
        }
      }
      else if (engine == "COMET")
      {
        if (hit.metaValueExists("MS:1002257")) // name: Comet:expectation value
        {
          return (-1) * log10(max((double)hit.getMetaValue("MS:1002257"), smallest_e_value_));
        }
        else if (hit.metaValueExists("expect"))
        {
          return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
        }
      }

      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No parameters for chosen search engine", "The chosen search engine is currently not supported");
    }

    map<String, vector<vector<double>>> PosteriorErrorProbabilityModel::extractAndTransformScores(
      const vector<ProteinIdentification> & protein_ids,
      const vector<PeptideIdentification> & peptide_ids,
      const bool split_charge,
      const bool top_hits_only,
      const bool target_decoy_available,
      const double fdr_for_targets_smaller)
    {
      std::set<Int> charges;
      const StringList search_engines = ListUtils::create<String>("XTandem,OMSSA,MASCOT,SpectraST,MyriMatch,SimTandem,MSGFPlus,MS-GF+,Comet,tide-search,OpenMS/ConsensusID_best,OpenMS/ConsensusID_worst,OpenMS/ConsensusID_average");

      if (split_charge)
      {  // determine different charges in data
        for (PeptideIdentification const & pep_id : peptide_ids)
        {
          const vector<PeptideHit>& hits = pep_id.getHits();
          for (PeptideHit const & hit : hits) { charges.insert(hit.getCharge()); }
        }
        if (charges.empty())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'split_charge' is set, but the list of charge states is empty");
        }
      }

      set<Int>::iterator charge_it = charges.begin(); // charges can be empty, no problem if split_charge is not set
      map<String, vector<vector<double> > > all_scores;
      char splitter = ','; // to split the engine from the charge state later on
      do
      {
        vector<double> scores, decoy, target;
        for (String supported_engine : search_engines)
        {
          supported_engine.toUpper();
          for (ProteinIdentification const & prot : protein_ids)
          {
            String search_engine = prot.getSearchEngine();
            if (search_engine.hasPrefix("OpenMS/ConsensusID"))
            {
              search_engine = prot.getMetaValue("ConsensusIDBaseSearch");
              search_engine = search_engine.prefix(':');
            }
            search_engine.toUpper();

            if (supported_engine == search_engine)
            {
              for (PeptideIdentification pep : peptide_ids)
              {
                // make sure we are comparing peptide and proteins of the same search run
                if (prot.getIdentifier() == pep.getIdentifier())
                {
                  pep.sort();
                  vector<PeptideHit>& hits = pep.getHits();
                  if (top_hits_only)
                  {
                    if (!hits.empty() && (!split_charge || hits[0].getCharge() == *charge_it))
                    {
                      double score = PosteriorErrorProbabilityModel::transformScore_(supported_engine, hits[0]);
                      if (!std::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                      {
                        scores.push_back(score);

                        if (target_decoy_available)
                        {
                          if (hits[0].getScore() < fdr_for_targets_smaller)
                          {
                            target.push_back(score);
                          }
                          else
                          {
                            decoy.push_back(score);
                          }
                        }
                      }
                    }
                  }
                  else
                  {
                    for (PeptideHit const & hit : hits)
                    {
                      if (!split_charge || (hit.getCharge() == *charge_it))
                      {
                        double score = PosteriorErrorProbabilityModel::transformScore_(supported_engine, hit);
                        if (!std::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                        {
                          scores.push_back(score);
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          if (scores.size() > 2)
          {
            vector<vector<double> > tmp;
            tmp.push_back(scores);
            tmp.push_back(target);
            tmp.push_back(decoy);

            if (split_charge)
            {
              String engine_with_charge_state = supported_engine + String(splitter) + String(*charge_it);
              all_scores.insert(make_pair(engine_with_charge_state, tmp));
            }
            else
            {
              all_scores.insert(make_pair(supported_engine, tmp));
            }
          }

          scores.clear();
          target.clear();
          decoy.clear();
        }

        if (split_charge) { ++charge_it; }
      } while (charge_it != charges.end());
      return all_scores;
    }


    void PosteriorErrorProbabilityModel::updateScores(
      const PosteriorErrorProbabilityModel & PEP_model,
      const String & search_engine,
      const Int charge,
      const bool prob_correct,
      const bool split_charge,
      vector<ProteinIdentification> & protein_ids,
      vector<PeptideIdentification> & peptide_ids,
      bool & unable_to_fit_data,
      bool & data_might_not_be_well_fit)
    {
      String engine(search_engine);
      unable_to_fit_data = false;
      data_might_not_be_well_fit = false;

      engine.toUpper();
      for (ProteinIdentification & prot : protein_ids)
      {
        String se = prot.getSearchEngine();
        se.toUpper();

        if (engine == se)
        {
          for (PeptideIdentification & pep : peptide_ids)
          {
            if (prot.getIdentifier() == pep.getIdentifier())
            {
              String score_type = pep.getScoreType() + "_score";
              vector<PeptideHit> hits = pep.getHits();
              for (PeptideHit & hit : hits)
              {
                if (!split_charge || (hit.getCharge() == charge))
                {
                  double score;
                  hit.setMetaValue(score_type, hit.getScore());
                  score = PosteriorErrorProbabilityModel::transformScore_(engine, hit);

                  //TODO they should be ignored during fitting already!
                  // and in this issue the -log(10^99) should actually be an acceptable value.
                  if (std::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                  {
                    score = 1.0;
                  }
                  else 
                  { 
                    score = PEP_model.computeProbability(score);

                    // invalid score? invalid fit!
                    if ((score < 0.0) || (score > 1.0)) unable_to_fit_data = true;
                    //TODO implement something to check the quality of fit and set data_might_not_be_well_fit
                  }
                  hit.setScore(score);
                  if (prob_correct)
                  {
                    hit.setScore(1.0 - score);
                  }
                  else
                  {
                    hit.setScore(score);
                  }
                }
              }
              pep.setHits(hits);
            }
            if (prob_correct)
            {
              pep.setScoreType("Posterior Probability");
              pep.setHigherScoreBetter(true);
            }
            else
            {
              pep.setScoreType("Posterior Error Probability");
              pep.setHigherScoreBetter(false);
            }
          }
        }
      }
    }
  } // namespace Math
} // namespace OpenMS
