// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <QDir>

#include <boost/math/special_functions/fpclassify.hpp>

#include <algorithm>



using namespace std;

namespace OpenMS
{
  namespace Math
  {
    PosteriorErrorProbabilityModel::PosteriorErrorProbabilityModel() :
      DefaultParamHandler("PosteriorErrorProbabilityModel"),
      incorrectly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      correctly_assigned_fit_param_(GaussFitter::GaussFitResult(-1, -1, -1)),
      negative_prior_(0.5), max_incorrectly_(0), max_correctly_(0), smallest_score_(0)
    {
      defaults_.setValue("out_plot", "", "If given, the some output files will be saved in the following manner: <out_plot>_scores.txt for the scores and <out_plot> which contains the fitted values for each step of the EM-algorithm, e.g., out_plot = /usr/home/OMSSA123 leads to /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123 will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123' just OMSSA123, the files will be written into the working directory.", ListUtils::create<String>("advanced,output file"));
      defaults_.setValue("number_of_bins", 100, "Number of bins used for visualization. Only needed if each iteration step of the EM-Algorithm will be visualized", ListUtils::create<String>("advanced"));
      defaults_.setValue("incorrectly_assigned", "Gumbel", "for 'Gumbel', the Gumbel distribution is used to plot incorrectly assigned sequences. For 'Gauss', the Gauss distribution is used.", ListUtils::create<String>("advanced"));
      defaults_.setValue("max_nr_iterations", 1000, "Bounds the number of iterations for the EM algorithm when convergence is slow.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("incorrectly_assigned", ListUtils::create<String>("Gumbel,Gauss"));
      defaultsToParam_();
      calc_incorrect_ = &PosteriorErrorProbabilityModel::getGumbel;
      calc_correct_ = &PosteriorErrorProbabilityModel::getGauss;
      getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
    }

    PosteriorErrorProbabilityModel::~PosteriorErrorProbabilityModel()
    {
    }

    bool PosteriorErrorProbabilityModel::fit(std::vector<double>& search_engine_scores)
    {
      if (search_engine_scores.empty())
      {
        return false;
      }
      //-------------------------------------------------------------
      // Initializing Parameters
      //-------------------------------------------------------------
      sort(search_engine_scores.begin(), search_engine_scores.end());

      smallest_score_ = search_engine_scores[0];
      vector<double> x_scores;
      x_scores.resize(search_engine_scores.size());
      std::vector<double>::iterator it = x_scores.begin();
      for (std::vector<double>::iterator iti = search_engine_scores.begin(); iti < search_engine_scores.end(); ++it, ++iti)
      {
        *it = *iti + fabs(smallest_score_) + 0.001;
      }
      negative_prior_ = 0.7;
      if (param_.getValue("incorrectly_assigned") == "Gumbel")
      {
        incorrectly_assigned_fit_param_.x0 = Math::mean(x_scores.begin(), x_scores.begin() + ceil(0.5 * x_scores.size())) + x_scores[0];
        incorrectly_assigned_fit_param_.sigma = Math::sd(x_scores.begin(), x_scores.end(), incorrectly_assigned_fit_param_.x0);
        incorrectly_assigned_fit_param_.A = 1   / sqrt(2 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2));
        //TODO: compute directly with gauss. Workaround:
        calc_incorrect_ = &PosteriorErrorProbabilityModel::getGauss;
        getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
      }
      else
      {
        incorrectly_assigned_fit_param_.x0 = Math::mean(x_scores.begin(), x_scores.begin() + ceil(0.5 * x_scores.size())) + x_scores[0];
        incorrectly_assigned_fit_param_.sigma = Math::sd(x_scores.begin(), x_scores.end(), incorrectly_assigned_fit_param_.x0);
        incorrectly_assigned_fit_param_.A = 1   / sqrt(2 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2));
        calc_incorrect_ = &PosteriorErrorProbabilityModel::getGauss;
        getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
      }
      getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
      calc_correct_ = &PosteriorErrorProbabilityModel::getGauss;
      Size x_score_start = std::min(x_scores.size() - 1, (Size) ceil(x_scores.size() * 0.7)); // if only one score is present, ceil(...) will yield 1, which is an invalid index
      correctly_assigned_fit_param_.x0 = Math::mean(x_scores.begin() + x_score_start, x_scores.end()) + x_scores[x_score_start]; //(gauss_scores.begin()->getX() + (gauss_scores.end()-1)->getX())/2;
      correctly_assigned_fit_param_.sigma = incorrectly_assigned_fit_param_.sigma;
      correctly_assigned_fit_param_.A = 1.0   / sqrt(2 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2));

      vector<double> incorrect_density;
      vector<double> correct_density;
      fillDensities(x_scores, incorrect_density, correct_density);

      double maxlike = computeMaxLikelihood(incorrect_density, correct_density);
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
          LOG_ERROR << "Could not navigate to output directory for plots from '" << String(dir.dirName()) << "'." << std::endl;
          return false;
        }
        if (!dir.exists() && !dir.mkpath("."))
        {
          LOG_ERROR << "Could not create output directory for plots '" << String(dir.dirName()) << "'." << std::endl;
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
      int delta = 6;
      int itns = 0;
      
      do
      {
        //E-STEP
        double one_minus_sum_posterior = one_minus_sum_post(incorrect_density, correct_density);
        double sum_posterior = sum_post(incorrect_density, correct_density);

        //new mean
        double sum_positive_x0 = sum_pos_x0(x_scores, incorrect_density, correct_density);
        double sum_negative_x0 = sum_neg_x0(x_scores, incorrect_density, correct_density);

        double positive_mean = sum_positive_x0 / one_minus_sum_posterior;
        double negative_mean = sum_negative_x0 / sum_posterior;

        //new standard deviation
        double sum_positive_sigma = sum_pos_sigma(x_scores, incorrect_density, correct_density, positive_mean);
        double sum_negative_sigma = sum_neg_sigma(x_scores, incorrect_density, correct_density, negative_mean);

        //update parameters
        correctly_assigned_fit_param_.x0 = positive_mean;
        if (sum_positive_sigma  != 0)
        {
          correctly_assigned_fit_param_.sigma = sqrt(sum_positive_sigma / one_minus_sum_posterior);
          correctly_assigned_fit_param_.A = 1 / sqrt(2 * Constants::PI * pow(correctly_assigned_fit_param_.sigma, 2));
        }

        incorrectly_assigned_fit_param_.x0 = negative_mean;
        if (sum_negative_sigma  != 0)
        {
          incorrectly_assigned_fit_param_.sigma = sqrt(sum_negative_sigma / sum_posterior);
          incorrectly_assigned_fit_param_.A = 1 / sqrt(2 * Constants::PI * pow(incorrectly_assigned_fit_param_.sigma, 2));
        }


        //compute new prior probabilities negative peptides
        fillDensities(x_scores, incorrect_density, correct_density);
        sum_posterior = sum_post(incorrect_density, correct_density);
        negative_prior_ = sum_posterior / x_scores.size();

        double new_maxlike(computeMaxLikelihood(incorrect_density, correct_density));
        if (boost::math::isnan(new_maxlike - maxlike) || new_maxlike < maxlike)
        {
          return false;
          //throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-PosteriorErrorProbability","Could not fit mixture model to data");
        }
        if ((new_maxlike - maxlike) < pow(10.0, -delta) || itns >= max_itns)
        {
          if (itns >= max_itns)
          {
            LOG_WARN << "Number of iterations exceeded. Convergence criterion not met. Last likelihood increase: " << (new_maxlike - maxlike) << endl;
            LOG_WARN << "Algorithm returns probabilites for suboptimal fit. You might want to try raising the max. number of iterations and have a look at the distribution." << endl;
          }
          stop_em_init = true;
          sum_posterior = sum_post(incorrect_density, correct_density);
          negative_prior_ = sum_posterior / x_scores.size();

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
      }
      while (!stop_em_init);
      //-------------------------------------------------------------
      // Finished fitting
      //-------------------------------------------------------------
      //!!Workaround:
      if (param_.getValue("incorrectly_assigned") == "Gumbel")
      {
        calc_incorrect_ = &PosteriorErrorProbabilityModel::getGumbel;
      }
      max_incorrectly_ = ((this)->*(calc_incorrect_))(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_);
      max_correctly_ = ((this)->*(calc_correct_))(correctly_assigned_fit_param_.x0, correctly_assigned_fit_param_);
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
      bool return_value;
      return_value = fit(search_engine_scores);
      if (!return_value)
        return false;

      probabilities.resize(search_engine_scores.size());
      vector<double>::iterator probs = probabilities.begin();
      for (vector<double>::iterator scores = search_engine_scores.begin(); scores != search_engine_scores.end(); ++scores, ++probs)
      {
        *probs = computeProbability(*scores);
      }
      return true;
    }

    void PosteriorErrorProbabilityModel::fillDensities(vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      if (incorrect_density.size() != x_scores.size())
      {
        incorrect_density.resize(x_scores.size());
        correct_density.resize(x_scores.size());
      }
      vector<double>::iterator incorrect = incorrect_density.begin();
      vector<double>::iterator correct = correct_density.begin();
      for (vector<double>::iterator scores = x_scores.begin(); scores != x_scores.end(); ++scores, ++incorrect, ++correct)
      {
        *incorrect = ((this)->*(calc_incorrect_))(*scores, incorrectly_assigned_fit_param_);
        *correct = ((this)->*(calc_correct_))(*scores, correctly_assigned_fit_param_);
      }
    }

    double PosteriorErrorProbabilityModel::computeMaxLikelihood(vector<double>& incorrect_density, vector<double>& correct_density)
    {
      double maxlike(0);
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect)
      {
        maxlike += log10(negative_prior_ * (*incorrect) + (1 - negative_prior_) * (*correct));
      }
      return maxlike;
    }

    double PosteriorErrorProbabilityModel::one_minus_sum_post(vector<double>& incorrect_density, vector<double>& correct_density)
    {
      double one_min(0);
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect)
      {
        one_min +=  1  - ((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)));
      }
      return one_min;
    }

    double PosteriorErrorProbabilityModel::sum_post(vector<double>& incorrect_density, vector<double>& correct_density)
    {
      double post(0);
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect)
      {
        post += ((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)));
      }
      return post;
    }

    double PosteriorErrorProbabilityModel::sum_pos_x0(vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      double pos_x0(0);
      vector<double>::iterator the_x = x_scores.begin();
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect, ++the_x)
      {
        pos_x0 += ((1  - ((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)))) * (*the_x));
      }
      return pos_x0;
    }

    double PosteriorErrorProbabilityModel::sum_neg_x0(vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density)
    {
      double neg_x0(0);
      vector<double>::iterator the_x = x_scores.begin();
      vector<double>::iterator correct = correct_density.begin();
      for (vector<double>::iterator incorrect = incorrect_density.begin(); incorrect < incorrect_density.end(); ++correct, ++incorrect, ++the_x)
      {
        neg_x0 += ((((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)))) * (*the_x));
      }
      return neg_x0;
    }

    double PosteriorErrorProbabilityModel::sum_pos_sigma(vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density, double positive_mean)
    {
      double pos_sigma(0);
      vector<double>::iterator the_x = x_scores.begin();
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect, ++the_x)
      {
        pos_sigma += ((1  - ((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)))) * pow((*the_x) - positive_mean, 2));
      }
      return pos_sigma;
    }

    double PosteriorErrorProbabilityModel::sum_neg_sigma(vector<double>& x_scores, vector<double>& incorrect_density, vector<double>& correct_density, double positive_mean)
    {
      double neg_sigma(0);
      vector<double>::iterator the_x = x_scores.begin();
      vector<double>::iterator incorrect = incorrect_density.begin();
      for (vector<double>::iterator correct = correct_density.begin(); correct < correct_density.end(); ++correct, ++incorrect, ++the_x)
      {
        neg_sigma += ((((negative_prior_ * (*incorrect)) / ((negative_prior_ * (*incorrect)) + (1 - negative_prior_) * (*correct)))) * pow((*the_x) - positive_mean, 2));
      }
      return neg_sigma;
    }

    double PosteriorErrorProbabilityModel::computeProbability(double score)
    {
      score = score + fabs(smallest_score_) + 0.001;
      double x_neg;
      double x_pos;
      // the score is smaller than the peak of incorrectly assigned sequences. To ensure that the probabilities wont rise again use the incorrectly assigned peak for computation
      if (score < incorrectly_assigned_fit_param_.x0)
      {
        x_neg = max_incorrectly_;
        x_pos = ((this)->*(calc_correct_))(score, correctly_assigned_fit_param_);
      }
      // same as above. However, this time to ensure that probabilities wont drop again.
      else if (score > correctly_assigned_fit_param_.x0)
      {
        x_neg = ((this)->*(calc_incorrect_))(score, incorrectly_assigned_fit_param_);
        x_pos = max_correctly_;
      }
      // if its in between use the normal formula
      else
      {
        x_neg = ((this)->*(calc_incorrect_))(score, incorrectly_assigned_fit_param_);
        x_pos = ((this)->*(calc_correct_))(score, correctly_assigned_fit_param_);
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
        LOG_WARN << "Target-Decoy plot was called, but '" << ListUtils::concatenate(empty, "' and '") << "' has no data! Unable to create a target-decoy plot." << std::endl;
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
      LOG_INFO << "Attempting to call 'gnuplot' ...";
      String cmd = String("gnuplot \"") + gp_file + "\"";
      if (system(cmd.c_str()))  // 0 is success!
      {
        LOG_WARN << "Calling 'gnuplot' on '" << gp_file << "' failed. Please create plots manually." << std::endl;
      }
      else LOG_INFO << " success!" << std::endl;

    }

  } // namespace Math
} // namespace OpenMS
