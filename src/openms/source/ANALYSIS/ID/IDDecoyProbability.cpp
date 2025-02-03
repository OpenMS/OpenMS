// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>

#include <boost/math/special_functions/gamma.hpp>
#include <fstream>

// #define IDDECOYPROBABILITY_DEBUG
// #undef  IDDECOYPROBABILITY_DEBUG

using namespace std;

namespace OpenMS
{
  IDDecoyProbability::IDDecoyProbability() :
    DefaultParamHandler("IDDecoyProbability")
  {
    defaults_.setValue("number_of_bins", 40, "Number of bins used for the fitting, if sparse datasets are used, this number should be smaller", {"advanced"});
    defaults_.setValue("lower_score_better_default_value_if_zero", 50.0, "This value is used if e.g. a E-value score is 0 and cannot be transformed in a real number (log of E-value)", {"advanced"});

#ifdef IDDECOYPROBABILITY_DEBUG
    defaults_.setValue("rev_filename", "", "bla", {"advanced"});
    defaults_.setValue("fwd_filename", "", "bla", {"advanced"});
#endif

    defaultsToParam_();
  }

  IDDecoyProbability::IDDecoyProbability(const IDDecoyProbability & rhs) = default;

  IDDecoyProbability::~IDDecoyProbability() = default;

  void IDDecoyProbability::apply(vector<PeptideIdentification> & ids)
  {
    double lower_score_better_default_value_if_zero(static_cast<double>(param_.getValue("lower_score_better_default_value_if_zero")));
    double lower_score_better_default_value_if_zero_exp = pow(10.0, -lower_score_better_default_value_if_zero);
    vector<double> rev_scores, fwd_scores, all_scores;

    // get the forward scores
    for (PeptideIdentification& pep : ids)
    {
      String score_type = pep.getScoreType();
      if (!pep.getHits().empty())
      {
        vector<PeptideHit> hits = pep.getHits();
        for (PeptideHit& pit : hits)
        {
          double score = pit.getScore();

          pit.setMetaValue(score_type + "_Score", score);

          if (!pep.isHigherScoreBetter())
          {
            if (score < lower_score_better_default_value_if_zero_exp)
            {
              score = lower_score_better_default_value_if_zero;
            }
            else
            {
              score = -log10(score);
            }
          }

          String target_decoy = (String)pit.getMetaValue("target_decoy");
          if (target_decoy == "target")
          {
            fwd_scores.push_back(score);
          }
          else if (target_decoy == "decoy")
          {
            rev_scores.push_back(score);
          }
          all_scores.push_back(score);
        }
        pep.setHits(hits);
      }
    }
#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << ids.size() << " " << rev_scores.size() << " " << fwd_scores.size() << " " << all_scores.size() << endl;
#endif
    apply_(ids, rev_scores, fwd_scores, all_scores);
    return;
  }

  void IDDecoyProbability::apply(vector<PeptideIdentification> & prob_ids, const vector<PeptideIdentification> & orig_fwd_ids, const vector<PeptideIdentification> & rev_ids)
  {
    double lower_score_better_default_value_if_zero((double)param_.getValue("lower_score_better_default_value_if_zero"));
    double lower_score_better_default_value_if_zero_exp = pow((double)10.0, -lower_score_better_default_value_if_zero);
    vector<PeptideIdentification> fwd_ids = orig_fwd_ids;
    vector<double> rev_scores, fwd_scores, all_scores;

    // get the forward scores
    for (PeptideIdentification& pep : fwd_ids)
    {
      String score_type = pep.getScoreType();
      if (!pep.getHits().empty())
      {
        vector<PeptideHit> hits = pep.getHits();
        for (PeptideHit& pit : hits)
        {
          double score = pit.getScore();

          pit.setMetaValue(score_type + "_Score", score);

          if (!pep.isHigherScoreBetter())
          {
            if (score < lower_score_better_default_value_if_zero_exp)
            {
              score = lower_score_better_default_value_if_zero;
            }
            else
            {
              score = -log10(score);
            }
          }
          fwd_scores.push_back(score);
          all_scores.push_back(score);
        }
        pep.setHits(hits);
      }
    }

    // get the reverse scores
    for (const PeptideIdentification& pep : rev_ids)
    {
      if (!pep.getHits().empty())
      {
        for (const PeptideHit& pit : pep.getHits())
        {
          double score = pit.getScore();
          if (!pep.isHigherScoreBetter())
          {
            if (score < lower_score_better_default_value_if_zero_exp)
            {
              score = lower_score_better_default_value_if_zero;
            }
            else
            {
              score = -log10(score);
            }
          }

          rev_scores.push_back(score);
          all_scores.push_back(score);
        }
      }
    }

    prob_ids = fwd_ids;
    apply_(prob_ids, rev_scores, fwd_scores, all_scores);
    return;
  }

  void IDDecoyProbability::apply_(vector<PeptideIdentification> & ids, const vector<double> & rev_scores, const vector<double> & fwd_scores, const vector<double> & all_scores)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));



    // normalize distribution to [0, 1]
    vector<double> fwd_scores_normalized(number_of_bins, 0.0), rev_scores_normalized(number_of_bins, 0.0), diff_scores(number_of_bins, 0.0), all_scores_normalized(number_of_bins, 0.0);
    Transformation_ rev_trafo, fwd_trafo, all_trafo;
    normalizeBins_(rev_scores, rev_scores_normalized, rev_trafo);
    normalizeBins_(fwd_scores, fwd_scores_normalized, fwd_trafo);
    normalizeBins_(all_scores, all_scores_normalized, all_trafo);

    // rev scores fitting
    vector<DPosition<2> > rev_data;

    for (Size i = 0; i < number_of_bins; ++i)
    {
      DPosition<2> pos;
      pos.setX(((double)i) / (double)number_of_bins + 0.0001);    // necessary????
      pos.setY(rev_scores_normalized[i]);
      rev_data.push_back(pos);
#ifdef IDDECOYPROBABILITY_DEBUG
      cerr << pos.getX() << " " << pos.getY() << endl;
#endif
    }

    Math::GammaDistributionFitter gdf;
    Math::GammaDistributionFitter::GammaDistributionFitResult result_gamma_1st (1.0, 3.0);
    gdf.setInitialParameters(result_gamma_1st);
    // TODO heuristic for good start parameters
    Math::GammaDistributionFitter::GammaDistributionFitResult result_gamma = gdf.fit(rev_data);

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << gdf.getGnuplotFormula() << endl;
    String rev_filename = param_.getValue("rev_filename");
    generateDistributionImage_(rev_scores_normalized, gdf.getGnuplotFormula(), rev_filename);
#endif

    // generate diffs of distributions
    // get the fwd and rev distribution, apply all_trafo and calculate the diff
    vector<Size> fwd_bins(number_of_bins, 0), rev_bins(number_of_bins, 0);
    double min(all_trafo.min_score), diff(all_trafo.diff_score);
    Size max_bin(0);
    for (vector<double>::const_iterator it = fwd_scores.begin(); it != fwd_scores.end(); ++it)
    {
      Size bin = (Size)((*it - min) / diff * (double)(number_of_bins - 1));
      ++fwd_bins[bin];
      if (fwd_bins[bin] > max_bin)
      {
        max_bin = fwd_bins[bin];
      }
    }

    Size max_reverse_bin(0), max_reverse_bin_value(0);
    //min = rev_trafo.min_score;
    //diff = rev_trafo.diff_score;
    for (vector<double>::const_iterator it = rev_scores.begin(); it != rev_scores.end(); ++it)
    {
      Size bin = (Size)((*it - min) / diff * (double)number_of_bins);
      ++rev_bins[bin];
      if (rev_bins[bin] > max_bin)
      {
        max_bin = rev_bins[bin];
      }
      if (rev_bins[bin] > max_reverse_bin_value)
      {
        max_reverse_bin = bin;
        max_reverse_bin_value = rev_bins[bin];
      }
    }

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "Trying to get diff scores" << endl;
#endif

    // get diff of fwd and rev
    for (Size i = 0; i < number_of_bins; ++i)
    {
      Size fwd(0), rev(0);
      fwd = fwd_bins[i];
      rev = rev_bins[i];
      if ((double)fwd > (double)(1.3 * rev) && max_reverse_bin < i)
      {
        diff_scores[i] = (double)(fwd - rev) / (double)max_bin;
      }
      else
      {
        diff_scores[i] = 0.0;
      }
    }
#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "Gauss Fitting values size of diff scores=" << diff_scores.size() << endl;
#endif
    // diff scores fitting
    vector<DPosition<2> > diff_data;
    double gauss_A(0), gauss_x0(0), norm_factor(0);
    for (Size i = 0; i < number_of_bins; ++i)
    {
      DPosition<2> pos;
      pos.setX((double)i / (double)number_of_bins);
      pos.setY(diff_scores[i]);

      if (pos.getY() > gauss_A)
      {
        gauss_A = pos.getY();
      }
      gauss_x0 += pos.getX() * pos.getY();
      norm_factor += pos.getY();


      diff_data.push_back(pos);
    }

    double gauss_sigma(0);
    gauss_x0 /= (double)diff_data.size();
    gauss_x0 /= norm_factor;

    for (Size i = 0; i <= number_of_bins; ++i)
    {
      gauss_sigma += fabs(gauss_x0 - (double)i / (double)number_of_bins);
    }

    gauss_sigma /= (double)diff_data.size();



#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "setting initial parameters: " << endl;
#endif
    Math::GaussFitter gf;
    Math::GaussFitter::GaussFitResult result_1st(gauss_A, gauss_x0, gauss_sigma);
    gf.setInitialParameters(result_1st);
#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "Initial Gauss guess: A=" << gauss_A << ", x0=" << gauss_x0 << ", sigma=" << gauss_sigma << endl;
#endif

    //TODO: fail-to-fit correction was done using the GNUPlotFormula. Seemed to be a hack.
    //Changed it to try-catch-block but I am not sure if this correction should be made
    //at all. Can someone please verify?
    Math::GaussFitter::GaussFitResult result_gauss (gauss_A, gauss_x0, gauss_sigma);
    try {
        result_gauss = gf.fit(diff_data);
    }
    catch(Exception::UnableToFit& /* e */)
    {
      result_gauss.A = gauss_A;
      result_gauss.x0 = gauss_x0;
      result_gauss.sigma = gauss_sigma;
    }

//    // fit failed?
//    if (gf.getGnuplotFormula() == "")
//    {
//      result_gauss.A = gauss_A;
//      result_gauss.x0 = gauss_x0;
//      result_gauss.sigma = gauss_sigma;
//    }

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << gf.getGnuplotFormula() << endl;
    String fwd_filename = param_.getValue("fwd_filename");
    if (gf.getGnuplotFormula() == "")
    {
      String formula("f(x)=" + String(gauss_A) + " * exp(-(x - " + String(gauss_x0) + ") ** 2 / 2 / (" + String(gauss_sigma) + ") ** 2)");
      generateDistributionImage_(diff_scores, formula, fwd_filename);
    }
    else
    {
      generateDistributionImage_(diff_scores, gf.getGnuplotFormula(), fwd_filename);
    }
#endif

#ifdef IDDECOYPROBABILITY_DEBUG
    //all_trafo.diff_score + all_trafo.min_score
    String gauss_formula("f(x)=" + String(result_gauss.A / all_trafo.max_intensity) + " * exp(-(x - " + String(result_gauss.x0 * all_trafo.diff_score + all_trafo.min_score) + ") ** 2 / 2 / (" + String(result_gauss.sigma * all_trafo.diff_score)   + ") ** 2)");

    String b_str(result_gamma.b), p_str(result_gamma.p);
    String gamma_formula = "g(x)=(" + b_str + " ** " + p_str + ") / gamma(" + p_str + ") * x ** (" + p_str + " - 1) * exp(- " + b_str + " * x)";

    generateDistributionImage_(all_scores_normalized, all_trafo, gauss_formula, gamma_formula, (String)param_.getValue("fwd_filename"));
#endif

    vector<PeptideIdentification> new_prob_ids;
    // calculate the probabilities and write them to the IDs
    for (const PeptideIdentification& pep : ids)
    {
      if (!pep.getHits().empty())
      {
        vector<PeptideHit> hits;
        String score_type = pep.getScoreType() + "_score";
        for (const PeptideHit& pit : pep.getHits())
        {
          PeptideHit hit = pit;
          double score = hit.getScore();
          if (!pep.isHigherScoreBetter())
          {
            score = -log10(score);
          }
          hit.setMetaValue(score_type, hit.getScore());
          hit.setScore(getProbability_(result_gamma, rev_trafo, result_gauss, fwd_trafo, score));
          hits.push_back(hit);
        }
        PeptideIdentification id = pep;
        id.setHigherScoreBetter(true);
        id.setScoreType(id.getScoreType() + "_DecoyProbability");
        id.setHits(hits);

        new_prob_ids.push_back(id);
      }
    }
    ids = new_prob_ids;
  }

  // normalize the bins to [0, 1]
  void IDDecoyProbability::normalizeBins_(const vector<double> & scores, vector<double> & binned, Transformation_ & trafo)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));
    // get the range of the scores
    double max(numeric_limits<double>::min()), min(numeric_limits<double>::max());
    for (vector<double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      if (*it > max)
      {
        max = *it;
      }
      if (*it < min)
      {
        min = *it;
      }
    }

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "Range is [" << min << ", " << max << "]" << endl;
#endif

    // perform the binning
    double diff = max - min;
    Size max_bin_number(0);
    double max_bin(0);
    for (vector<double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      Size bin = (Size)((*it - min) / diff * (double)(number_of_bins - 1));
      binned[bin] += 1;

      if (binned[bin] > max_bin)
      {
        max_bin = binned[bin];
        max_bin_number = bin;
      }
    }


    // normalize to \sum = 1
    for (vector<double>::iterator it = binned.begin(); it != binned.end(); ++it)
    {
      *it /= (double)max_bin / 4.0;   // 4 is best value for the gamma distribution
    }


    // store the transformation
    trafo.max_intensity = 4.0 / (double)max_bin;
    trafo.diff_score = diff;
    trafo.min_score = min;
    trafo.max_intensity_bin = max_bin_number;
    trafo.max_score = max;

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "TRAFO: max_intensity=" << trafo.max_intensity << ", diff_score=" << trafo.diff_score << ", min_score=" << trafo.min_score << ", max_intensity_bin=" << trafo.max_intensity_bin << ", max_score=" << trafo.max_score << endl;
#endif
  }

  double IDDecoyProbability::getProbability_(const Math::GammaDistributionFitter::GammaDistributionFitResult & result_gamma,
                                                 const Transformation_ & gamma_trafo,
                                                 const Math::GaussFitter::GaussFitResult & result_gauss,
                                                 const Transformation_ & gauss_trafo,
                                                 double score)
  {
    double rho_rev(0), rho_fwd(0);
    Size number_of_bins(param_.getValue("number_of_bins"));

    // first transform the score into a background distribution density value
    double score_rev_trans = (score - gamma_trafo.min_score) / gamma_trafo.diff_score;
    if (score_rev_trans < gamma_trafo.max_intensity_bin / (double)number_of_bins)
    {
      rho_rev = 1.0 / gamma_trafo.max_intensity;
    }
    else
    {
      rho_rev = pow(result_gamma.b, result_gamma.p) / std::tgamma(result_gamma.p) * pow(score_rev_trans, result_gamma.p - 1) * exp(-result_gamma.b * score_rev_trans);
    }

    // second transform the score into a 'correct' distribution density value
    double score_fwd_trans = (score - gauss_trafo.min_score) / gauss_trafo.diff_score;

#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "score=" << score << ", score_rev_trans=" << score_rev_trans << ", score_fwd_trans=" << score_fwd_trans << ", rho_rev=" << rho_rev << ", gauss_trafor.max_score=" << gauss_trafo.max_score;
#endif

    if (score_fwd_trans < result_gauss.x0)
    {
#ifdef IDDECOYPROBABILITY_DEBUG
      cerr << "(score_fwd_trans > gauss_trafo.max_score, " << score_fwd_trans << " " << gauss_trafo.max_score << " -> 1)" << endl;
#endif
      rho_fwd = result_gauss.A * exp(-pow(score_fwd_trans - result_gauss.x0, 2) / 2.0 / pow(result_gauss.sigma, 2));
    }
    else
    {
      rho_fwd = 1;
    }


#ifdef IDDECOYPROBABILITY_DEBUG
    cerr << "rho_fwd=" << rho_fwd << endl;
#endif

    // calc P using Bayes theorem
    return rho_fwd / (rho_fwd + rho_rev);
  }

  void IDDecoyProbability::generateDistributionImage_(const vector<double> & ids, const String & formula, const String & filename)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));

    // write distribution to file
    ofstream o((filename + "_dist_tmp.dat").c_str());
    for (Size i = 0; i < number_of_bins; ++i)
    {
      o << (double)i / (double)number_of_bins << " " << ids[i] << endl;
    }
    o.close();

    ofstream os((filename + "_gnuplot.gpl").c_str());
    os << "set terminal png" << endl;
    os << "set output '" << filename << "_distribution.png'" << endl;
    os << formula << endl;
    os << "plot f(x), '" << filename << "_dist_tmp.dat' w boxes" << endl;
    os.close();

#ifdef IDDECOYPROBABILITY_DEBUG
    Int syscalret = system(("gnuplot " + filename + "_gnuplot.gpl").c_str());
    if (syscalret != 0)
    {
      cerr << "gnuplot system call failed!" << endl;
    }
#endif


    return;
  }

  void IDDecoyProbability::generateDistributionImage_(const vector<double> & all_ids, const Transformation_ & all_trans, const String & fwd_formula, const String & rev_formula, const String & filename)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));

    ofstream all_output((filename + "_all_tmp.dat").c_str());
    for (Size i = 0; i < number_of_bins; ++i)
    {
      all_output << (double)i / (double)number_of_bins * all_trans.diff_score + all_trans.min_score << " " << all_ids[i] / all_trans.max_intensity << endl;
    }
    all_output.close();

    ofstream os((filename + "_both_gnuplot.gpl").c_str());
    os << "set terminal png" << endl;
    os << "set output '" << filename << "_both_distributions.png'" << endl;
    os << fwd_formula << endl;
    os << rev_formula << endl;
    //os << "plot f(x), '" << filename << "_fwd_tmp.dat' w boxes, g(x), '" << filename << "_rev_tmp.dat' w boxes, '" << filename << "_all_tmp.dat' w i" << endl;
    os << "plot f(x), g(x), '" << filename << "_all_tmp.dat' w i" << endl;
    os.close();

#ifdef IDDECOYPROBABILITY_DEBUG
    Int syscalret = system(("gnuplot " + filename + "_both_gnuplot.gpl").c_str());
    if (syscalret != 0)
    {
      cerr << "gnuplot system call failed!" << endl;
    }
#endif

    return;
  }

} // namespace OpenMS
