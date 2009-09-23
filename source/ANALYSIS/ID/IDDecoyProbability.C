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
// $Maintainer: Andreas Bertsch, Sven Nahnsen $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <boost/math/special_functions/gamma.hpp>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>

#include <fstream>

#define IDDECOYPROBABILITY_DEBUG
#undef  IDDECOYPROBABILITY_DEBUG

using namespace std;

namespace OpenMS 
{
	IDDecoyProbability::IDDecoyProbability()
		: DefaultParamHandler("IDDecoyProbability")
  {
		defaults_.setValue("number_of_bins", 40, "Number of bins used for the fitting, if sparse datasets are used, this number should be smaller", StringList::create("advanced"));
		defaults_.setValue("lower_score_better_default_value_if_zero", 50.0, "This value is used if e.g. a E-value score is 0 and cannot be transformed in a real number (log of E-value)", StringList::create("advanced"));

#ifdef IDDECOYPROBABILITY_DEBUG
		defaults_.setValue("rev_filename", "", "bla", StringList::create("advanced"));
		defaults_.setValue("fwd_filename", "", "bla", StringList::create("advanced"));
#endif

		defaultsToParam_();
  }

	IDDecoyProbability::IDDecoyProbability(const IDDecoyProbability& rhs)
		: DefaultParamHandler(rhs)
	{

	}

	IDDecoyProbability::~IDDecoyProbability()
	{
	}

	void IDDecoyProbability::apply(vector<PeptideIdentification>& prob_ids, const vector<PeptideIdentification>& orig_fwd_ids, const vector<PeptideIdentification>& rev_ids)
	{
		Size number_of_bins(param_.getValue("number_of_bins"));
		double lower_score_better_default_value_if_zero((double)param_.getValue("lower_score_better_default_value_if_zero"));

		vector<PeptideIdentification> fwd_ids = orig_fwd_ids;
  	vector<double> rev_scores, fwd_scores, all_scores;

		// get the forward scores
  	for (vector<PeptideIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
  	{
			String score_type = it->getScoreType();
    	if (it->getHits().size() > 0)
    	{
				vector<PeptideHit> hits = it->getHits();
     		for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
     		{
       		double score = pit->getScore();

					pit->setMetaValue(score_type+"_Score", score);
					
       		if (!it->isHigherScoreBetter())
       		{
						if (score == 0)
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
				it->setHits(hits);
   		}
 		}

		// get the reverse scores
		for (vector<PeptideIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
  	{
    	if (it->getHits().size() > 0)
    	{
      	for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      	{
        	double score = pit->getScore();
        	if (!it->isHigherScoreBetter())
        	{
						if (score == 0)
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
    	pos.setX(((double)i) / (double)number_of_bins + 0.0001);  // necessary????
    	pos.setY(rev_scores_normalized[i]);
    	rev_data.push_back(pos);

			cerr << pos.getX() << " " << pos.getY() << endl;
  	}

		Math::GammaDistributionFitter gdf;
		Math::GammaDistributionFitter::GammaDistributionFitResult result_gamma_1st;
		result_gamma_1st.b = 1.0;
		result_gamma_1st.p = 3.0;
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
    	if (fwd > rev && max_reverse_bin < i)
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
		double gauss_A(0), gauss_x0(0);
  	for (Size i = 0; i < number_of_bins; ++i)
  	{
    	DPosition<2> pos;
    	pos.setX((double)i / (double)number_of_bins);
    	pos.setY(diff_scores[i]);

			if (pos.getY() > gauss_A)
			{
				gauss_A = pos.getY();
			}
			gauss_x0 += pos.getX();
			
    	diff_data.push_back(pos);
  	}

		double gauss_sigma(0);
		gauss_x0 /= (double)diff_data.size();
		
		for (Size i = 0; i <= number_of_bins; ++i)
		{
			gauss_sigma += fabs(gauss_x0 - (double)i / (double)number_of_bins);
		}

		gauss_sigma /= (double)diff_data.size();

	
		
#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "setting initial parameters: " << endl;
#endif
  	Math::GaussFitter gf;
  	Math::GaussFitter::GaussFitResult result_1st;
  	result_1st.A = gauss_A; //0.06;
  	result_1st.x0 = gauss_x0; //0.7;
  	result_1st.sigma = gauss_sigma; //0.5;
		gf.setInitialParameters(result_1st);
#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "Initial Gauss guess: A=" << gauss_A << ", x0=" << gauss_x0 << ", sigma=" << gauss_sigma << endl;
#endif
  	
		Math::GaussFitter::GaussFitResult result_gauss = gf.fit(diff_data);

		// fit failed?
		if (gf.getGnuplotFormula() == "")
		{
			result_gauss.A = gauss_A;
			result_gauss.x0 = gauss_x0;
			result_gauss.sigma = gauss_sigma;
		}

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

		// calculate the probabilities and write them to the IDs
		for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      if (it->getHits().size() > 0)
      {
				vector<PeptideHit> hits;
				String score_type = it->getScoreType() + "_score";
				for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
        {
					PeptideHit hit = *pit;
					double score = hit.getScore();
					if (!it->isHigherScoreBetter())
					{
						score = -log10(score);
					}
					hit.setMetaValue(score_type, hit.getScore());
					hit.setScore(getProbability_(result_gamma, rev_trafo, result_gauss, fwd_trafo, score));
					hits.push_back(hit);
				}
				PeptideIdentification id = *it;
				id.setHigherScoreBetter(true);
				id.setScoreType(id.getScoreType() + "_DecoyProbability");
				id.setHits(hits);
				
				prob_ids.push_back(id);
			}
		}
	}
	
	// normalize the bins to [0, 1]
	void IDDecoyProbability::normalizeBins_(const vector<double>& scores, vector<double>& binned, Transformation_& trafo)
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
    	*it /= (double)max_bin / 4.0; // 4 is best value for the gamma distribution
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

	double IDDecoyProbability::getProbability_(const Math::GammaDistributionFitter::GammaDistributionFitResult& result_gamma,
																						const Transformation_& gamma_trafo,
																						const Math::GaussFitter::GaussFitResult& result_gauss,
																						const Transformation_& gauss_trafo,
																						double score)
	{
  	double rho_rev(0), rho_fwd(0);
		Size number_of_bins(param_.getValue("number_of_bins"));
	
  	// first transform the score into a background distribution density value
  	double score_rev_trans = (score - gamma_trafo.min_score) / gamma_trafo.diff_score;
  	if (score_rev_trans < gamma_trafo.max_intensity_bin/(double)number_of_bins)
  	{
    	rho_rev = 1.0 / gamma_trafo.max_intensity;
  	}
  	else
  	{
			rho_rev = pow(result_gamma.b, result_gamma.p) / boost::math::tgamma(result_gamma.p) * pow(score_rev_trans, result_gamma.p - 1) * exp(- result_gamma.b * score_rev_trans);
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
			rho_fwd = result_gauss.A * exp(- pow(score_fwd_trans - result_gauss.x0, 2)/2.0/pow(result_gauss.sigma, 2));
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

	void IDDecoyProbability::generateDistributionImage_(const vector<double>& ids, const String& formula, const String& filename)
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

  	system(("gnuplot " + filename + "_gnuplot.gpl").c_str());

  	return;
	}	

	void IDDecoyProbability::generateDistributionImage_(const vector<double>& all_ids, const Transformation_& all_trans, const String& fwd_formula, const String& rev_formula, const String& filename)
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

    system(("gnuplot " + filename + "_both_gnuplot.gpl").c_str());

    return;
	}
} // namespace OpenMS
