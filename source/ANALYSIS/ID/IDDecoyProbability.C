// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

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
		defaults_.setValue("number_of_bins", 40.0, "Number of bins used for the fitting, if sparse datasets are used, this number should be smaller", StringList::create("advanced"));
		defaults_.setValue("lower_score_better_default_value_if_zero", 50.0, "This value is used if e.g. a E-value score is 0 and cannot be transformed in a real number (log of E-value)", StringList::create("advanced"));

#ifdef IDDECOYPROBABILITY_DEBUG
		defaults_.setValue("rev_filename", "", StringList::create("advanced"));
		defaults_.setValue("fwd_filename", "", StringList::create("advanced"));
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
		double number_of_bins((double)param_.getValue("number_of_bins"));
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
  	Map<double, double> fwd_scores_normalized, rev_scores_normalized, diff_scores, all_scores_normalized;
  	Transformation_ rev_trafo, fwd_trafo, all_trafo;
#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "rev-bins" << endl;
#endif
  	normalizeBins_(rev_scores, rev_scores_normalized, rev_trafo);
#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "fwd-bins" << endl;
#endif
  	normalizeBins_(fwd_scores, fwd_scores_normalized, fwd_trafo);
#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "all-bins" << endl;
#endif
  	normalizeBins_(all_scores, all_scores_normalized, all_trafo);

  	// rev scores fitting
  	vector<DPosition<2> > rev_data;

  	for (UInt i = 0; i != rev_scores_normalized.size(); ++i)
  	{
    	DPosition<2> pos;
    	pos.setX(((double)i + 1.0) / number_of_bins);
    	pos.setY(rev_scores_normalized[(double)i / number_of_bins]);
    	rev_data.push_back(pos);
  	}


		Math::GammaDistributionFitter gdf;
  	// TODO heuristic for good start parameters
  	Math::GammaDistributionFitter::GammaDistributionFitResult result_gamma = gdf.fit(rev_data);

#ifdef IDDECOYPROBABILITY_DEBUG
  	cerr << gdf.getGnuplotFormula() << endl;
		String rev_filename = param_.getValue("rev_filename");
		generateDistributionImage(rev_scores_normalized, gdf.getGnuplotFormula(), rev_filename);
		cerr << "blubb" << endl;
#endif
		
  	// generate diffs of distributions
  	// get the fwd and rev distribution, apply all_trafo and calculate the diff
  	Map<UInt, UInt> fwd_bins, rev_bins;
  	double min(rev_trafo.x_shift);
  	double diff(rev_trafo.x_factor);
  	UInt max_bin(0);
		UInt max_reverse_bin(0);
		UInt max_reverse_bin_value(0);
  	for (vector<double>::const_iterator it = fwd_scores.begin(); it != fwd_scores.end(); ++it)
  	{
    	UInt bin = (UInt)((*it - min) / diff * number_of_bins);
    	if (fwd_bins.has(bin))
    	{
      	++fwd_bins[bin];
    	}
    	else
    	{
      	fwd_bins[bin] = 1;
    	}
    	if (fwd_bins[bin] > max_bin)
    	{
      	max_bin = fwd_bins[bin];
    	}
  	}
		
  	for (vector<double>::const_iterator it = rev_scores.begin(); it != rev_scores.end(); ++it)
  	{
    	UInt bin = (UInt)((*it - min) / diff * number_of_bins);
    	if (rev_bins.has(bin))
    	{
      	++rev_bins[bin];
    	}
    	else
    	{
      	rev_bins[bin] = 1;
    	}
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
  	for (UInt i = 0; i < number_of_bins; ++i)
  	{
    	UInt fwd(fwd_bins[i]), rev(rev_bins[i]);
    	if (fwd > rev && max_reverse_bin < i)
    	{
      	diff_scores[(double)i / number_of_bins] = (double)(fwd - rev) / (double)max_bin * 4.0;
    	}
    	else
    	{
      	diff_scores[(double)i / number_of_bins] = 0.0;
    	}
  	}

#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "Gauss Fitting values size of diff scores=" << diff_scores.size() << endl;
#endif
  	// diff scores fitting
  	vector<DPosition<2> > diff_data;
		double gauss_A(0);
		double gauss_x0(0);
  	for (UInt i = 0; i != diff_scores.size(); ++i)
  	{
    	DPosition<2> pos;
    	pos.setX(((double)i + 1.0) / (double)number_of_bins);
    	pos.setY(diff_scores[(double)i / (double)number_of_bins]);

			if (pos.getY() > gauss_A)
			{
				gauss_A = pos.getY();
			}
			gauss_x0 += pos.getX();
			
#ifdef IDDECOYPROBABILITY_DEBUG
			cerr << pos.getX() << " " << pos.getY() << endl;
#endif
			
    	diff_data.push_back(pos);
  	}

		double gauss_sigma(0);
		gauss_x0 /= (double)diff_scores.size();
		
		for (UInt i = 0; i != diff_scores.size(); ++i)
		{
			gauss_sigma += fabs(gauss_x0 - ((double)i + 1.0) / (double)number_of_bins);
		}

		gauss_sigma /= (double)diff_scores.size();

#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "setting initial parameters: " << endl;
#endif
  	Math::GaussFitter gf;
  	Math::GaussFitter::GaussFitResult result_1st;
  	result_1st.A = gauss_A; //0.06;
  	result_1st.x0 = gauss_x0; //0.7;
  	result_1st.sigma = gauss_sigma; //0.5;
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
			generateDistributionImage(diff_scores, formula, fwd_filename); 
		}
		else
		{
			generateDistributionImage(diff_scores, gf.getGnuplotFormula(), fwd_filename);
		}
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
	void IDDecoyProbability::normalizeBins_(const vector<double>& scores, Map<double, double>& binned, Transformation_& trafo)
	{
		double number_of_bins((double)param_.getValue("number_of_bins"));
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
  	Map<UInt, UInt> bins;
  	UInt max_bin(0), max_bin_number(0);
  	for (vector<double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
  	{
    	UInt bin = (UInt)((*it - min) / diff * number_of_bins);
    	if (bins.has(bin))
    	{
      	++bins[bin];
    	}
    	else
    	{
      	bins[bin] = 1;
    	}
    	if (bins[bin] > max_bin)
    	{
      	max_bin = bins[bin];
      	max_bin_number = bin;
    	}
  	}

  	// normalize to \sum = 1
  	for (Map<UInt, UInt>::ConstIterator it = bins.begin(); it != bins.end(); ++it)
  	{
    	binned[it->first/number_of_bins] = it->second / (double)max_bin * 4.0; // 4 is best value for the gamma distribution
#ifdef IDDECOYPROBABILITY_DEBUG
			cerr << it->first/number_of_bins * diff + min << " " << it->second << endl;
#endif
  	}

		// store the transformation
  	trafo.y_factor = 4.0 / (double)max_bin;
  	trafo.x_factor = diff;
  	trafo.x_shift = min;
  	trafo.y_max_bin = max_bin_number;
  	trafo.x_max = max;

#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "TRAFO: y_factor=" << trafo.y_factor << ", x_factor=" << trafo.x_factor << ", x_shift=" << trafo.x_shift << ", y_max_bin=" << trafo.y_max_bin << ", x_max=" << trafo.x_max << endl;
#endif
	}

	double IDDecoyProbability::getProbability_(const Math::GammaDistributionFitter::GammaDistributionFitResult& result_gamma,
																						const Transformation_& gamma_trafo,
																						const Math::GaussFitter::GaussFitResult& result_gauss,
																						const Transformation_& gauss_trafo,
																						double score)
	{
  	double rho_rev(0), rho_fwd(0);
		double number_of_bins((double)param_.getValue("number_of_bins"));
	
  	// first transform the score into a background distribution density value
  	double score_rev_trans = (score - gamma_trafo.x_shift) / gamma_trafo.x_factor;
  	if (score_rev_trans < gamma_trafo.y_max_bin/number_of_bins)
  	{
    	rho_rev = 1.0 / gamma_trafo.y_factor;
  	}
  	else
  	{
    	rho_rev = pow(result_gamma.b, result_gamma.p) / tgamma(result_gamma.p) * pow(score_rev_trans, result_gamma.p - 1) * exp(- result_gamma.b * score_rev_trans);
  	}

  	// second transform the score into a 'correct' distribution density value
  	double score_fwd_trans = (score - gauss_trafo.x_shift) / gauss_trafo.x_factor;

#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "score=" << score << ", score_rev_trans=" << score_rev_trans << ", score_fwd_trans=" << score_fwd_trans << ", rho_rev=" << rho_rev << " ";
#endif
		
  	if (score_fwd_trans > gauss_trafo.x_max)
  	{
#ifdef IDDECOYPROBABILITY_DEBUG
			cerr << "(score_fwd_trans > gauss_trafo.x_max, " << score_fwd_trans << " " << gauss_trafo.x_max << " -> 1)" << endl;
#endif
    	return 1;
  	}
		
  	rho_fwd = result_gauss.A * exp(- pow(score_fwd_trans - result_gauss.x0, 2)/2.0/pow(result_gauss.sigma, 2));

#ifdef IDDECOYPROBABILITY_DEBUG
		cerr << "rho_fwd=" << rho_fwd << endl;
#endif
		
  	// calc P using Bayes theorem
  	return rho_fwd / (rho_fwd + rho_rev);
	}

	void IDDecoyProbability::generateDistributionImage(const Map<double, double>& ids, const String& formula, const String& filename)
	{
		double number_of_bins((double)param_.getValue("number_of_bins"));
		
  	// write distribution to file
  	ofstream o((filename + "_dist_tmp.dat").c_str());
  	for (UInt i = 0; i < number_of_bins; ++i)
  	{
			if (ids.has((double)i/number_of_bins))
			{
    		o << (double)i / number_of_bins << " " << ids[(double)i/number_of_bins] << endl;
			}
			else
			{
				cerr << "Error: some bins are not present!" << endl;
			}
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
} // namespace OpenMS
