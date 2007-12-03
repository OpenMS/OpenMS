// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

using namespace std;

namespace OpenMS 
{
	IDDecoyProbability::IDDecoyProbability()
		: DefaultParamHandler("IDDecoyProbability")
  {
    defaults_.setValue("bin_size", 100, true);
		defaults_.setValue("base", 20.0, true);
		defaults_.setValue("discretization", 100.0, true);
		defaults_.setValue("number_of_bins", 40.0, true);

		defaultsToParam_();
  }

	IDDecoyProbability::IDDecoyProbability(const IDDecoyProbability& rhs)
		: DefaultParamHandler(rhs)
	{

	}

	IDDecoyProbability::~IDDecoyProbability()
	{
	}

	void IDDecoyProbability::apply(vector<PeptideIdentification>& prob_ids, const vector<PeptideIdentification>& fwd_ids, const vector<PeptideIdentification>& rev_ids) throw (Exception::MissingInformation)
	{
		UInt bin_size((UInt)param_.getValue("bin_size"));
		double base((double)param_.getValue("base"));
		double discretization((double)param_.getValue("discretization"));
		double number_of_bins((double)param_.getValue("number_of_bins"));


		HashMap<UInt, UInt> fwd_scores_bins, rev_scores_bins;

  	for (UInt i = 0; i != bin_size; ++i)
  	{
    	fwd_scores_bins[i] = 0;
    	rev_scores_bins[i] = 0;
  	}

  	vector<double> rev_scores, fwd_scores, all_scores;

		// get the forward scores
  	for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
  	{
    	if (it->getHits().size() > 0)
    	{
     		for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
     		{
       		double score = pit->getScore();
       		if (!it->isHigherScoreBetter())
       		{
         		score = -log10(score);
       		}
       		fwd_scores.push_back(score);
       		all_scores.push_back(score);
       		UInt bin = min(bin_size - 1, (UInt)((score + base)/double(discretization) * (double)bin_size));
       		if (fwd_scores_bins.has(bin))
       		{
         		fwd_scores_bins[bin]++;
       		}
       		else
       		{
         		fwd_scores_bins[bin] = 1;
       		}
        	//break; // TODO
     		}
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
        	if (!it->isHigherScoreBetter() && score != 0) // TODO, what to do with score 0????
        	{
          	score = -log10(score);
        	}

        	rev_scores.push_back(score);
        	all_scores.push_back(score);

        	UInt bin = min(bin_size - 1, (UInt)((score + base)/double(discretization) * (double)bin_size));
        	if (rev_scores_bins.has(bin))
        	{
          	rev_scores_bins[bin]++;
        	}
        	else
        	{
          	rev_scores_bins[bin] = 1;
        	}
      	}
    	}
  	}

  	// normalize distribution to [0, 1]
  	HashMap<double, double> fwd_scores_normalized, rev_scores_normalized, diff_scores, all_scores_normalized;
  	Transformation_ rev_trafo, fwd_trafo, all_trafo;
  	normalizeBins_(rev_scores, rev_scores_normalized, rev_trafo);
  	normalizeBins_(fwd_scores, fwd_scores_normalized, fwd_trafo);
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


		GammaDistributionFitter gdf;
  	// TODO heuristic for good start parameters
  	GammaDistributionFitter::GammaDistributionFitResult result_gamma = gdf.fit(rev_data);

#ifdef IDDECOYPROBABILITY_DEBUG
  	cerr << gdf.getGnuplotFormula() << endl;
#endif
		
  	// generate diffs of distributions
  	// get the fwd and rev distribution, apply all_trafo and calculate the diff
  	HashMap<UInt, UInt> fwd_bins, rev_bins;
  	double min(rev_trafo.x_shift);
  	double diff(rev_trafo.x_factor);
  	UInt max_bin(0);
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
  	}

  	// get diff of fwd and rev
  	for (UInt i = 0; i < number_of_bins; ++i)
  	{
    	UInt fwd(fwd_bins[i]), rev(rev_bins[i]);
    	if (fwd > rev)
    	{
      	diff_scores[(double)i / number_of_bins] = (double)(fwd - rev) / (double)max_bin * 4.0;
    	}
    	else
    	{
      	diff_scores[(double)i / number_of_bins] = 0.0;
    	}
  	}

  	// diff scores fitting
  	vector<DPosition<2> > diff_data;
  	for (UInt i = 0; i != diff_scores.size(); ++i)
  	{
    	DPosition<2> pos;
    	pos.setX(((double)i + 1.0) / number_of_bins);
    	pos.setY(diff_scores[(double)i / number_of_bins]);
    	diff_data.push_back(pos);
  	}

  	GaussFitter gf;
  	GaussFitter::GaussFitResult result_1st = gf.getInitialParameters();
		// TODO find good start parameters
  	result_1st.A = 0.06;
  	result_1st.x0 = 0.7;
  	result_1st.sigma = 0.5;
  	GaussFitter::GaussFitResult result_gauss = gf.fit(diff_data);


		// calculate the probabilities and write them to the IDs
		for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      if (it->getHits().size() > 0)
      {
				vector<PeptideHit> hits;
				for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
        {
					PeptideHit hit = *pit;
					hit.setScore(getProbability_(result_gamma, rev_trafo, result_gauss, fwd_trafo, hit.getScore()));
					hits.push_back(hit);
				}
				PeptideIdentification id = *it;
				id.setScoreType(id.getScoreType() + "_DecoyProbability");
				id.setHits(hits);
				prob_ids.push_back(id);
			}
		}
	}
	
	// normalize the bins to [0, 1]
	void IDDecoyProbability::normalizeBins_(const vector<double>& scores, HashMap<double, double>& binned, Transformation_& trafo)
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

  	// perform the binning
  	double diff = max - min;
  	HashMap<UInt, UInt> bins;
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
  	for (HashMap<UInt, UInt>::ConstIterator it = bins.begin(); it != bins.end(); ++it)
  	{
    	binned[it->first/number_of_bins] = it->second / (double)max_bin * 4.0; // 4 is best value for the gamma distribution
  	}

		// store the transformation
  	trafo.y_factor = 4.0 / (double)max_bin;
  	trafo.x_factor = diff;
  	trafo.x_shift = min;
  	trafo.y_max_bin = max_bin_number;
  	trafo.x_max = max;

	}

	double IDDecoyProbability::getProbability_(const GammaDistributionFitter::GammaDistributionFitResult& result_gamma,
																						const Transformation_& gamma_trafo,
																						const GaussFitter::GaussFitResult& result_gauss,
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
  	if (score_fwd_trans > gamma_trafo.x_max)
  	{
    	return 1;
  	}

  	rho_fwd = result_gauss.A * exp(- pow(score_fwd_trans - result_gauss.x0, 2)/2.0/pow(result_gauss.sigma, 2));

  	// calc P using Bayes theorem
  	return rho_fwd / (rho_fwd + rho_rev);
	}

	
} // namespace OpenMS
