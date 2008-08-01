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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <algorithm>

#define FALSE_DISCOVERY_RATE_DEBUG
#undef  FALSE_DISCOVERY_RATE_DEBUG

using namespace std;

namespace OpenMS 
{
	FalseDiscoveryRate::FalseDiscoveryRate()
		: DefaultParamHandler("FalseDiscoveryRate")
	{
		
		defaultsToParam_();
	}
	
	void FalseDiscoveryRate::apply(vector<PeptideIdentification>& fwd_ids, vector<PeptideIdentification>& rev_ids)
	{
		if (fwd_ids.size() == 0 || rev_ids.size() == 0)
		{
			return;
		}
		vector<double> fwd_scores, rev_scores;
		// get the scores of all peptide hits
		for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
		{
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				fwd_scores.push_back(pit->getScore());
			}
		}
		for (vector<PeptideIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        rev_scores.push_back(pit->getScore());
      }
    }

		bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());
		
		// sort the scores
		if (higher_score_better)
		{
			sort(fwd_scores.rbegin(), fwd_scores.rend());
			sort(rev_scores.rbegin(), rev_scores.rend());
		}
		else
		{
			sort(fwd_scores.begin(), fwd_scores.end());
			sort(rev_scores.begin(), rev_scores.end());
		}

#ifdef FALSE_DISCOVERY_RATE_DEBUG
		cerr << "fwd-scores: " << endl;
		for (vector<double>::const_iterator it = fwd_scores.begin(); it != fwd_scores.end(); ++it)
		{
			cerr << *it << ", ";
		}
		cerr << endl;
		cerr << "rev-scores: " << endl;
    for (vector<double>::const_iterator it = rev_scores.begin(); it != rev_scores.end(); ++it)
    {
      cerr << *it << ", ";
    }
    cerr << endl;
#endif

		// calculate fdr for the forward scores
		Map<double, double> score_to_fdr;
		UInt j = 0;
		for (UInt i = 0; i != fwd_scores.size(); ++i)
		{
			while (j != rev_scores.size() && 
						 ((fwd_scores[i] <= rev_scores[j] && higher_score_better) ||
						 (fwd_scores[i] >= rev_scores[j] && !higher_score_better)))
			{
				++j;
			}

#ifdef FALSE_DISCOVERY_RATE_DEBUG
			cerr << fwd_scores[i] << " " << rev_scores[j] << " " << i << " " << j << " ";
#endif
			
			double fdr(0);

			if (j != 0 && i != 0)
			{
				fdr = (double)j / (double)(i + j);
			}
#ifdef FALSE_DISCOVERY_RATE_DEBUG
			cerr << fdr << endl;
#endif
			score_to_fdr[fwd_scores[i]] = fdr;
		}

		// annotate fdr 
		String score_type = fwd_ids.begin()->getScoreType() + "_score";
		for (vector<PeptideIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
		{
			it->setScoreType("FDR");
			it->setHigherScoreBetter(false);
			vector<PeptideHit> hits = it->getHits();
			for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
			{
#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << pit->getScore() << " " << score_to_fdr[pit->getScore()] << endl;
#endif
				pit->setMetaValue(score_type, pit->getScore());
				pit->setScore(score_to_fdr[pit->getScore()]);
			}
			it->setHits(hits);
		}
		
		return;
	}

	void FalseDiscoveryRate::apply(vector<ProteinIdentification>& fwd_ids, vector<ProteinIdentification>& rev_ids)
	{
    if (fwd_ids.size() == 0 || rev_ids.size() == 0)
    {
      return;
    }
    vector<double> fwd_scores, rev_scores;
    // get the scores of all peptide hits
    for (vector<ProteinIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        fwd_scores.push_back(pit->getScore());
      }
    }
    for (vector<ProteinIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        rev_scores.push_back(pit->getScore());
      }
    }

    bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());

    // sort the scores
    if (higher_score_better)
    {
      sort(fwd_scores.rbegin(), fwd_scores.rend());
      sort(rev_scores.rbegin(), rev_scores.rend());
    }
    else
    {
      sort(fwd_scores.begin(), fwd_scores.end());
      sort(rev_scores.begin(), rev_scores.end());
    }

   	// calculate fdr for the forward scores
    Map<double, double> score_to_fdr;
    UInt j = 0;
    for (UInt i = 0; i != fwd_scores.size(); ++i)
    {
      while (j != rev_scores.size() &&
             ((fwd_scores[i] <= rev_scores[j] && higher_score_better) ||
             (fwd_scores[i] >= rev_scores[j] && !higher_score_better)))
      {
        ++j;
      }
     	double fdr(0);

      if (j != 0 && i != 0)
      {
        fdr = (double)j / (double)(i + j);
      }
      score_to_fdr[fwd_scores[i]] = fdr;
    }

    // annotate fdr
    String score_type = fwd_ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      it->setScoreType("FDR");
      it->setHigherScoreBetter(false);
      vector<ProteinHit> hits = it->getHits();
      for (vector<ProteinHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
      {
        pit->setMetaValue(score_type, pit->getScore());
        pit->setScore(score_to_fdr[pit->getScore()]);
      }
      it->setHits(hits);
    }

		return;
	}

} // namespace OpenMS
