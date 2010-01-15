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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
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
 		defaults_.setValue("q_value", "true", "if true, the q-values will be calculated instead of the FDRs");
		defaults_.setValidStrings("q_value", StringList::create("true,false"));
		defaults_.setValue("use_all_hits", "false", "if true not only the first hit, but all are used");
		defaults_.setValidStrings("use_all_hits", StringList::create("true,false"));
		defaults_.setValue("decoy_string", "_rev", "String which is appended at the accession of the protein to indicate that it is a decoy protein.");
		defaultsToParam_();
	}

	void FalseDiscoveryRate::apply(vector<PeptideIdentification>& ids)
	{
		bool q_value = param_.getValue("q_value").toBool();
		bool use_all_hits = param_.getValue("use_all_hits").toBool();
		if (ids.size() == 0)
		{
			return;
		}

		// get the scores of all peptide hits
		vector<DoubleReal> target_scores, decoy_scores;
		for (vector<PeptideIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
		{
			for (Size i = 0; i < it->getHits().size(); ++i)
			{
				if (!it->getHits()[i].metaValueExists("target_decoy"))
				{
					cerr << "FalseDiscoveryRate: error, meta value 'target_decoy' does not exists, reindex the idXML file with PeptideIndexer first!" << endl;
					continue;
				}

				String target_decoy(it->getHits()[i].getMetaValue("target_decoy"));
				if (!use_all_hits && i > 0)
				{
					break;
				}

				if (target_decoy == "target")
				{
					target_scores.push_back(it->getHits()[i].getScore());
				}
				else 
				{
					if (target_decoy == "decoy" || target_decoy == "target+decoy")
					{
						decoy_scores.push_back(it->getHits()[i].getScore());
					}
					else
					{
						if (target_decoy != "")
						{
							cerr << "FalseDiscoveryRate: error, unknown value of meta value 'target_decoy': '" << target_decoy << "'!" << endl;
						}
					}
				}
			}
		}
		Size number_of_target_scores = target_scores.size();

		
		//cerr << "FalseDiscoveryRate: #target sequences=" << target_scores.size() << ", #decoy sequences=" << decoy_scores.size() << endl;

		if (decoy_scores.size() == 0)
		{
			cerr << "FalseDiscoveryRate: #decoy sequences is zero! Cannot proceed" << endl;
		}

		if (target_scores.size() == 0)
		{
			cerr << "FalseDiscoveryRate: #target sequences is zero! Cannot proceed" << endl;
		}

		if (target_scores.size() == 0 || decoy_scores.size() == 0)
		{
			return;
		}

    bool higher_score_better(ids.begin()->isHigherScoreBetter());

    // sort the scores
    if (higher_score_better && !q_value)
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else if (!higher_score_better && !q_value)
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.begin(), decoy_scores.end());
    }
    else if (higher_score_better)
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.begin(), decoy_scores.end());
    }

		// calculate fdr for the forward scores
    Map<DoubleReal, DoubleReal> score_to_fdr;
    Size j = 0;
    DoubleReal minimal_fdr = 1.;

    if (q_value)
    {
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        if (i == 0 && j == 0)
        {
          while (j != decoy_scores.size()
                 &&((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
                    (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
          {
            ++j;
          }
        }
        else
        {
          if (j == decoy_scores.size())
          {
            j--;
          }
          while (j != 0
                 &&((target_scores[i] > decoy_scores[j] && higher_score_better) ||
                    (target_scores[i] < decoy_scores[j] && !higher_score_better)))
          {
            --j;
          }
          // Since j has to be equal to the number of fps above the threshold we add one
          if ((target_scores[i] <= decoy_scores[j] && higher_score_better)
              || (target_scores[i] >= decoy_scores[j] && !higher_score_better))
          {
            ++j;
          }
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif

        DoubleReal fdr = 0.;

        if (minimal_fdr >= (DoubleReal)j / (number_of_target_scores - i))
        {
          minimal_fdr = (DoubleReal)j / (number_of_target_scores - i);
        }
        fdr = minimal_fdr;

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << fdr << endl;
#endif
        score_to_fdr[target_scores[i]] = fdr;

      }
    }
    else
    {
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        while (j != decoy_scores.size() &&
               ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
               (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
        {
          ++j;
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif

        DoubleReal fdr(0);

        fdr = (DoubleReal)j / (DoubleReal)(i + 1);

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << fdr << endl;
#endif
        score_to_fdr[target_scores[i]] = fdr;
      }
    }


		// annotate fdr
    String score_type = ids.begin()->getScoreType() + "_score";
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      if (q_value)
      {
        it->setScoreType("q-value");
      }
      else
      {
        it->setScoreType("FDR");
      }

      it->setHigherScoreBetter(false);
      vector<PeptideHit> hits;
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
				PeptideHit hit = *pit;
				if (hit.metaValueExists("target_decoy"))
				{
					String meta_value = (String)hit.getMetaValue("target_decoy");
					if (meta_value == "decoy" || meta_value == "target+decoy")
					{
						continue;
					}
				}
        hit.setMetaValue(score_type, pit->getScore());
        hit.setScore(score_to_fdr[pit->getScore()]);
				hits.push_back(hit);
      }
      it->setHits(hits);
    }
    return;
	}

	void FalseDiscoveryRate::apply(vector<PeptideIdentification>& fwd_ids, vector<PeptideIdentification>& rev_ids)
	{
		bool q_value = param_.getValue("q_value").toBool();
		
		if (fwd_ids.size() == 0 || rev_ids.size() == 0)
		{
			return;
		}
		vector<DoubleReal> target_scores, decoy_scores;
		// get the scores of all peptide hits
		for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
		{
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				target_scores.push_back(pit->getScore());
			}
		}
		Size number_of_target_scores = target_scores.size();

		for (vector<PeptideIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        decoy_scores.push_back(pit->getScore());
      }
    }

		bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());
		
		// sort the scores
		if (higher_score_better && !q_value)
		{
			sort(target_scores.rbegin(), target_scores.rend());
			sort(decoy_scores.rbegin(), decoy_scores.rend());
		}
		else if (!higher_score_better && !q_value)
		{
			sort(target_scores.begin(), target_scores.end());
			sort(decoy_scores.begin(), decoy_scores.end());
		}
		else if (higher_score_better)
		{
			sort(target_scores.begin(), target_scores.end());
			sort(decoy_scores.rbegin(), decoy_scores.rend());
		}
		else
		{
			sort(target_scores.rbegin(), target_scores.rend());
			sort(decoy_scores.begin(), decoy_scores.end());
		}

#ifdef FALSE_DISCOVERY_RATE_DEBUG
		cerr << "fwd-scores: " << endl;
		for (vector<DoubleReal>::const_iterator it = target_scores.begin(); it != target_scores.end(); ++it)
		{
			cerr << *it << ", ";
		}
		cerr << endl;
		cerr << "rev-scores: " << endl;
    for (vector<DoubleReal>::const_iterator it = decoy_scores.begin(); it != decoy_scores.end(); ++it)
    {
      cerr << *it << ", ";
    }
    cerr << endl;
#endif

		// calculate fdr for the forward scores
		Map<DoubleReal, DoubleReal> score_to_fdr;
		Size j = 0;
		DoubleReal minimal_fdr = 1.;
		
		if (q_value)
		{
			for (Size i = 0; i != target_scores.size(); ++i)
			{
				if (i == 0 && j == 0)
				{
					while (j != decoy_scores.size()
								 &&((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
								 		(target_scores[i] >= decoy_scores[j] && !higher_score_better)))
					{
						++j;
					}
				}
				else
				{
					if (j == decoy_scores.size())
					{
						j--;
					}
					while (j != 0
								 &&((target_scores[i] > decoy_scores[j] && higher_score_better) ||
								 		(target_scores[i] < decoy_scores[j] && !higher_score_better)))
					{
						--j;
					}
					// Since j has to be equal to the number of fps above the threshold we add one
					if ((target_scores[i] <= decoy_scores[j] && higher_score_better)
							|| (target_scores[i] >= decoy_scores[j] && !higher_score_better))
					{
						++j;					
					}
				}	

#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif

				DoubleReal fdr = 0.;

				if (minimal_fdr >= (DoubleReal)j / (number_of_target_scores - i))
				{
					minimal_fdr = (DoubleReal)j / (number_of_target_scores - i);
				}
				fdr = minimal_fdr;

#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << fdr << endl;
#endif
				score_to_fdr[target_scores[i]] = fdr;

			}					
		}
		else
		{
			for (Size i = 0; i != target_scores.size(); ++i)
			{
				while (j != decoy_scores.size() && 
							 ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
							 (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
				{
					++j;
				}
	
#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif
				
				DoubleReal fdr(0);
	
				fdr = (DoubleReal)j / (DoubleReal)(i + 1);

#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << fdr << endl;
#endif
				score_to_fdr[target_scores[i]] = fdr;
			}
		}
		// annotate fdr 
		String score_type = fwd_ids.begin()->getScoreType() + "_score";
		for (vector<PeptideIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
		{
			if (q_value)
			{
				it->setScoreType("q-value");
			}
			else
			{
				it->setScoreType("FDR");
			}
				
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

	void FalseDiscoveryRate::apply(vector<ProteinIdentification>& ids)
	{
		if (ids.size() == 0) 
		{
			return;
		}

		vector<DoubleReal> target_scores, decoy_scores;
		String decoy_string = (String)param_.getValue("decoy_string");
		for (vector<ProteinIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
		{
			for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				if (pit->getAccession().hasSubstring(decoy_string))
				{
					decoy_scores.push_back(pit->getScore());
				}
				else
				{
					target_scores.push_back(pit->getScore());
				}
			}
		}

    bool higher_score_better(ids.begin()->isHigherScoreBetter());

    // sort the scores
    if (higher_score_better)
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.begin(), decoy_scores.end());
    }

    // calculate fdr for the forward scores
    Map<DoubleReal, DoubleReal> score_to_fdr;
    Size j = 0;
    for (Size i = 0; i != target_scores.size(); ++i)
    {
      while (j != decoy_scores.size() &&
             ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
             (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
      {
        ++j;
      }
      DoubleReal fdr(0);
      fdr = (DoubleReal)j / ((DoubleReal)(i) + 1);
      score_to_fdr[target_scores[i]] = fdr;
			if (j < decoy_scores.size())
			{
				score_to_fdr[decoy_scores[j]] = fdr;
			}
    }

    // annotate fdr
    String score_type = ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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

	void FalseDiscoveryRate::apply(vector<ProteinIdentification>& fwd_ids, vector<ProteinIdentification>& rev_ids)
	{
    if (fwd_ids.size() == 0 || rev_ids.size() == 0)
    {
      return;
    }
    vector<DoubleReal> target_scores, decoy_scores;
    // get the scores of all peptide hits
    for (vector<ProteinIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        target_scores.push_back(pit->getScore());
      }
    }
    for (vector<ProteinIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        decoy_scores.push_back(pit->getScore());
      }
    }

    bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());

    // sort the scores
    if (higher_score_better)
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.begin(), decoy_scores.end());
    }

   	// calculate fdr for the forward scores
    Map<DoubleReal, DoubleReal> score_to_fdr;
    Size j = 0;
    for (Size i = 0; i != target_scores.size(); ++i)
    {
      while (j != decoy_scores.size() &&
             ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
             (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
      {
        ++j;
      }
     	DoubleReal fdr(0);
      fdr = (DoubleReal)j / ((DoubleReal)(i) + 1);
      score_to_fdr[target_scores[i]] = fdr;
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
