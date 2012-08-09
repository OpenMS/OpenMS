// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>

#define FALSE_DISCOVERY_RATE_DEBUG
#undef  FALSE_DISCOVERY_RATE_DEBUG

using namespace std;

namespace OpenMS 
{
	FalseDiscoveryRate::FalseDiscoveryRate()
		: DefaultParamHandler("FalseDiscoveryRate")
	{		
 		defaults_.setValue("q_value", "true", "If 'true', the q-values will be calculated instead of the FDRs");
		defaults_.setValidStrings("q_value", StringList::create("true,false"));
		defaults_.setValue("use_all_hits", "false", "If 'true' not only the first hit, but all are used (peptides only)");
		defaults_.setValidStrings("use_all_hits", StringList::create("true,false"));
		defaults_.setValue("split_charge_variants", "false", "If set to 'true' charge variants are treated separately (for peptides of combined target/decoy searches only).");
		defaults_.setValidStrings("split_charge_variants", StringList::create("true,false"));
		defaults_.setValue("treat_runs_separately", "false", "If set to 'true' different search runs are treated separately (for peptides of combined target/decoy searches only).");
		defaults_.setValidStrings("treat_runs_separately", StringList::create("true,false"));
		defaults_.setValue("decoy_string", "_rev", "String which is appended at the accession of the protein to indicate that it is a decoy protein (for proteins only).");
		defaults_.setValue("add_decoy_peptides","false", "If set to true, decoy peptides will be written to output file, too. The q-value is set to the closest target score.");
		defaults_.setValidStrings("add_decoy_peptides", StringList::create("true,false"));
		defaultsToParam_();
	}

	void FalseDiscoveryRate::apply(vector<PeptideIdentification>& ids)
	{
		bool q_value = param_.getValue("q_value").toBool();
		bool use_all_hits = param_.getValue("use_all_hits").toBool();
		bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
		bool split_charge_variants = param_.getValue("split_charge_variants").toBool();
	  bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
#ifdef FALSE_DISCOVERY_RATE_DEBUG
		cerr << "Parameters: q_value=" << q_value << ", use_all_hits=" << use_all_hits << ", treat_runs_separately=" << treat_runs_separately << ", split_charge_variants=" << split_charge_variants << endl;
#endif


		if (ids.empty())
		{
      LOG_WARN << "No peptide identifications given to FalseDiscoveryRate! No calculation performed.\n";
			return;
		}

    // first search for all identifiers and charge variants
    set<String> identifiers;
		set<SignedSize> charge_variants;
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      identifiers.insert(it->getIdentifier());
			it->assignRanks();

			if (!use_all_hits)
			{
				vector<PeptideHit> hits = it->getHits();
				hits.resize(1);
				it->setHits(hits);
			}
		
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				charge_variants.insert(pit->getCharge());
			}
    }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
		cerr << "#id-runs: " << identifiers.size() << " ";
		for (set<String>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it)
		{
			cerr << "," << *it;
		}
		cerr << endl;


		cerr << "#of charge states: " << charge_variants.size() << " ";
		for (set<SignedSize>::const_iterator it = charge_variants.begin(); it != charge_variants.end(); ++it)
		{
			cerr << "," << *it;
		}
		cerr << endl;
#endif

		for (set<SignedSize>::const_iterator zit = charge_variants.begin(); zit != charge_variants.end(); ++zit)
		{
#ifdef FALSE_DISCOVERY_RATE_DEBUG
			cerr << "Charge variant=" << *zit << endl;
#endif

			// for all identifiers
			for (set<String>::const_iterator iit = identifiers.begin(); iit != identifiers.end(); ++iit)
			{
				if (!treat_runs_separately && iit != identifiers.begin())
				{
					continue;
				}

#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << "Id-run: " << *iit << endl;
#endif
				// get the scores of all peptide hits
				vector<DoubleReal> target_scores, decoy_scores;
				for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
				{
					// if runs should be treated separately, the identifiers must be the same
					if (treat_runs_separately && it->getIdentifier() != *iit)
					{
						continue;
					}

					for (Size i = 0; i < it->getHits().size(); ++i)
					{
						if (split_charge_variants && it->getHits()[i].getCharge() != *zit)
						{
							continue;
						}

						if (!it->getHits()[i].metaValueExists("target_decoy"))
						{
							LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' first (run-id='" << it->getIdentifier() << ", rank=" << i+1 << " of " << it->getHits().size() << ")!" << endl;
              throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Meta value 'target_decoy' does not exist!");
						}

						String target_decoy(it->getHits()[i].getMetaValue("target_decoy"));
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
									LOG_FATAL_ERROR << "Unknown value of meta value 'target_decoy': '" << target_decoy << "'!" << endl;
								}
							}
						}
					}
				}
		
#ifdef FALSE_DISCOVERY_RATE_DEBUG
				cerr << "#target-scores=" << target_scores.size() << ", #decoy-scores=" << decoy_scores.size() << endl;
#endif

				// check decoy scores
				if (decoy_scores.empty())
				{
					String error_string = "FalseDiscoveryRate: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0! ";
					if (split_charge_variants || treat_runs_separately)
					{
						error_string += "(";
						if (split_charge_variants)
						{
							error_string += "charge_variant=" + String(*zit) + " ";
						}
						if (treat_runs_separately)
						{
							error_string += "run-id=" + *iit;
						}
						error_string += ")";
					}
					LOG_ERROR << error_string << std::endl;
				}

				// check target scores
				if (target_scores.empty())
				{	
					String error_string = "FalseDiscoveryRate: #target sequences is zero! Ignoring. ";
          if (split_charge_variants || treat_runs_separately)
          {
            error_string += "(";
            if (split_charge_variants)
            {
              error_string += "charge_variant=" + String(*zit) + " ";
            }
            if (treat_runs_separately)
            {
              error_string += "run-id=" + *iit;
            }
            error_string += ")";
          }
          LOG_ERROR << error_string << std::endl;
				}

				if (target_scores.empty() || decoy_scores.empty())
				{
					// no remove the the relevant entries, or put 'pseudo-scores' in
	        for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
  	      {
    	      // if runs should be treated separately, the identifiers must be the same
      	    if (treat_runs_separately && it->getIdentifier() != *iit)
        	  {
          	  continue;
          	}

						vector<PeptideHit> hits(it->getHits()), new_hits;
	          for (Size i = 0; i < hits.size(); ++i)
  	        {
    	        if (split_charge_variants && hits[i].getCharge() != *zit)
      	      {
								new_hits.push_back(hits[i]);
        	      continue;
          	  }

	            if (!hits[i].metaValueExists("target_decoy"))
  	          {
    	          LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' first (run-id='" << it->getIdentifier() << ", rank=" << i+1 << " of " << hits.size() << ")!" << endl;
                throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Meta value 'target_decoy' does not exist!");
        	    }

       	     	String target_decoy(hits[i].getMetaValue("target_decoy"));
        	  	if (target_decoy == "target")
          	  {
								// if it is a target hit, there are now decoys, fdr/q-value should be zero then
								new_hits.push_back(hits[i]);
								String score_type = it->getScoreType() + "_score";
								new_hits.back().setMetaValue(score_type, new_hits.back().getScore());
								new_hits.back().setScore(0);
            	}
            	else
            	{
              	if (target_decoy != "decoy" && target_decoy != "target+decoy")
              	{
                  LOG_FATAL_ERROR << "Unknown value of meta value 'target_decoy': '" << target_decoy << "'!" << endl;
              	}
							}
						}
						it->setHits(new_hits);
					}
					continue;
				}

				// calculate fdr for the forward scores
 		   	bool higher_score_better(ids.begin()->isHigherScoreBetter());
 		   	Map<DoubleReal, DoubleReal> score_to_fdr;
				calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

				// annotate fdr
 		   	for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
   		 	{
  	   	  // if runs should be treated separately, the identifiers must be the same
   	     	if (treat_runs_separately && it->getIdentifier() != *iit)
   	     	{
   	      	continue;
        	}

 		   		String score_type = it->getScoreType() + "_score";
      		vector<PeptideHit> hits;
      		for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      		{
						PeptideHit hit = *pit;

						if (split_charge_variants && pit->getCharge() != *zit)
						{
							hits.push_back(*pit);
							continue;
						}
						if (hit.metaValueExists("target_decoy"))
						{
							String meta_value = (String)hit.getMetaValue("target_decoy");
							if ((meta_value == "decoy" || meta_value == "target+decoy") && !add_decoy_peptides)
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
			}
			if (!split_charge_variants)
			{
				break;
			}
		}

		// higher-score-better can be set now, calculations are finished
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
			if (q_value)
			{
				if (it->getScoreType() != "q-value")
				{
					it->setScoreType("q-value");
				}
			}
			else
			{
				if (it->getScoreType() != "FDR")
				{
					it->setScoreType("FDR");
				}
			}
      it->setHigherScoreBetter(false);
			it->assignRanks();
		}

    return;
	}

	void FalseDiscoveryRate::apply(vector<PeptideIdentification>& fwd_ids, vector<PeptideIdentification>& rev_ids)
	{
		if (fwd_ids.empty() || rev_ids.empty())
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

		for (vector<PeptideIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        decoy_scores.push_back(pit->getScore());
      }
    }

		bool q_value(param_.getValue("q_value").toBool());
		bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());
	  bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
		// calculate fdr for the forward scores
		Map<DoubleReal, DoubleReal> score_to_fdr;
		calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

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
		//write as well decoy peptides
		if(add_decoy_peptides)
		{
			score_type = rev_ids.begin()->getScoreType() + "_score";
			for (vector<PeptideIdentification>::iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
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
		}
		
		return;
	}

	void FalseDiscoveryRate::apply(vector<ProteinIdentification>& ids)
	{
		if (ids.empty()) 
		{
      LOG_WARN << "No protein identifications given to FalseDiscoveryRate! No calculation performed.\n";
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

		bool q_value(param_.getValue("q_value").toBool());
    bool higher_score_better(ids.begin()->isHigherScoreBetter());

    // calculate fdr for the forward scores
    Map<DoubleReal, DoubleReal> score_to_fdr;
		calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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
    if (fwd_ids.empty() || rev_ids.empty())
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

		bool q_value(param_.getValue("q_value").toBool());
    bool higher_score_better(fwd_ids.begin()->isHigherScoreBetter());
   	// calculate fdr for the forward scores
    Map<DoubleReal, DoubleReal> score_to_fdr;
		calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = fwd_ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
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


	void FalseDiscoveryRate::calculateFDRs_(Map<DoubleReal, DoubleReal>& score_to_fdr, vector<DoubleReal>& target_scores, vector<DoubleReal>& decoy_scores, bool q_value, bool higher_score_better)
	{
		Size number_of_target_scores = target_scores.size();
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

    Size j = 0;
    DoubleReal minimal_fdr = 1.;

    if (q_value)
    {
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        if (decoy_scores.empty())
        {
          // set FDR to 0 (done below automatically)
        }
        else if (i == 0 && j == 0)
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
	

   // assign	q-value of decoy_score to closest target_score
   for (Size i = 0; i != decoy_scores.size(); ++i)
   {	
	   Size closest_idx = 0;	
     for (Size j = 0; j != target_scores.size(); ++j)
     {
		 	 if (fabs(decoy_scores[i] - target_scores[j]) < fabs(decoy_scores[i] - target_scores[closest_idx]))
			 {
			   closest_idx = j;				 
			 }       
	   }
		 score_to_fdr[decoy_scores[i]] = score_to_fdr[target_scores[closest_idx]];
	 }

	}

} // namespace OpenMS
