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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

#include <map>

using namespace std;

namespace OpenMS 
{
	ConsensusID::ConsensusID()
		: DefaultParamHandler("ConsensusID")
	{
		defaults_.setValue("algorithm","ranked","Allowed algorithm names are 'ranked', 'merge' and 'average'.\n"
											 "merge -- merges the runs with respect to their score. The score is not modified. Make sure to use PeptideIdentifications with the same score type only!\n"
										   "ranked -- reorders the hits according to a consensus score computed from the ranks in the input runs. The score is normalized to the interval (0,100). The PeptideIdentifications do not need to have the same score type.\n"
										   "average -- reorders the hits according to the average score of the input runs. Make sure to use PeptideIdentifications with the same score type only!"
										   , false);
  	
		defaults_.setValue("considered_hits",10,"The number of top hits that are used for the consensus scoring.", false);
		defaults_.setValue("number_of_runs",0,"The number of runs used as input. This information is used in 'Ranked' and 'Average' to compute the new scores. If not given, the number of input identifications is taken.", false);
		
		defaultsToParam_();
	}
	
	void ConsensusID::apply(vector<PeptideIdentification>& ids) throw (Exception::InvalidValue)
	{
		//Abort if no IDs present
		if (ids.size()==0)
		{
			return;
		}
		
		String algorithm = param_.getValue("algorithm");
		
		if (algorithm == "ranked")
		{
			ranked_(ids);
			ids[0].assignRanks();
		}
		else if (algorithm == "merge")
		{	
			merge_(ids);
			ids[0].assignRanks();
		}
		else if (algorithm == "average")
		{	
			average_(ids);
			ids[0].assignRanks();
		}
		else
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There is no such ConsensusID algorithm!",algorithm);
		}
		
#ifdef DEBUG_ID_CONSENSUS
		const vector<PeptideHit>& hits2 = ids[0].getHits();
		for (UInt i=0; i< hits2.size(); ++i)
		{
			cout << "  " << hits2[i].getSequence() << " " << hits2[i].getScore() << endl;
		}
#endif
	}

	void ConsensusID::ranked_(vector<PeptideIdentification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
		UInt number_of_runs = (UInt)(param_.getValue("number_of_runs"));
		
		//iterate over the different ID runs
		for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - ID run" << endl;
#endif
			//make sure that the ranks are present
			id->assignRanks();
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
			{
				if (scores.find(hit->getSequence())==scores.end())
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - New hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
					scores.insert(make_pair(hit->getSequence(),( considered_hits + 1 - hit->getRank())));  
				}
				else
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - Added hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
					scores[hit->getSequence()] += ( considered_hits + 1 - hit->getRank());  
				}
				++hit_count;
			}
		}
		//divide score by the maximum possible score and multiply with 100 => %
		UInt max_score;
		if (number_of_runs==0)
		{
			max_score = ids.size()*considered_hits;
		}
		else
		{
			max_score = number_of_runs*considered_hits;
		}
		for (map<String,Real>::iterator it = scores.begin(); it != scores.end(); ++it)
		{
			it->second = (it->second * 100.0f / max_score);
		}

		//Replace IDs by consensus
		ids.clear();
		ids.resize(1);
		ids[0].setScoreType("Consensus_averaged");

		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second);
			ids[0].insertHit(hit);
		}

	}

	void ConsensusID::merge_(vector<PeptideIdentification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
				
		//store the score type (to make sure only IDs of the same type are merged)
		String score_type = ids[0].getScoreType();
		bool higher_better = ids[0].isHigherScoreBetter();
		
		//iterate over the different ID runs
		for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
			//check the score type
			if (id->getScoreType()!=score_type)
			{
				cerr << "Warning: You are merging different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
			}
			if (id->isHigherScoreBetter()!=higher_better)
			{
				cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
			}
			//make sure that the ranks are present
			id->assignRanks();
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
			{
				if (scores.find(hit->getSequence())==scores.end())
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - Added hit: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
					scores.insert(make_pair(hit->getSequence(),hit->getScore()));  
				}
				else if(scores[hit->getSequence()] < hit->getScore())
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - Updated hit: " << hit->getSequence() << " " << scores[hit->getSequence()] << " -> " << hit->getScore() << endl;
#endif
					scores[hit->getSequence()] = hit->getScore();
				}
				++hit_count;
			}
		}

		//Replace IDs by consensus
		ids.clear();
		ids.resize(1);
		ids[0].setScoreType(String("Consensus_merged (") + score_type +")");
		ids[0].setHigherScoreBetter(higher_better);
		
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second);
			ids[0].insertHit(hit);
		}
	}

	void ConsensusID::average_(vector<PeptideIdentification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
		UInt number_of_runs = (UInt)(param_.getValue("number_of_runs"));
		
		//store the score type (to make sure only IDs of the same type are averaged)
		String score_type = ids[0].getScoreType();
		bool higher_better = ids[0].isHigherScoreBetter();
				
		//iterate over the different ID runs
		for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
#ifdef DEBUG_ID_CONSENSUS
			cout << " - ID run" << endl;
#endif
			//make sure that the ranks are present
			id->assignRanks();
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
			{
				//check the score type
				if (id->getScoreType()!=score_type)
				{
					cerr << "Warning: You are averaging different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
				}
				if (id->isHigherScoreBetter()!=higher_better)
				{
					cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
				}
				if (scores.find(hit->getSequence())==scores.end())
				{
#ifdef DEBUG_ID_CONSENSUS
					cout << " - New hit: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
					scores.insert(make_pair(hit->getSequence(),hit->getScore()));  
				}
				else
				{
#ifdef DEBUG_ID_CONSENSUS
					cout << " - Summed up: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
					scores[hit->getSequence()] += hit->getScore();  
				}
				++hit_count;
			}
		}
		//normalize score by number of id runs
		for (map<String,Real>::iterator it = scores.begin(); it != scores.end(); ++it)
		{
			if (number_of_runs==0)
			{
				it->second = (it->second / ids.size());
			}
			else
			{
				it->second = (it->second / number_of_runs);
			}
		}

		//Replace IDs by consensus
		ids.clear();
		ids.resize(1);
		ids[0].setScoreType(String("Consensus_averaged (") + score_type +")");
		ids[0].setHigherScoreBetter(higher_better);
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second);
			ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
			cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

		}
	}

} // namespace OpenMS
