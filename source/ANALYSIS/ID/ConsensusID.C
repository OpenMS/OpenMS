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
		: DefaultParamHandler("ConsensusID"),
			inverse_order_(false)
	{
		defaults_.setValue("Algorithm","Ranked");
		defaults_.setValue("ConsideredHits","10");
		defaults_.setValue("NumberOfRuns","0");
		defaults_.setValue("InverseOrder","0");
		defaults_.setValue("MinOutputScore","0");
		
		defaultsToParam_();
	}
	
	void ConsensusID::apply(vector<Identification>& ids) throw (Exception::InvalidValue)
	{
		//Abort if no IDs present
		if (ids.size()==0)
		{
			return;
		}
		
		//check order
		if ((UInt)(param_.getValue("InverseOrder"))!=0)
		{
			inverse_order_ = true;
		}
		
		String algorithm = param_.getValue("Algorithm");
		
		if (algorithm == "Ranked")
		{
			ranked_(ids);
			ids[0].assignRanks();
		}
		else if (algorithm == "Merge")
		{	
			merge_(ids);
			ids[0].assignRanks(inverse_order_);
		}
		else if (algorithm == "Average")
		{	
			average_(ids);
			ids[0].assignRanks(inverse_order_);
		}
		else
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There is no such ConsensusID algorithm!",algorithm);
		}
		
#ifdef DEBUG_ID_CONSENSUS
		vector<PeptideHit>& hits2 = ids[0].getPeptideHits();
		for (UInt i=0; i< hits2.size(); ++i)
		{
			cout << "  " << hits2[i].getSequence() << " " << hits2[i].getScore() << endl;
		}
#endif
	}

	void ConsensusID::ranked_(vector<Identification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("ConsideredHits"));
		UInt number_of_runs = (UInt)(param_.getValue("NumberOfRuns"));
		
		//iterate over the different ID runs
		for (vector<Identification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - ID run" << endl;
#endif
			//make sure that the ranks are present
			id->assignRanks(inverse_order_);
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getPeptideHits().begin(); hit != id->getPeptideHits().end() && hit_count <= considered_hits; ++hit)
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
		Real min_score = (Real)(param_.getValue("MinOutputScore"));
		ids.clear();
		ids.resize(1);
		vector<PeptideHit>& hits = ids[0].getPeptideHits();
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			if (it->second >= min_score)
			{
				PeptideHit hit;
				hit.setScoreType("Consensus_averaged");
				hit.setSequence(it->first);
				hit.setScore(it->second);
				hits.push_back(hit);
			}
		}

	}

	void ConsensusID::merge_(vector<Identification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("ConsideredHits"));
				
		//store the score type (to make sure only IDs of the same type are merged)
		String score_type = ids[0].getPeptideHits().begin()->getScoreType();
		//iterate over the different ID runs
		for (vector<Identification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
			//check the score type
			if (id->getPeptideHits().begin()->getScoreType()!=score_type)
			{
				cerr << "Warning: You are merging differnt types of scores: '" << score_type << "' and '" << id->getPeptideHits().begin()->getScoreType() << "'" << endl;
			}
			//make sure that the ranks are present
			id->assignRanks(inverse_order_);
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getPeptideHits().begin(); hit != id->getPeptideHits().end() && hit_count <= considered_hits; ++hit)
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
		Real min_score = (Real)(param_.getValue("MinOutputScore"));
		ids.clear();
		ids.resize(1);
		vector<PeptideHit>& hits = ids[0].getPeptideHits();
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			if (it->second >= min_score)
			{
				PeptideHit hit;
				hit.setScoreType("Consensus_averaged");
				hit.setSequence(it->first);
				hit.setScore(it->second);
				hits.push_back(hit);
			}
		}
	}

	void ConsensusID::average_(vector<Identification>& ids)
	{
		map<String,Real> scores;		
		UInt considered_hits = (UInt)(param_.getValue("ConsideredHits"));
		UInt number_of_runs = (UInt)(param_.getValue("NumberOfRuns"));
		
		//iterate over the different ID runs
		for (vector<Identification>::iterator id = ids.begin(); id != ids.end(); ++id)
		{
#ifdef DEBUG_ID_CONSENSUS
			//cout << " - ID run" << endl;
#endif
			//make sure that the ranks are present
			id->assignRanks(inverse_order_);
			//iterate over the hits
			UInt hit_count = 1;
			for (vector<PeptideHit>::const_iterator hit = id->getPeptideHits().begin(); hit != id->getPeptideHits().end() && hit_count <= considered_hits; ++hit)
			{
				if (scores.find(hit->getSequence())==scores.end())
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - New hit: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
					scores.insert(make_pair(hit->getSequence(),hit->getScore()));  
				}
				else
				{
#ifdef DEBUG_ID_CONSENSUS
					//cout << " - Summed up: " << hit->getSequence() << " " << hit->getScore() << endl;
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
		Real min_score = (Real)(param_.getValue("MinOutputScore"));
		ids.clear();
		ids.resize(1);
		vector<PeptideHit>& hits = ids[0].getPeptideHits();
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			if (it->second >= min_score)
			{
				PeptideHit hit;
				hit.setScoreType("Consensus_averaged");
				hit.setSequence(it->first);
				hit.setScore(it->second);
				hits.push_back(hit);
			}
		}
	}

} // namespace OpenMS
