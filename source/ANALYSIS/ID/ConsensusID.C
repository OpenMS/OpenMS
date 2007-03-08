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
#include <OpenMS/KERNEL/Feature.h>

#include <map>

using namespace std;

namespace OpenMS 
{
	ConsensusID::ConsensusID()
		: DefaultParamHandler("ConsensusID")
	{
		defaults_.setValue("Algorithm","Ranked");
		defaults_.setValue("ConsideredHits","10");
		
		defaultsToParam_();
	}
	
	void ConsensusID::apply(Feature& feature) throw (Exception::InvalidValue)
	{
#ifdef DEBUG_ID
		cout << "Feature " << feature.getRT() << " / " << feature.getMZ() << " -> " << feature.getIdentifications().size() << endl;
#endif
		//Abort if no IDs present
		if (feature.getIdentifications().size()==0)
		{
			return;
		}
		
		String algorithm = param_.getValue("Algorithm");

		//temporary datastructure to accumulate scores
		map<String,Real> scores;
		
		UInt considered_hits = (UInt)(param_.getValue("ConsideredHits"));
		
		if (algorithm == "Ranked")
		{
			//iterate over the different ID runs
			vector<Identification>& ids = feature.getIdentifications();
			for (vector<Identification>::iterator id = ids.begin(); id != ids.end(); ++id)
			{
				//make sure that the ranks are present
				id->assignRanks();
				//iterate over the hits
				UInt hit_count = 1;
				for (vector<PeptideHit>::const_iterator hit = id->getPeptideHits().begin(); hit != id->getPeptideHits().end() && hit_count <= considered_hits; ++hit)
				{
					if (scores.find(hit->getSequence())==scores.end())
					{
						scores.insert(make_pair(hit->getSequence(),( considered_hits + 1 - hit->getRank())));  
					}
					else
					{
						scores[hit->getSequence()] += ( considered_hits + 1 - hit->getRank());  
					}
					++hit_count;
				}
			}
		}
		else if (algorithm == "Merge")
		{	
			//store the score type (to make sure only IDs of the same type are merged)
			String score_type = feature.getIdentifications()[0].getPeptideHits().begin()->getScoreType();
			//iterate over the different ID runs
			vector<Identification>& ids = feature.getIdentifications();
			for (vector<Identification>::iterator id = ids.begin(); id != ids.end(); ++id)
			{
#ifdef DEBUG_ID
				cout << "Id run with " << id->getPeptideHits().size() << " hits" << endl;
#endif
				//check the score type
				if (id->getPeptideHits().begin()->getScoreType()!=score_type)
				{
					cerr << "Warning: You are merging differnt types of scores: '" << score_type << "' and '" << id->getPeptideHits().begin()->getScoreType() << "'" << endl;
				}
				//make sure that the ranks are present
				id->assignRanks();
				//iterate over the hits
				UInt hit_count = 1;
				for (vector<PeptideHit>::const_iterator hit = id->getPeptideHits().begin(); hit != id->getPeptideHits().end() && hit_count <= considered_hits; ++hit)
				{
					if (scores.find(hit->getSequence())==scores.end())
					{
#ifdef DEBUG_ID
						cout << " - Added hit: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
						scores.insert(make_pair(hit->getSequence(),hit->getScore()));  
					}
					else if(scores[hit->getSequence()] < hit->getScore())
					{
#ifdef DEBUG_ID
						cout << " - Updated hit: " << hit->getSequence() << " " << scores[hit->getSequence()] << " -> " << hit->getScore() << endl;
#endif
						scores[hit->getSequence()] = hit->getScore();
					}
					++hit_count;
				}
			}
		}
		else
		{
			throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There is no such ConsensusID algorithm!",algorithm);
		}
		
		//Replace IDs by consensus
		feature.getIdentifications().clear();
		feature.getIdentifications().resize(1);
		vector<PeptideHit>& hits = feature.getIdentifications()[0].getPeptideHits();
		hits.resize(scores.size());
		UInt hit_count = 0;
#ifdef DEBUG_ID
		cout << "Writing " << scores.size() << " hits" << endl;
#endif
		for (map<String,Real>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			if (algorithm == "Ranked")
			{
				hits[hit_count].setScoreType("Consensus_ranked");
			}
			else if (algorithm == "Merge")
			{
				hits[hit_count].setScoreType("Consensus_merge");
			}
			else
			{
				throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"No Score type defined for this algorithm!",algorithm);
			}
			hits[hit_count].setSequence(it->first);
			hits[hit_count].setScore(it->second);
			++hit_count;
		}
		feature.getIdentifications()[0].assignRanks();
	}
	
} // namespace OpenMS
