// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen and others $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>


#include <map>
#include <cmath>

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

namespace OpenMS 
{
	ConsensusID::ConsensusID()
		: DefaultParamHandler("ConsensusID")
	{
		defaults_.setValue("algorithm","ranked","Algorithm used for the consensus scoring.\n"
											 "merge -- merges the runs with respect to their score. The score is not modified. Make sure to use PeptideIdentifications with the same score type only!\n"
										   "ranked -- reorders the hits according to a consensus score computed from the ranks in the input runs. The score is normalized to the interval (0,100). The PeptideIdentifications do not need to have the same score type.\n"
										   "average -- reorders the hits according to the average score of the input runs. Make sure to use PeptideIdentifications with the same score type only!\n"
										   "PEP -- calculates a consensus score based on posterior error probabilities. Make sure to use PeptideIdentifications with score types converted to PEPs only!\n");
  	defaults_.setValidStrings("algorithm",StringList::create("ranked,merge,average,PEP"));
		defaults_.setValue("considered_hits",10,"The number of top hits that are used for the consensus scoring.");
		defaults_.setMinInt("considered_hits",1);
		defaults_.setValue("number_of_runs",0,"The number of runs used as input. This information is used in 'Ranked' and 'Average' to compute the new scores. If not given, the number of input identifications is taken.");
		defaults_.setMinInt("number_of_runs",0);
		defaults_.setValue("mz_tolerance", 0.01, "give a a mz tolerance");
		defaults_.setValue("rt_tolerance", 0.01, "give an rt tolerance");
		defaults_.setValue("min_number_of_engines", 2, "The minimum number of search engines used to generate peptide lists");
		defaults_.setMinInt("min_number_of_engines",2);
//nice value for mz and RT tolerances 0.01
		
		defaultsToParam_();
	}
	
	void ConsensusID::apply(vector<PeptideIdentification>& ids)
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
		else if (algorithm == "PEP")
		{	
			PEP_(ids);
			ids[0].assignRanks();
		}
		else if (algorithm == "majority")
		{
			majority_(ids);
			ids[0].assignRanks();
		}
		
#ifdef DEBUG_ID_CONSENSUS
		const vector<PeptideHit>& hits2 = ids[0].getHits();
		for (Size i=0; i< hits2.size(); ++i)
		{
			cout << "  " << hits2[i].getSequence() << " " << hits2[i].getScore() << endl;
		}
#endif
	}

	void ConsensusID::ranked_(vector<PeptideIdentification>& ids)
	{
		cout << "ranked_" << endl;
		map<AASequence, DoubleReal> scores;		
		UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
		UInt number_of_runs = (UInt)(param_.getValue("number_of_runs"));
		
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
				if (scores.find(hit->getSequence())==scores.end())
				{
#ifdef DEBUG_ID_CONSENSUS
					cout << " - New hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
					scores.insert(make_pair(hit->getSequence(),DoubleReal( considered_hits + 1 - hit->getRank())));  
				}
				else
				{
#ifdef DEBUG_ID_CONSENSUS
					cout << " - Added hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
					scores[hit->getSequence()] += ( considered_hits + 1 - hit->getRank());  
				}
				++hit_count;
			}
		}
		//divide score by the maximum possible score and multiply with 100 => %
		Size max_score;
		if (number_of_runs==0)
		{
			max_score = ids.size()*considered_hits;
		}
		else
		{
			max_score = number_of_runs*considered_hits;
		}
		for (map<AASequence,DoubleReal>::iterator it = scores.begin(); it != scores.end(); ++it)
		{
			it->second = (it->second * 100.0f / max_score);
		}

		//Replace IDs by consensus
		ids.clear();
		ids.resize(1);
		ids[0].setScoreType("Consensus_averaged");

		for (map<AASequence,DoubleReal>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second);
			ids[0].insertHit(hit);
		}

	}

	void ConsensusID::merge_(vector<PeptideIdentification>& ids)
	{
		map<AASequence,DoubleReal> scores;		
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
		
		for (map<AASequence,DoubleReal>::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second);
			ids[0].insertHit(hit);
		}
	}

	void ConsensusID::average_(vector<PeptideIdentification>& ids)
	{
		map<AASequence,DoubleReal> scores;		
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
				if (scores.find(hit->getSequence())==scores.end())						//.end zeigt auf ein Element nach dem letzten
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
		for (map<AASequence,DoubleReal>::iterator it = scores.begin(); it != scores.end(); ++it)
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
		for (map<AASequence,DoubleReal>::const_iterator it = scores.begin(); it != scores.end(); ++it)
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
	
	
	
void ConsensusID::PEP_(vector<PeptideIdentification>& ids)
	{
		Map<AASequence,vector<DoubleReal> > scores;	

		UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
        //UInt number_of_runs = (UInt)(param_.getValue("numberOfRuns"));
			

		String score_type = ids[0].getScoreType();
		bool higher_better = ids[0].isHigherScoreBetter();

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
					cerr << "Warning: You are working with different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
				}
				if (id->isHigherScoreBetter()!=higher_better)
				{
					cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
				}
				DoubleReal a_score=(double)hit->getMetaValue("PEP");
				DoubleReal a_sim=1;
					
					
					set<String> myset;
					for(vector<PeptideHit>::const_iterator t = id->getHits().begin(); t != id->getHits().end(); ++t)
					{
						if(myset.find(t->getMetaValue("scoring"))==myset.end() & hit->getMetaValue("scoring") != t->getMetaValue("scoring"))
						{
							DoubleReal z=0;
							DoubleReal a=0;
							
							//find the same or most similar peptide sequence in lists from other search engines
							for(vector<PeptideHit>::const_iterator tt = id->getHits().begin(); tt != id->getHits().end(); ++tt)
							{
								PeptideHit k = *tt;
								if (hit->getMetaValue("scoring") != t->getMetaValue("scoring") & tt->getMetaValue("scoring") == t->getMetaValue("scoring"))
								{
									//use SEQAN similarity scoring
									AASequence S1,S2;
									String SS1, SS2;
									const char *pc1, *pc2;
									DoubleReal c;
									S1 =tt->getSequence();
									SS1= S1.toUnmodifiedString(); 
									pc1 = SS1.c_str();
									S2 =hit->getSequence();
									SS2= S2.toUnmodifiedString();
									pc2 = SS2.c_str();
									typedef::seqan::String< ::seqan::AminoAcid > TSequence;
									TSequence seq1=pc1;
									TSequence seq2=pc2;
									::seqan::Score<int, ::seqan::Pam<> > pam(110, -10, -10);
									::seqan::Align<TSequence, ::seqan::ArrayGaps> align, self1, self2;
									::seqan::resize(rows(align), 2);
									::seqan::resize(rows(self1), 2);
									::seqan::resize(rows(self2), 2);
									::seqan::assignSource(row(align, 0), seq1);
									::seqan::assignSource(row(align, 1), seq2);
									::seqan::assignSource(row(self1, 0), seq1);
									::seqan::assignSource(row(self1, 1), seq1);
									::seqan::assignSource(row(self2, 0), seq2);
									::seqan::assignSource(row(self2, 1), seq2); 
	
									vector<DoubleReal> temp;
									temp.push_back(globalAlignment(self1, pam, ::seqan::NeedlemanWunsch()));
									temp.push_back(globalAlignment(self2, pam, ::seqan::NeedlemanWunsch()));
									DoubleReal b;
									c = (DoubleReal)globalAlignment(align, pam,::seqan::NeedlemanWunsch());
									b = *( min_element( temp.begin(), temp.end() ) );
									c=c/b;
									if(c<0){c=0;}
									if(c>a)									
									{										
										a=c;
										z=(double)tt->getMetaValue("PEP")*a;
									}	
								}	
							}
							a_score+=z;
							a_sim+=a;
							myset.insert(t->getMetaValue("scoring"));
						}
					}
					//the meta value similarity corresponds to the sum of the similarities. Note that if similarity equals the number of search engines, the 
					//same peptide has been assigned by all engines
					::std::cout <<hit->getSequence()<<" a_score="<<a_score<<" a_sim="<< a_sim <<::std::endl;
					vector<DoubleReal> ScoreSim;
					ScoreSim.push_back(a_score/a_sim);
					ScoreSim.push_back(a_sim);				
					scores.insert(make_pair(hit->getSequence(),ScoreSim));
					++hit_count;
			}
		}

		//Replace IDs by consensus
		ids.clear();
		ids.resize(1);
		ids[0].setScoreType(String("Consensus_PEP (") + score_type +")");
		ids[0].setHigherScoreBetter(FALSE);
		for (Map<AASequence,vector<DoubleReal> >::const_iterator it = scores.begin(); it != scores.end(); ++it)
		{
			PeptideHit hit;
			hit.setSequence(it->first);
			hit.setScore(it->second[0]);
			hit.setMetaValue("similarity", it->second[1]);
			ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
			cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

		}
	}

	
	
	
	
	
	void ConsensusID::majority_(vector<PeptideIdentification>& ids)
	{
		UInt min_number_of_engines = (UInt)param_.getValue("min_number_of_engines");
		UInt number_of_runs = (UInt)param_.getValue("number_of_runs");
		map<AASequence, vector<PeptideHit> > hits;

		for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
		{
			// collect the PeptideHits
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				PeptideHit hit = *pit;
				hit.metaRegistry().registerName(it->getScoreType(), "Score type");
				hit.setMetaValue(it->getScoreType(), pit->getScore());									
				hits[pit->getSequence()].push_back(hit);
			}
		}


		// now process the hits
		vector<PeptideHit> new_hits;
		for (map<AASequence, vector<PeptideHit> >::iterator pit = hits.begin(); pit != hits.end(); ++pit)
		{
			if (pit->second.size() >= min_number_of_engines)
			{
				if (pit->second.size() == number_of_runs)
				{
					double score(1);
					PeptideHit new_hit = *pit->second.begin();
					for (vector<PeptideHit>::const_iterator hit = pit->second.begin(); hit != pit->second.end(); ++hit)
					{
						score *= hit->getScore();
						vector<String> meta_keys;
            			hit->getKeys(meta_keys);
            			for (vector<String>::const_iterator kit = meta_keys.begin(); kit != meta_keys.end(); ++kit)
            			{
              				new_hit.setMetaValue(*kit, hit->getMetaValue(*kit));
            			}
					}
					new_hit.setScore(pow(score, 1.0/(double)pit->second.size()));
					new_hit.setProteinAccessions(vector<String>());
					new_hits.push_back(new_hit);
				}
				else
				{
					if (pit->second.size() < number_of_runs)
					{
						double score(1);
						PeptideHit new_hit = *pit->second.begin();
						for (vector<PeptideHit>::const_iterator hit = pit->second.begin(); hit != pit->second.end(); ++hit)
						{
							score *= hit->getScore();
							vector<String> meta_keys;
							hit->getKeys(meta_keys);
							for (vector<String>::const_iterator kit = meta_keys.begin(); kit != meta_keys.end(); ++kit)
							{
								new_hit.setMetaValue(*kit, hit->getMetaValue(*kit));
							}
						}
						new_hit.setScore(pow(score, 1.0/(double)pit->second.size())); // if 2 engines got the ID, 3 is number of runs -> score ** 3/2
						new_hit.setProteinAccessions(vector<String>());
						new_hits.push_back(new_hit);
					}
					else
					{
						cerr << "Not defined what happens if number_of_runs < peptide hit size of one peptide?!?!? " << endl;
					}
				}
			}
			else
			{
				cerr << "Peptide only identified by " << pit->second.size() << ", need to have a least " << min_number_of_engines << "!" << endl;
			}
		}

		PeptideIdentification id = *ids.begin();
		id.setHits(new_hits);
		id.assignRanks();
		ids.resize(1);
		ids[0] = id;		
		return;
	}



//THIS SHOULD ALREADY BE DONE IN /APPLICATION/TOPP/ConsensusID.C
	void ConsensusID::mapIdentifications_(vector<PeptideIdentification >& sorted_ids  , const vector<PeptideIdentification>& ids)
	{
		DoubleReal mz_delta = 0.01;
		DoubleReal rt_delta = 0.01;
		for (vector<PeptideIdentification>::const_iterator it1 = ids.begin(); it1 != ids.end(); ++it1)
		{
			DoubleReal rt1(it1->getMetaValue("RT"));
			DoubleReal mz1(it1->getMetaValue("MZ"));
			PeptideIdentification new_ids;
			for (vector<PeptideIdentification>::const_iterator it2 = it1 + 1; it2 != ids.end(); ++it2)
			{
				DoubleReal rt2(it2->getMetaValue("RT"));
				DoubleReal mz2(it2->getMetaValue("MZ"));
				if (fabs(rt1 - rt2) < rt_delta && fabs(mz1 - mz2) < mz_delta)
				{
					if (new_ids.empty() == TRUE)
					{
						new_ids=(*it1);
					}
					//inset all hits from *it2 into the existing Peptide Identification
					else
					{
						vector<PeptideHit> hits;
						for (vector<PeptideHit>::const_iterator pit = it2->getHits().begin(); pit != it2->getHits().end();++pit)
							{
								new_ids.insertHit(*pit);
							}
					}
				}
			}
			sorted_ids.push_back(new_ids);
		}
		return;
	}


} // namespace OpenMS
