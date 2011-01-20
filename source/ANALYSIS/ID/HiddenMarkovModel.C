// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>

#define SIMPLE_DEBUG2
#undef  SIMPLE_DEBUG2

#define HIDDEN_MARKOV_MODEL_DEBUG
#undef HIDDEN_MARKOV_MODEL_DEBUG

#define STATE_DEBUG
#undef STATE_DEBUG

#define EVALUATE_DEBUG
#undef EVALUATE_DEBUG

#ifdef EVALUATE_DEBUG
#include <gsl/gsl_statistics.h>
#endif

using namespace std;

namespace OpenMS 
{
	// HMMState implementation
	HMMState::HMMState()
		: hidden_(true)
	{
	}

	HMMState::HMMState(const String& name, bool hidden)
		: hidden_(hidden),
			name_(name)
	{
		#ifdef STATE_DEBUG
		cerr << "State: " << name << ", hidden=" << hidden << endl;
		#endif
	}

	HMMState::HMMState(const HMMState& state)
		: hidden_(state.hidden_),
			name_(state.name_)
	{
	}
	/*		pre_states_(state.pre_states_),
			succ_states_(state.succ_states_)
	{
		
	}*/
	
	HMMState::~HMMState()
	{
	}
	
	void HMMState::setName(const String& name)
	{
		name_ = name;
	}

	const String& HMMState::getName() const
	{
		return name_;
	}

	HMMState& HMMState::operator = (const HMMState& state)
	{
		hidden_ = state.hidden_;
		name_ = state.name_;
		pre_states_.clear();
		succ_states_.clear();
		//pre_states_ = state.pre_states_;
		//succ_states_ = state.succ_states_;
		return *this;
	}

	void HMMState::addSuccessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
		cerr << "'" << state->getName() << "'" << endl;
#endif
		succ_states_.insert(state);
	}

	void HMMState::deleteSuccessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
		succ_states_.erase(state);
	}

	void HMMState::addPredecessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
		pre_states_.insert(state);
	}

	void HMMState::deletePredecessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
		pre_states_.erase(state);
	}

	const set<HMMState*>& HMMState::getPredecessorStates() const
	{
		return pre_states_;
	}

	const set<HMMState*>& HMMState::getSuccessorStates() const
	{
		return succ_states_;
	}

	bool HMMState::isHidden() const
	{
		return hidden_;
	}

	void HMMState::setHidden(bool hidden)
	{
		hidden_ = hidden;
	}

	// The hidden markov model
	HiddenMarkovModel::HiddenMarkovModel() 
		: pseudo_counts_(0.0)
	{
	}

	HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel& hmm)
	{
		copy_(hmm);
	}
	
	HiddenMarkovModel::~HiddenMarkovModel()
	{
		clear();
	}

	HiddenMarkovModel& HiddenMarkovModel::operator = (const HiddenMarkovModel& hmm)
	{
		if (&hmm != this)
		{
			clear();
			copy_(hmm);
		}
		return *this;
	}

	HMMState* HiddenMarkovModel::getState(const String& name)
	{
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_[name];
	}

	const HMMState* HiddenMarkovModel::getState(const String& name) const
	{
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_.find(name)->second;
	}
	
	void HiddenMarkovModel::writeGraphMLFile(const String& filename)
	{
		set<HMMState*> states = states_;
		Map<HMMState*, vector<HMMState*> > transitions;
		
		ofstream out(filename.c_str());
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;

		out << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns/graphml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" " <<
         "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/graphml http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd\" " <<
				 "xmlns:y=\"http://www.yworks.com/xml/graphml\">" << endl;

		out << "<key id=\"d0\" for=\"node\" yfiles.type=\"nodegraphics\"/>" << endl;
  	out << "<key id=\"d1\" for=\"edge\" yfiles.type=\"edgegraphics\"/>" << endl;
    out << "  <graph id=\"G\" edgedefault=\"directed\">" << endl;
		for (set<HMMState*>::const_iterator it = states.begin(); it != states.end(); ++it)
		{
			out << "    <node id=\"" << (*it)->getName() << "\">" << endl;
			out << "      <data key=\"d0\">" << endl;
			out << "        <y:ShapeNode>" << endl;
			out << "          <y:NodeLabel>" << (*it)->getName() << "</y:NodeLabel>" << endl;
			out << "        </y:ShapeNode>" << endl;
			out << "      </data>" << endl;
			out << "    </node>" << endl;

			set<HMMState*> succ = (*it)->getSuccessorStates();
			for (set<HMMState*>::const_iterator sit = succ.begin(); sit != succ.end(); ++sit)
			{
				transitions[*it].push_back(*sit);
			}
		}


		
		for (Map<HMMState*, vector<HMMState*> >::const_iterator it = transitions.begin(); it != transitions.end(); ++it)
		{
			for (vector<HMMState*>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				out << "    <edge source=\"" << it->first->getName() << "\" target=\"" << (*it1)->getName() << "\" directed=\"true\">" << endl;
				out << "      <data key=\"d1\">" << endl;
				out << "        <y:PolyLineEdge>" << endl;
				out << "          <y:EdgeLabel>" << getTransitionProbability_(it->first, *it1) << "</y:EdgeLabel>" << endl;
				out << "        </y:PolyLineEdge>" << endl;
				out << "      </data>" << endl;
				out << "    </edge>" << endl;
			}
		}

		out << "  </graph>" << endl;
		out << "</graphml>" << endl;
		
	}

	void HiddenMarkovModel::write(ostream& out) const
	{
		// write states
		for (set<HMMState*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			out << "State " << (*it)->getName();

			if (!(*it)->isHidden())
			{
				out << " false";
			}
			out << endl;
		}

		// write transitions
		for (Map<HMMState*, Map<HMMState*, DoubleReal> >::const_iterator it1 = trans_.begin(); it1 != trans_.end(); ++it1)
		{
			for (Map<HMMState*, DoubleReal>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				out << "Transition " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second << endl;
			}
		}
	
		// write synonyms
		/*
		for (Map<String, Map<String, pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (Map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				out << "Synonym " << it->first << " " << it2->first << " " << it2->second.first << " " << it2->second.second << endl;
			}
		}*/

    for (Map<HMMState*, Map<HMMState*, std::pair<HMMState*, HMMState*> > >::const_iterator it1 = synonym_trans_.begin(); it1 != synonym_trans_.end(); ++it1)
    {
      for (Map<HMMState*, std::pair<HMMState*, HMMState*> >::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
				//cerr << "Synonym " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second.first->getName() << " " << it2->second.second->getName() << endl;
        out << "Synonym " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second.first->getName() << " " << it2->second.second->getName() << endl;
      }
    }

		return;
	}

	void HiddenMarkovModel::clear()
	{
		for (set<HMMState*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
    {
      delete *it;
    }
    trans_.clear();
    count_trans_.clear();
		train_count_trans_all_.clear();
    forward_.clear();
    backward_.clear();
    name_to_state_.clear();
    train_emission_prob_.clear();
    init_prob_.clear();
    states_.clear();
    trained_trans_.clear();
    synonym_trans_.clear();
    synonym_trans_names_.clear();	
		var_modifications_ = StringList::create("");
		return;
	}
	
	Size HiddenMarkovModel::getNumberOfStates() const
	{
		return states_.size();
	}

	void HiddenMarkovModel::addNewState(HMMState* s)
	{
		states_.insert(s);
		if (name_to_state_.find(s->getName()) == name_to_state_.end())
		{
			name_to_state_[s->getName()] = s;
		}
		else
		{
			cerr << "HiddenMarkovModel: state name '" << s->getName() << "' (" << s << ") already used!" << endl;
		}
	}

	void HiddenMarkovModel::addNewState(const String& name)
	{
		HMMState* s = new HMMState(name);
		states_.insert(s);
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			name_to_state_[name] = s;
		}
		else
		{
			cerr << "HiddenMarkovModel: state name '" << name << "' (" << name_to_state_[name] << ") already used!" << endl;
		}
		return;
	}

	void HiddenMarkovModel::addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2)
	{
#ifdef SIMPLE_DEBUG2
		cerr << "addSynonymTransition: '" << name1 << "' -> '" << name2 << "' :: '" << synonym1 << "' -> '" << synonym2 << "'" << endl;
#endif
		if (name_to_state_.find(name1) == name_to_state_.end())
		{
			cerr << "state '" << name1 << "' unknown" << endl;
		}
    if (name_to_state_.find(name2) == name_to_state_.end())
    {
      cerr << "state '" << name2 << "' unknown" << endl;
    }
    if (name_to_state_.find(synonym1) == name_to_state_.end())
    {
      cerr << "state '" << synonym1 << "' unknown" << endl;
    }
    if (name_to_state_.find(synonym2) == name_to_state_.end())
    {
      cerr << "state '" << synonym2 << "' unknown" << endl;
    }
		synonym_trans_names_[synonym1][synonym2] = make_pair(name1, name2);

		synonym_trans_[name_to_state_[synonym1]][name_to_state_[synonym2]] = make_pair(name_to_state_[name1], name_to_state_[name2]);

/*
		for (Map<String, Map<String, std::pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
    {
      for (Map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      {
        synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] =
            make_pair(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
      }
    }
	*/	
	}

	void HiddenMarkovModel::setTransitionProbability_(HMMState * s1, HMMState * s2, DoubleReal trans_prob)
	{
		trans_[s1][s2] = trans_prob;
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
		training_steps_count_[s1][s2] = 0;
	}

	void HiddenMarkovModel::setTransitionProbability(const String& s1, const String& s2, DoubleReal trans_prob)
	{
    OPENMS_PRECONDITION(
		name_to_state_.find(s1) != name_to_state_.end() && name_to_state_.find(s2) != name_to_state_.end(),
		String("HiddenMarkovModel::setTransitionProbability(" + String(s1) + ", " + String(s2) + ", " +  String(trans_prob) + "), no suchstate!").c_str());
#ifdef SIMPLE_DEBUG2
    cerr << "setTransitionProbability: '" << s1 << "' -> '" << s2 << "'" << " " << trans_prob << endl;
#endif

		trans_[name_to_state_[s1]][name_to_state_[s2]] = trans_prob;
		name_to_state_[s1]->addSuccessorState(name_to_state_[s2]);
		name_to_state_[s2]->addPredecessorState(name_to_state_[s1]);
		enabled_trans_[name_to_state_[s1]].insert(name_to_state_[s2]);
		training_steps_count_[name_to_state_[s1]][name_to_state_[s2]] = 0;
	}

	DoubleReal HiddenMarkovModel::getTransitionProbability(const String& s1, const String& s2) const
	{
		if (name_to_state_.find(s1) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, s1);
		}
		HMMState* state1 = name_to_state_[s1];
		if (name_to_state_.find(s2) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, s2);
		}
		HMMState* state2 = name_to_state_[s2];
		return getTransitionProbability_(state1, state2);
	}
	
	DoubleReal HiddenMarkovModel::getTransitionProbability_(HMMState* s1, HMMState* s2) const
	{
		HMMState* state1 = s1;
		HMMState* state2 = s2;

#ifdef SIMPLE_DEBUG2
		cerr << "getTransitionProbability(" << s1->getName() << ", " << s2->getName() << ")" << endl;
#endif
		
		// check if transition is a synonym transition
		if (synonym_trans_.find(state1) != synonym_trans_.end() && 
				synonym_trans_.find(state1)->second.find(state2) != synonym_trans_.find(state1)->second.end())
		{
			HMMState* tmp = synonym_trans_.find(state1)->second.find(state2)->second.first;
			state2 = synonym_trans_.find(state1)->second.find(state2)->second.second; //name_to_state_[p.second];
			state1 = tmp;
		}
		
#ifdef SIMPLE_DEBUG2
		cerr << "getTransitionProbability: " << state1->getName() << " " << state2->getName() << endl;
#endif
		
		// get the transition prob 
		if (trans_.find(state1) != trans_.end() && trans_.find(state1)->second.find(state2) != trans_.find(state1)->second.end())
		{
		/*
			// TODO !!!! ???? path length correction
			if (state2->getName().hasSubstring("next"))
			{
				return sqrt(trans_[state1][state2]);
			}
			else
			{*/
				return trans_.find(state1)->second.find(state2)->second;
				/*
			}
			*/
		}
		else
		{
			return 0;
		}
	}

	void HiddenMarkovModel::train()
	{
#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "HiddenMarkovModel::train()" << endl;
#endif
		trained_trans_.clear();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc forward part" << endl;
		#endif
		
		calculateForwardPart_();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc backward part" << endl;
		#endif
		
		calculateBackwardPart_();
	

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc px" << endl;
		#endif
	
		// calc p_x from forward part
		DoubleReal px(0);

		for (Map<HMMState*, DoubleReal>::const_iterator it1 = train_emission_prob_.begin(); it1 != train_emission_prob_.end(); ++it1)
		{
			for (set<HMMState*>::const_iterator it2 = it1->first->getPredecessorStates().begin(); it2 != it1->first->getPredecessorStates().end(); ++it2)
			{
				px += getForwardVariable_(*it2);
			}
		}

		DoubleReal num_px(0);
		if (px != 0)
		{
			num_px = 1.0/px;
		}

#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "px=" << num_px << endl;
		cerr << "add contributions to count_trans" << endl;
#endif

		// add contributions to count_trans_
		for (set<pair<HMMState*, HMMState*> >::const_iterator it = trained_trans_.begin(); it != trained_trans_.end(); ++it)
		{
			DoubleReal tmp(0);
			tmp = num_px * getForwardVariable_(it->first) * getBackwardVariable_(it->second) * getTransitionProbability_(it->first, it->second);
			tmp += pseudo_counts_;
			HMMState* s1 = it->first;
			HMMState* s2 = it->second;

			if (synonym_trans_.find(s1) != synonym_trans_.end() && synonym_trans_[s1].find(s2) != synonym_trans_[s1].end())
			{
				HMMState* tmp = synonym_trans_[s1][s2].first;
				s2 = synonym_trans_[s1][s2].second;
				s1 = tmp;
			}
			
			train_count_trans_all_[s1][s2].push_back(tmp);
			if (count_trans_.find(s1) != count_trans_.end() && count_trans_[s1].find(s2) != count_trans_[s1].end())
			{
				count_trans_[s1][s2] += tmp;
			}
			else
			{
				count_trans_[s1][s2] = tmp;
			}
			training_steps_count_[s1][s2]++;
		}
	}

	void HiddenMarkovModel::calculateForwardPart_()
	{
		forward_.clear();
		set<HMMState*> succ;
		for (Map<HMMState*, DoubleReal>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			//cerr << it->first << " " << it->second << endl;
			//cerr << it->first->getName() << endl;
			forward_[it->first] = it->second;
		}
		
		for (Map<HMMState*, DoubleReal>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());
		
			while (succ.size() != 0)
			{
				set<HMMState*> succ_new;
				for (set<HMMState*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
				{
					set<HMMState*> pre = (*it)->getPredecessorStates();
					DoubleReal sum(0);
					for (set<HMMState*>::const_iterator it2 = pre.begin(); it2 != pre.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getName() << " " << getForwardVariable_(*it2) << " * " << getTransitionProbability(*it2, *it) << endl;
						#endif
						sum += getForwardVariable_(*it2) * getTransitionProbability_(*it2, *it);
						trained_trans_.insert(make_pair(*it2, *it));
					}
					forward_[*it] = sum;
					succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "f: " << (*it)->getName() << "\t" << sum << endl;
					#endif
				}
				succ = succ_new;
			}
		}
	}

	void HiddenMarkovModel::calculateBackwardPart_()
	{
		backward_.clear();
		set<HMMState*> pre;
		for (Map<HMMState*, DoubleReal>::iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			backward_[it->first] = it->second;
		}

		for (Map<HMMState*, DoubleReal>::iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			#ifdef HIDDEN_MARKOV_MODEL_DEBUG
			cerr << "b:" << it->first << " " <<  it->first->getName() << " " << it->second << endl;
			#endif
			pre.insert(it->first->getPredecessorStates().begin(), it->first->getPredecessorStates().end());
		
			while (pre.size() != 0)
			{
				set <HMMState*> pre_new;
				for (set<HMMState*>::const_iterator it = pre.begin(); it != pre.end(); ++it)
				{
					set<HMMState*> succ = (*it)->getSuccessorStates();
					DoubleReal sum(0);
					for (set<HMMState*>::const_iterator it2 = succ.begin(); it2 != succ.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getName() << " " << getBackwardVariable_(*it2) << " * " << getTransitionProbability(*it, *it2) << endl;
						#endif
						sum += getBackwardVariable_(*it2) * getTransitionProbability_(*it, *it2);
						trained_trans_.insert(make_pair(*it, *it2));
					}
					backward_[*it] = sum;
					pre_new.insert((*it)->getPredecessorStates().begin(), (*it)->getPredecessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "b: " << (*it)->getName() << "\t" << sum << endl;
					#endif
				}
				pre = pre_new;
			}
		}
	}

	DoubleReal HiddenMarkovModel::getForwardVariable_(HMMState* state)
	{
		return forward_.find(state) != forward_.end() ? forward_[state] : 0;
	}

	DoubleReal HiddenMarkovModel::getBackwardVariable_(HMMState* state)
	{
		return backward_.find(state) != backward_.end() ? backward_[state] : 0;
	}

	void HiddenMarkovModel::evaluate()
	{
		for (Map<HMMState*, Map<HMMState*, DoubleReal> >::const_iterator it1 = count_trans_.begin(); it1 != count_trans_.end(); ++it1)
		{
#ifdef EVALUATE_DEBUG
			cerr <<  it1->first->getName() << endl;
#endif
			DoubleReal sum(0);
			for (Map<HMMState*, DoubleReal>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				if (count_trans_.find(it1->first) != count_trans_.end() && 
						count_trans_[it1->first].find(it2->first) != count_trans_[it1->first].end())
				{
					sum += count_trans_[it1->first][it2->first];
#ifdef EVALUATE_DEBUG
					cerr << it1->first->getName() << " " << it2->first->getName() << " ";
					
					//<< count_trans_[it1->first][it2->first] << endl;
					for (vector<DoubleReal>::const_iterator it = train_count_trans_all_[it1->first][it2->first].begin(); it != train_count_trans_all_[it1->first][it2->first].end(); ++it)
					{
						cerr << *it << " ";
					}
					vector<DoubleReal> data = train_count_trans_all_[it1->first][it2->first];
					std::sort(data.begin(),data.end());
		      DoubleReal mean = gsl_stats_mean(&data.front(),1,data.size());
		      DoubleReal variance = gsl_stats_variance_m(&data.front(),1,data.size(),mean);
					cerr << "mean=" << mean << ", variance=" << variance << endl;
#endif
				}
			}

			if (sum != 0)
			{
				for (Map<HMMState*, DoubleReal>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					if (count_trans_.find(it1->first) != count_trans_.end() && 
							count_trans_[it1->first].find(it2->first) != count_trans_[it1->first].end())
					{
						trans_[it1->first][it2->first] = count_trans_[it1->first][it2->first] / sum;
					}
				}
			}
		}
	}

	void HiddenMarkovModel::setInitialTransitionProbability(const String& state, DoubleReal prob)
	{
		OPENMS_PRECONDITION(name_to_state_.find(state) != name_to_state_.end(), String("HiddenMarkovModel::setInitialTransitionProbability(" + state + ", " + String(prob) + "), no suchstate!").c_str());
		//cerr << state << " " << prob << endl;
		init_prob_[name_to_state_[state]] = prob;
	}

	void HiddenMarkovModel::clearInitialTransitionProbabilities()
	{
		init_prob_.clear();
	}

	void HiddenMarkovModel::setTrainingEmissionProbability(const String& state, DoubleReal prob)
	{
#ifdef SIMPLE_DEBUG2
		cerr << "setTrainingEmissionProbability(" << state << "(" << name_to_state_[state] << "), " << prob << ")" << endl;
#endif

		OPENMS_PRECONDITION(name_to_state_.find(state) != name_to_state_.end(), String("HiddenMarkovModel::setTrainingEmissionProbability(" + state + ", " + String(prob) + "), no such state!").c_str());
		train_emission_prob_[name_to_state_[state]] = prob;
	}

	void HiddenMarkovModel::clearTrainingEmissionProbabilities()
	{
		train_emission_prob_.clear();
	}

	void HiddenMarkovModel::enableTransition(const String& s1, const String& s2)
	{
#ifdef SIMPLE_DEBUG2
		cerr << "enableTransition: '" << s1 << "' -> '" << s2 << "'" << endl;
#endif
		enableTransition_(name_to_state_[s1], name_to_state_[s2]);
	}

	void HiddenMarkovModel::enableTransition_(HMMState* s1, HMMState* s2)
	{
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
	}

	void HiddenMarkovModel::disableTransition(const String& s1, const String& s2)
	{
		disableTransition_(name_to_state_[s1], name_to_state_[s2]);
	}
	
	void HiddenMarkovModel::disableTransition_(HMMState* s1, HMMState* s2)
	{
		s1->deleteSuccessorState(s2);
		s2->deletePredecessorState(s1);
		enabled_trans_[s1].erase(s2);
	}

	void HiddenMarkovModel::disableTransitions()
	{
		for (Map<HMMState*, set<HMMState*> >::const_iterator it = enabled_trans_.begin(); it != enabled_trans_.end(); ++it)
		{
			for (set<HMMState*>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				it->first->deleteSuccessorState(*it2);
				(*it2)->deletePredecessorState(it->first);
			}
		}
		enabled_trans_.clear();
	}

	void HiddenMarkovModel::calculateEmissionProbabilities(Map<HMMState*, DoubleReal>& emission_probs)
	{
		Map<HMMState*, DoubleReal> states = init_prob_;

		while (states.size() != 0)
		{
			Map<HMMState*, DoubleReal> tmp = states;
			for (Map<HMMState*, DoubleReal>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
			{
				for (set<HMMState*>::const_iterator it2 = it->first->getSuccessorStates().begin(); it2 != it->first->getSuccessorStates().end(); ++it2)
				{
					//cerr << "->" << (*it2)->getName() << "=(from " << it->first->getName() << "): "; // << states[it->first] << " * " << getTransitionProbability_(it->first, *it2) << " ";
					if (states.find(*it2) != states.end())
					{
						//cerr << " += " << states[it->first] * getTransitionProbability_(it->first, *it2);
						states[*it2] += states[it->first] * getTransitionProbability_(it->first, *it2);
					}
					else
					{
						//cerr << " = " << states[it->first] * getTransitionProbability_(it->first, *it2);
						states[*it2] = states[it->first] * getTransitionProbability_(it->first, *it2);
					}
					if (!(*it2)->isHidden())
					{
						if (emission_probs.find(*it2) != emission_probs.end())
						{
							//cerr << " emission: += " << states[it->first] * getTransitionProbability_(it->first, *it2);
							emission_probs[*it2] += states[it->first] * getTransitionProbability_(it->first, *it2);
						}
						else
						{	
							//cerr << " emission: += " << states[it->first] * getTransitionProbability_(it->first, *it2);
							emission_probs[*it2] = states[it->first] * getTransitionProbability_(it->first, *it2);
						}
					}
					//cerr << endl;
				}
				states.erase(it->first);
			}
		}
		
		return;
	}

	void HiddenMarkovModel::dump()
	{
		cerr << "dump of transitions: " << endl;
		for (Map<HMMState*, Map<HMMState*, DoubleReal> >::const_iterator it = trans_.begin(); it != trans_.end(); ++it)
		{
			for (Map<HMMState*, DoubleReal>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cout << it->first->getName() << " -> " << it1->first->getName() << " " << it1->second << " " << training_steps_count_[it->first][it1->first] << ": ";
				vector<DoubleReal> all_trans = train_count_trans_all_[it->first][it1->first];
				
				if (all_trans.size() != 0)
				{
					DoubleReal sum = accumulate(all_trans.begin(), all_trans.end(), 0.0);
					DoubleReal avg(sum/DoubleReal(all_trans.size()));
					DoubleReal rsd(0);
					for (Size i = 0; i != all_trans.size(); ++i)
					{
						cout << all_trans[i] << " ";
						rsd += abs(all_trans[i] - avg);
					}
					cout << "rsd=" << rsd / DoubleReal(all_trans.size()) / avg;
					cout << ", avg=" << avg;
				}
				
				cout << endl;
			}
		}
		
		cerr << "dump completed" << endl;
	}

	void HiddenMarkovModel::forwardDump()
	{
   set<HMMState*> succ;
    for (Map<HMMState*, DoubleReal>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
    {
      succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());

      while (succ.size() != 0)
      {
        set<HMMState*> succ_new;
        for (set<HMMState*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
        {
					cerr << (*it)->getName() << endl;
          succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
        }
        succ = succ_new;
      }
    }
	}

	void HiddenMarkovModel::estimateUntrainedTransitions()
	{
		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "Transition training stats" << endl;
		for (Map<HMMState*, Map<HMMState*, Size> >::const_iterator it = training_steps_count_.begin(); it != training_steps_count_.end(); ++it)
		{
			for (Map<HMMState*, Size>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cerr << it->first->getName() << " " << it1->first->getName() << " " << it1->second << endl;
			}
		}

		cerr << "Estimation" << endl;
		#endif


		set<const Residue*> residues(ResidueDB::getInstance()->getResidues("Natural20"));
    for (StringList::const_iterator it = var_modifications_.begin(); it != var_modifications_.end(); ++it)
    {
      residues.insert(ResidueDB::getInstance()->getModifiedResidue(*it));
    }

		// pathways axyz and bxyz and the first two explicitely modeled ones
		HMMState* s2 = 0;
		HMMState* end_state = name_to_state_["end"];
		StringList pathways = StringList::create("axyz,axyz1,axyz1,bxyz,bxyz1,bxyz2");	
		for (StringList::const_iterator pathway_it = pathways.begin(); pathway_it != pathways.end(); ++pathway_it)
		{
			String pathway = *pathway_it;
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				s2 = name_to_state_[pathway];
				for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
				{
					AASequence first_aa, second_aa;
					first_aa += *it;
					second_aa += *jt;
					String aa1(first_aa.toString()), aa2(second_aa.toString());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "Estimating: " << aa1 << " -> " << aa2 << " (" << training_steps_count_[name_to_state_[aa1 + aa2 + "_" + pathway]][s2] << ") :: " ;
					#endif

					if (training_steps_count_[name_to_state_[aa1 + aa2 + "_" + pathway]][s2] == 0)
					{
						Size count(0);
						DoubleReal sum(0);
						// "rows" of the amino acid matrix
						for (set<const Residue*>::const_iterator kt = residues.begin(); kt != residues.end(); ++kt)
						{
							AASequence third_aa;
							third_aa += *kt;
							String aa3(third_aa.toString());
							if (training_steps_count_[name_to_state_[aa1 + aa3 + "_" + pathway]][s2] != 0)
							{
								sum += trans_[name_to_state_[aa1 + aa3 + "_" + pathway]][s2];
								count++;
							}
						}
						// "columns" of the amino acid matrix
						for (set<const Residue*>::const_iterator kt = residues.begin(); kt != residues.end(); ++kt)
						{
							AASequence third_aa;
							third_aa += *kt;
							String aa3(third_aa.toString());

							if (training_steps_count_[name_to_state_[aa3 + aa2 + "_" + pathway]][s2] != 0)
							{
								sum += trans_[name_to_state_[aa3 + aa2 + "_" + pathway]][s2];
								count++;
							}
						}

						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << trans_[name_to_state_[aa1 + aa2 + "_" + pathway]][s2] << " to ";
						#endif
						if (count != 0)
						{
							#ifdef HIDDEN_MARKOV_MODEL_DEBUG
							cerr << "setting transitions of " << aa1 << aa2 << "_" << pathway << " -> " << pathway << " to " << sum/DoubleReal(count) << endl;
							#endif
							trans_[name_to_state_[aa1 + aa2 + "_" + pathway]][s2] = sum/DoubleReal(count);
							trans_[name_to_state_[aa1 + aa2 + "_" + pathway]][end_state] = 1 - sum/DoubleReal(count);
						}
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << sum/DoubleReal(count) << endl;
					}
					else
					{
						cerr << "not needed" << endl;
						#endif
					}
				}
			}
		}

		// sc and cr
		String sc_residues("HKRDE");
		for (String::ConstIterator it = sc_residues.begin(); it != sc_residues.end(); ++it)
		{
			String sc_res = *it;

			for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
			{
        AASequence second_aa;
				second_aa += *jt;
				String aa2(second_aa.toString());
			
				s2 = name_to_state_[sc_res];
				if (training_steps_count_[name_to_state_[aa2 + "_" + sc_res]][s2] == 0)
				{
					Size count(0);
					DoubleReal sum(0);
					for (set<const Residue*>::const_iterator kt = residues.begin(); kt != residues.end(); ++kt)
					{
						AASequence third_aa;
						third_aa += *kt;
						String aa3(third_aa.toString());

						HMMState* s1 = name_to_state_[aa3 + "_" + sc_res];
						if (training_steps_count_[s1][s2] != 0)
						{
							sum += trans_[s1][s2];
							count++;
						}
					}

					if (count != 0)
					{
						trans_[name_to_state_[aa2 + "_" + sc_res]][s2] = sum/DoubleReal(count);
						trans_[name_to_state_[aa2 + "_" + sc_res]][end_state] = 1 - sum/DoubleReal(count);
					}
				}
			}
		}

		StringList bk_pathways = StringList::create("bk-1,bk-2");

		for (StringList::const_iterator pathway_it = bk_pathways.begin(); pathway_it != bk_pathways.end(); ++pathway_it)
		{
			String pathway = *pathway_it;
			s2 = name_to_state_[pathway];
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				AASequence first_aa;
				first_aa += *it;
				String aa1(first_aa.toString());
				if (training_steps_count_[name_to_state_[aa1 + "_" + pathway]][s2] == 0)
				{
					Size count(0);
					DoubleReal sum(0);
					for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
					{
						AASequence second_aa;
						second_aa += *jt;
						String aa2(second_aa.toString());
						HMMState* s1 = name_to_state_[aa2 + "_" + pathway];
						if (training_steps_count_[s1][s2] != 0)
						{
							sum += trans_[s1][s2];
							count++;
						}
						//cerr << "Estimating transition of '" << aa1 << pathway << "' -> '" << pathway << "' to " << sum/(DoubleReal)count << endl;
						if (count != 0)
						{
							trans_[name_to_state_[aa1 + "_" + pathway]][s2] = sum/(DoubleReal)count;
							trans_[name_to_state_[aa1 + "_" + pathway]][end_state] = 1 - sum/(DoubleReal)count;
						}
					}
				}
			}
		}

	}

	void HiddenMarkovModel::setPseudoCounts(DoubleReal pseudo_counts)
	{
		pseudo_counts_ = pseudo_counts;
		return;
	}

	DoubleReal HiddenMarkovModel::getPseudoCounts() const
	{
		return pseudo_counts_;
	}

	void HiddenMarkovModel::setVariableModifications(const StringList& modifications)
	{
		var_modifications_ = modifications;
	}

	void HiddenMarkovModel::copy_(const HiddenMarkovModel& source)
	{
		Map<HMMState*, HMMState*> old_to_new;
    for (set<HMMState*>::const_iterator it = source.states_.begin(); it != source.states_.end(); ++it)
    {
      HMMState* s = new HMMState(**it);
      states_.insert(s);
      name_to_state_[s->getName()] = s;
      old_to_new[*it] = s;
    }
		
		// trans_
		for (Map<HMMState*, Map<HMMState*, DoubleReal> >::const_iterator it1 = source.trans_.begin(); it1 != source.trans_.end(); ++it1)
		{
			for (Map<HMMState*, DoubleReal>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
			}
		}

		// count_trans_
    for (Map<HMMState*, Map<HMMState*, DoubleReal> >::const_iterator it1 = source.count_trans_.begin(); it1 != source.count_trans_.end(); ++it1)
    {
      for (Map<HMMState*, DoubleReal>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        count_trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

		for (Map<HMMState*, Map<HMMState*, std::vector<DoubleReal> > >::const_iterator it1 = source.train_count_trans_all_.begin(); it1 != source.train_count_trans_all_.end(); ++it1)
		{
			for (Map<HMMState*, vector<DoubleReal> >::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        train_count_trans_all_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }

		}


		for (Map<HMMState*, Map<HMMState*, Size> >::const_iterator it1 = source.training_steps_count_.begin(); it1 != source.training_steps_count_.end(); ++it1)
    {
      for (Map<HMMState*, Size>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        training_steps_count_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

		// forward and backward are just temporary objects

		for (Map<HMMState*, DoubleReal>::const_iterator it = source.train_emission_prob_.begin(); it != source.train_emission_prob_.end(); ++it)
		{
			train_emission_prob_[old_to_new[it->first]] = it->second;
		}

		for (Map<HMMState*, DoubleReal>::const_iterator it = source.init_prob_.begin(); it != source.init_prob_.end(); ++it)
		{
			init_prob_[old_to_new[it->first]] = it->second;
		}

		for (std::set<std::pair<HMMState*, HMMState*> >::const_iterator it = source.trained_trans_.begin(); it != source.trained_trans_.end(); ++it)
		{
			trained_trans_.insert(make_pair(old_to_new[it->first], old_to_new[it->second]));
		}

		synonym_trans_names_ = source.synonym_trans_names_;
		pseudo_counts_ = source.pseudo_counts_;
		var_modifications_ = source.var_modifications_;


    for (Map<String, Map<String, std::pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
    {
      for (Map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      {
        synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] =
            make_pair(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
      }
    }


		for (Map<HMMState*, set<HMMState*> >::const_iterator it1 = source.enabled_trans_.begin(); it1 != source.enabled_trans_.end(); ++it1)
		{
			for (set<HMMState*>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				enabled_trans_[old_to_new[it1->first]].insert(old_to_new[*it2]);
			}
		}

	}
}


